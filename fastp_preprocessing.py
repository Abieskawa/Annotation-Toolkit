#!/usr/bin/env python3
"""
run_fastp.py

Process paired-end FASTQ files with fastp (RNA-only pipeline).

- Detects pairs named either <base>_1 / <base>_2 or <base>_R1 / <base>_R2
  with extensions .fq/.fastq[.gz]
- Optional *front* trimming ONLY:
    * one global integer (e.g., --trim-front 10)
    * or a 2-column file mapping sample_key -> value, where sample_key ends
      with _1/_2 or _R1/_R2, e.g.:
         SRR26245751_1  8
         SRR26245751_2  6
         4202P-HPGonad_R1 12
         4202P-HPGonad_R2 12
- Writes cleaned reads to --outdir as <base>_<R1or1>.cleaned.fastq.gz
- Generates per-sample HTML/JSON reports

Notes:
- For PE data, adapter trimming is robust via overlap; you can add
  --detect-adapter-pe for ultra-cleaning.
- --trim_poly_x is configurable by the user (default ON).
"""

import sys
import argparse
from pathlib import Path
from typing import Dict, Tuple, List, Optional, Set
import subprocess
import re

# ----------------------------
# Pair discovery (repo import preferred, local fallback otherwise)
# ----------------------------

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from Utils.utils import discover_pairs_by_patterns  # type: ignore
except Exception:
    _FASTQ_R1_RE = re.compile(
        r"^(?P<pre>.+?)(?P<tag>_R1|_1)(?P<post>(?:_.+)?)(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$",
        re.IGNORECASE,
    )

    def _discover_pairs_paths(
        wd: Path,
        r1_patterns: List[str],
        require_r2: bool = True,
        sort_paths: bool = False,
    ) -> List[Tuple[str, str, Path, Optional[Path]]]:
        r1_paths: List[Path] = []
        seen: Set[Path] = set()

        for pat in r1_patterns:
            for p in wd.glob(pat):
                if p.is_file() and p not in seen:
                    seen.add(p)
                    r1_paths.append(p)

        if sort_paths:
            r1_paths.sort(key=lambda x: x.name)

        pairs: List[Tuple[str, str, Path, Optional[Path]]] = []
        for r1p in r1_paths:
            m = _FASTQ_R1_RE.match(r1p.name)
            if not m:
                continue

            pre = m.group("pre")
            tag = m.group("tag")
            post = m.group("post")
            ext = m.group("ext")

            style = "R" if tag.upper().startswith("_R") else "1"
            mate_tag = tag[:-1] + "2"
            r2_name = f"{pre}{mate_tag}{post}{ext}"
            r2p = wd / r2_name

            base = f"{pre}{post}"

            if r2p.exists():
                pairs.append((base, style, r1p, r2p))
            elif not require_r2:
                pairs.append((base, style, r1p, None))

        return pairs

    def discover_pairs_by_patterns(wd: Path, r1_patterns: List[str]) -> List[Tuple[str, str, str, str]]:
        pairs = _discover_pairs_paths(wd, r1_patterns, require_r2=True, sort_paths=False)
        return [
            (base, style, r1p.name, r2p.name)
            for base, style, r1p, r2p in pairs
            if r2p is not None
        ]


# ----------------------------
# fastp runner (RNA-only)
# ----------------------------

R1_PATTERNS = [
    "*_1.fastq", "*_1.fq", "*_1.fastq.gz", "*_1.fq.gz",
    "*_R1.fastq", "*_R1.fq", "*_R1.fastq.gz", "*_R1.fq.gz",
    "*_1_*.fastq", "*_1_*.fq", "*_1_*.fastq.gz", "*_1_*.fq.gz",
    "*_R1_*.fastq", "*_R1_*.fq", "*_R1_*.fastq.gz", "*_R1_*.fq.gz",
]

FASTP_MAX_THREADS = 16


def load_trim_file(trim_file: str) -> Tuple[Dict[str, int], Dict[str, int]]:
    r1: Dict[str, int] = {}
    r2: Dict[str, int] = {}
    with open(trim_file) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                print(f"[warn] Bad trim line (skip): {line}", file=sys.stderr)
                continue
            key, val = parts
            try:
                iv = int(val)
            except ValueError:
                print(f"[warn] Non-integer trim value (skip): {line}", file=sys.stderr)
                continue
            if key.endswith(("_1", "_R1")):
                r1[key] = iv
            elif key.endswith(("_2", "_R2")):
                r2[key] = iv
            else:
                print(f"[warn] Key has no _1/_2 or _R1/_R2 suffix (skip): {key}", file=sys.stderr)
    return r1, r2


def discover_pairs(wd: Path) -> List[Tuple[str, str, str, str]]:
    return discover_pairs_by_patterns(wd, R1_PATTERNS)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Batch run fastp on paired-end FASTQs (RNA-only; front-only trimming).",
        formatter_class=argparse.RawTextHelpFormatter
    )
    p.add_argument("-d", "--wd", required=True, help="Working directory containing input FASTQs")
    p.add_argument("-o", "--outdir", required=True, help="Output directory (will be created if missing)")
    p.add_argument(
        "-w", "--threads",
        type=int, default=16,
        help="Threads for fastp (default: 16, max: 64)"
    )
    p.add_argument("-l", "--min-length", type=int, default=30, help="--length_required (default: 30)")
    p.add_argument("-q", "--qualified-phred", type=int, default=20, help="Qualified base threshold Q (default: 20)")
    p.add_argument("-e", "--average-qual", type=int, default=0, help="Drop reads with avg qual < this (0=off)")
    p.add_argument(
        "--trim-front",
        help="Global integer (front trim) OR path to 2-col file mapping sample_key -> value",
        default=""
    )
    p.add_argument(
        "--detect-adapter-pe",
        action="store_true",
        help="Enable adapter detection for PE (adds --detect_adapter_for_pe)"
    )
    p.add_argument("--report-dir", default="", help="Directory for HTML/JSON reports (default: outdir)")

    # Allow user to control polyX trimming; default ON for RNA pipeline.
    g = p.add_mutually_exclusive_group()
    g.add_argument(
        "--trim_poly_x", "--trim-poly-x",
        dest="trim_poly_x",
        action="store_true",
        help="Enable fastp --trim_poly_x (default: ON)"
    )
    g.add_argument(
        "--no-trim-poly-x",
        dest="trim_poly_x",
        action="store_false",
        help="Disable fastp --trim_poly_x"
    )
    p.set_defaults(trim_poly_x=True)

    return p


def run_fastp(
    wd: Path,
    outdir: Path,
    threads: int,
    min_length: int,
    qualified_phred: int,
    average_qual: int,
    trim_front: str,
    detect_adapter_pe: bool,
    trim_poly_x: bool,
    report_dir: Optional[Path]
) -> None:

    if not wd.exists():
        raise FileNotFoundError(f"Working dir not found: {wd}")

    outdir.mkdir(parents=True, exist_ok=True)

    if report_dir and not report_dir.exists():
        raise FileNotFoundError(f"Report dir not found: {report_dir}")

    if threads < 1:
        print(f"[warn] fastp threads must be >= 1; using 1 (requested {threads})", file=sys.stderr)
        threads = 1
    elif threads > FASTP_MAX_THREADS:
        raise ValueError(f"--threads exceeds max supported ({FASTP_MAX_THREADS}); got {threads}")

    # parse front trimming
    per_r1: Dict[str, int] = {}
    per_r2: Dict[str, int] = {}
    global_front: Optional[int] = None
    if trim_front:
        tf = Path(trim_front)
        if tf.exists():
            per_r1, per_r2 = load_trim_file(str(tf))
        else:
            try:
                global_front = int(trim_front)
            except ValueError:
                print(f"[warn] --trim-front must be int or file. Ignored: {trim_front}", file=sys.stderr)

    pairs = discover_pairs(wd)
    if not pairs:
        print(f"[info] No paired FASTQs found in {wd}", file=sys.stderr)
        return

    for base, style, r1, r2 in pairs:
        rtag1, rtag2 = ("R1", "R2") if style == "R" else ("1", "2")

        out1 = outdir / f"{base}_{rtag1}.cleaned.fastq.gz"
        out2 = outdir / f"{base}_{rtag2}.cleaned.fastq.gz"

        key_r1 = f"{base}_{rtag1}"
        key_r2 = f"{base}_{rtag2}"
        alt_r1 = f"{base}_{'1' if style == 'R' else 'R1'}"
        alt_r2 = f"{base}_{'2' if style == 'R' else 'R2'}"

        trim_r1 = per_r1.get(key_r1, per_r1.get(alt_r1, global_front))
        trim_r2 = per_r2.get(key_r2, per_r2.get(alt_r2, global_front))

        repdir = report_dir if report_dir else outdir
        html = repdir / f"{base}.fastp.html"
        jsn = repdir / f"{base}.fastp.json"

        cmd = [
            "fastp",
            "-w", str(threads),
            "-i", str(wd / r1),
            "-I", str(wd / r2),
            "-o", str(out1),
            "-O", str(out2),
            "-l", str(min_length),
            "-q", str(qualified_phred),
            "-h", str(html),
            "-j", str(jsn),
        ]

        if trim_poly_x:
            cmd += ["--trim_poly_x"]

        if average_qual and average_qual > 0:
            cmd += ["-e", str(average_qual)]

        if detect_adapter_pe:
            cmd += ["--detect_adapter_for_pe"]

        if isinstance(trim_r1, int):
            cmd += ["-f", str(trim_r1)]
        if isinstance(trim_r2, int):
            cmd += ["-F", str(trim_r2)]

        print("[info] Running:", " ".join(cmd), flush=True)
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[error] fastp failed for {base}: {e}", file=sys.stderr)
            raise

    print("[info] All samples processed.", flush=True)


def main() -> None:
    ap = build_argparser()
    a = ap.parse_args()

    run_fastp(
        wd=Path(a.wd).resolve(),
        outdir=Path(a.outdir).resolve(),
        threads=a.threads,
        min_length=a.min_length,
        qualified_phred=a.qualified_phred,
        average_qual=a.average_qual,
        trim_front=a.trim_front,
        detect_adapter_pe=a.detect_adapter_pe,
        trim_poly_x=a.trim_poly_x,
        report_dir=Path(a.report_dir).resolve() if a.report_dir else None
    )


if __name__ == "__main__":
    main()
