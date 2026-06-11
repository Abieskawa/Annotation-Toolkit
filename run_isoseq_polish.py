import argparse
import subprocess
from pathlib import Path


def run(cmd, log):
    print("[CMD]", " ".join(map(str, cmd)), flush=True)
    with open(log, "w") as f:
        subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)


def find_fastq(raw, name):
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        f = raw / f"{name}{ext}"
        if f.is_file() and f.stat().st_size > 0:
            return f
    raise FileNotFoundError(f"Missing FASTQ for {name} in {raw}")


def illumina_csv(raw, illumina_ids):
    files = []
    for srr in illumina_ids.split(","):
        srr = srr.strip()
        files.append(find_fastq(raw, f"{srr}_1"))
        files.append(find_fastq(raw, f"{srr}_2"))
    return ",".join(map(str, files))


def process_sample(pb, label, illumina_ids, args):
    raw = Path(args.raw)
    outdir = Path(args.outdir)

    flnc = outdir / "01_flnc_RNA_PacBio"
    cluster = outdir / "02_clustered_RNA_PacBio"
    polish = outdir / "03_polished_RNA_PacBio"
    tmp = outdir / "tmp"
    logs = outdir / "logs"

    for d in [flnc, cluster, polish, tmp, logs]:
        d.mkdir(parents=True, exist_ok=True)

    pb_fq = find_fastq(raw, pb)
    sr_fq = illumina_csv(raw, illumina_ids)

    trimmed_fq = tmp / f"{pb}.trimmed.fastq.gz"
    ubam = flnc / f"{pb}.flnc.bam"
    cbam = cluster / f"{pb}.clustered.bam"
    hq_gz = cluster / f"{pb}.clustered.hq.fasta.gz"
    hq = cluster / f"{pb}.clustered.hq.fasta"
    final = polish / f"{pb}.polished.hq.fasta"

    print(f"\n[{pb}] {label}", flush=True)
    print(f"PacBio FASTQ : {pb_fq}", flush=True)
    print(f"Final FASTA  : {final}", flush=True)

    if trimmed_fq.exists():
        trimmed_fq.unlink()

    cutadapt_cmd = [
        "cutadapt",
        "-g", "X" + args.primer5,
        "-a", args.primer3 + "X",
        "--poly-a",
        "--trimmed-only",
        "--minimum-length", str(args.min_len),
        "--cores", str(args.threads),
        "-o", str(trimmed_fq),
        str(pb_fq),
    ]

    cutadapt_log = logs / f"{pb}.cutadapt.log"
    try:
        run(cutadapt_cmd, cutadapt_log)

        if ubam.exists():
            ubam.unlink()

        run([
            "picard", "FastqToSam",
            f"F1={trimmed_fq}",
            f"O={ubam}",
            f"SM={pb}",
            "PL=PACBIO",
            "VALIDATION_STRINGENCY=SILENT",
        ], logs / f"{pb}.picard.log")

    finally:
        if trimmed_fq.exists():
            trimmed_fq.unlink()

    run([
        "isoseq3", "cluster",
        str(ubam),
        str(cbam),
        "--use-qvs",
        "--verbose",
        "-j", str(args.threads),
    ], logs / f"{pb}.isoseq3.log")

    if hq_gz.is_file() and hq_gz.stat().st_size > 0:
        run(["gunzip", "-f", str(hq_gz)], logs / f"{pb}.gunzip.log")
    if not hq.is_file() or hq.stat().st_size == 0:
        raise FileNotFoundError(f"IsoSeq3 HQ FASTA not found: {hq} (nor {hq_gz})")

    run([
        "lordec-correct",
        "-i", str(hq),
        "-2", sr_fq,
        "-o", str(final),
        "-T", str(args.threads),
        "-k", str(args.kmer),
        "-s", str(args.solid),
    ], logs / f"{pb}.lordec.log")

    if ubam.exists():
        ubam.unlink()

    print(f"[DONE] {final}", flush=True)


def main():
    p = argparse.ArgumentParser(
        description="Iso-Seq FASTQ preprocessing: 5p/3p primer trimming, IsoSeq3 clustering, LoRDEC polishing."
    )

    p.add_argument("-s", "--samples", required=True)
    p.add_argument("-r", "--raw", required=True)
    p.add_argument("-o", "--outdir", required=True)

    p.add_argument("-t", "--threads", type=int, required=True)
    p.add_argument("--primer5", required=True)
    p.add_argument("--primer3", required=True)
    p.add_argument("--min-len", type=int, required=True)
    p.add_argument("-k", "--kmer", type=int, required=True)
    p.add_argument("--solid", type=int, required=True)

    args = p.parse_args()

    with open(args.samples) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            pb, label, illumina = line.rstrip("\n").split("\t")
            process_sample(pb, label, illumina, args)

    print("\nFinal polished FASTA files:", flush=True)
    for f in sorted((Path(args.outdir) / "03_polished_RNA_PacBio").glob("*.polished.hq.fasta")):
        print(f, flush=True)


if __name__ == "__main__":
    main()
