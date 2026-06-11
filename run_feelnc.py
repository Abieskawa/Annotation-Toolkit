#!/usr/bin/env python3

import argparse
import gzip
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path


def check_tool(tool):
    """
    Check whether an external command exists in PATH.
    Example: STAR, stringtie, FEELnc_filter.pl.
    """
    if shutil.which(tool) is None:
        sys.exit(f"ERROR: required tool not found in PATH: {tool}")


def run(cmd, log_file=None, stdout_file=None, stderr_file=None):
    """
    Run one external command.

    cmd is a list, not one shell string.
    Example:
      ["STAR", "--runThreadN", "128", "--genomeDir", "index"]

    This avoids shell quoting problems with long paths.
    """
    cmd = [str(x) for x in cmd]

    print("RUN:", " ".join(cmd), flush=True)

    if log_file:
        with open(log_file, "a") as log:
            log.write(" ".join(cmd) + "\n")

    if stdout_file and stderr_file:
        with open(stdout_file, "w") as out, open(stderr_file, "w") as err:
            subprocess.run(cmd, check=True, stdout=out, stderr=err)

    elif stdout_file:
        with open(stdout_file, "w") as out:
            subprocess.run(cmd, check=True, stdout=out)

    elif stderr_file:
        with open(stderr_file, "w") as err:
            subprocess.run(cmd, check=True, stderr=err)

    else:
        subprocess.run(cmd, check=True)


def compute_genomeSAindex_nbases(genome_path):
    """
    Compute STAR --genomeSAindexNbases from total genome length.

    This follows the logic from your previous STAR mapping script:
      floor(log2(total_genome_length) / 2) - 1
    then cap at 14.
    """
    total_length = 0

    with open(genome_path, "r") as fh:
        for line in fh:
            if not line.startswith(">"):
                total_length += len(line.strip())

    computed = int(math.log2(total_length) / 2) - 1
    return min(14, computed)


def copy_or_decompress(infile, outfile):
    """
    Copy a plain file or decompress a .gz file.

    We use this to convert braker_utr.gtf.gz into a plain GTF for STAR,
    StringTie, and FEELnc.
    """
    infile = Path(infile)
    outfile = Path(outfile)

    if str(infile).endswith(".gz"):
        with gzip.open(infile, "rt") as fin, open(outfile, "w") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(infile, outfile)


def discover_paired_fastqs(fastq_dir):
    """
    Find paired-end cleaned FASTQ files.

    Supported examples:
      1128AJ2-PF_L1_1.cleaned.fq
      1128AJ2-PF_L1_2.cleaned.fq

      SRR30938054_1.cleaned.fastq
      SRR30938054_2.cleaned.fastq

      sample_1.cleaned.fq.gz
      sample_2.cleaned.fq.gz

    Ignored automatically:
      *_fastqc.html
      *_fastqc.zip
      multiqc_report.html
      RNA_seqkit_summary.tsv

    Returns:
      [(sample_name, read1_path, read2_path), ...]
    """
    fastq_dir = Path(fastq_dir)

    pattern = re.compile(
        r"^(?P<sample>.+)_(?P<read>[12])\.cleaned\.(?:fq|fastq)(?:\.gz)?$",
        re.IGNORECASE,
    )

    samples = {}

    for path in fastq_dir.iterdir():
        if not path.is_file():
            continue

        match = pattern.match(path.name)
        if not match:
            continue

        sample = match.group("sample")
        read = match.group("read")

        if sample not in samples:
            samples[sample] = {}

        samples[sample][read] = path

    paired = []

    for sample in sorted(samples):
        if "1" not in samples[sample] or "2" not in samples[sample]:
            raise SystemExit(
                f"ERROR: sample {sample} does not have both R1 and R2 FASTQ files"
            )

        paired.append((sample, samples[sample]["1"], samples[sample]["2"]))

    if not paired:
        raise SystemExit(f"ERROR: no paired cleaned FASTQ files found in {fastq_dir}")

    return paired


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Remap paired-end RNA-seq with annotation-aware STAR, "
            "then run StringTie and FEELnc for lncRNA annotation."
        )
    )

    parser.add_argument(
        "--genome",
        required=True,
        help="Genome FASTA file."
    )

    parser.add_argument(
        "--braker_gtf",
        required=True,
        help="BRAKER UTR GTF, usually braker_utr.gtf.gz."
    )

    parser.add_argument(
        "--fastq_dir",
        required=True,
        help="Directory containing cleaned paired-end RNA-seq FASTQ files."
    )

    parser.add_argument(
        "--out_dir",
        default=".",
        help=(
            "Parent output directory. The script will create a FEELnc "
            "directory beneath it. Default: current directory."
        )
    )

    parser.add_argument(
        "--prefix",
        default="JPeel6_5_FEELnc",
        help="Prefix for FEELnc output files."
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=32,
        help="Threads for STAR, StringTie, and FEELnc_filter."
    )

    parser.add_argument(
        "--monoex",
        default="0",
        help="FEELnc --monoex value. Default 0 removes monoexonic candidates."
    )

    parser.add_argument(
        "--make_gff3",
        action="store_true",
        help="Convert final FEELnc lncRNA GTF to GFF3 with AGAT."
    )

    parser.add_argument(
        "--force_star_index",
        action="store_true",
        help="Rebuild STAR index even if it already exists."
    )

    args = parser.parse_args()

    check_tool("STAR")
    check_tool("stringtie")
    check_tool("FEELnc_filter.pl")
    check_tool("FEELnc_codpot.pl")
    check_tool("FEELnc_classifier.pl")

    if args.make_gff3:
        check_tool("agat_convert_sp_gxf2gxf.pl")

    genome = Path(args.genome).resolve()
    braker_gtf = Path(args.braker_gtf).resolve()
    fastq_dir = Path(args.fastq_dir).resolve()

    parent_out_dir = Path(args.out_dir).resolve()
    out_dir = parent_out_dir / "FEELnc"

    if not genome.exists():
        sys.exit(f"ERROR: genome FASTA not found: {genome}")

    if not braker_gtf.exists():
        sys.exit(f"ERROR: BRAKER GTF not found: {braker_gtf}")

    if not fastq_dir.exists():
        sys.exit(f"ERROR: FASTQ directory not found: {fastq_dir}")

    input_dir = out_dir / "00_inputs"
    star_index_dir = out_dir / "01_star_index_braker_utr"
    star_bam_dir = out_dir / "02_star_mapping_braker_utr"
    stringtie_dir = out_dir / "03_stringtie"
    feelnc_dir = out_dir / "04_feelnc"
    codpot_dir = feelnc_dir / "codpot"
    final_dir = out_dir / "05_final"
    log_dir = out_dir / "logs"

    for directory in [
        input_dir,
        star_index_dir,
        star_bam_dir,
        stringtie_dir,
        feelnc_dir,
        codpot_dir,
        final_dir,
        log_dir,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    command_log = log_dir / "commands.log"

    pc_gtf = input_dir / "braker_utr.gtf"

    print(f"Main FEELnc directory: {out_dir}", flush=True)
    print(f"Preparing annotation GTF: {pc_gtf}", flush=True)

    copy_or_decompress(braker_gtf, pc_gtf)

    saindex = star_index_dir / "SAindex"

    if args.force_star_index or not saindex.exists():
        genome_saindex_nbases = compute_genomeSAindex_nbases(genome)

        run(
            [
                "STAR",
                "--runMode", "genomeGenerate",
                "--genomeFastaFiles", genome,
                "--genomeSAindexNbases", genome_saindex_nbases,
                "--runThreadN", args.threads,
                "--genomeDir", star_index_dir,
                "--sjdbGTFfile", pc_gtf,
            ],
            log_file=command_log,
            stderr_file=log_dir / "STAR_genomeGenerate.stderr",
        )

    else:
        print(f"STAR index already exists, skipping: {star_index_dir}", flush=True)

    samples = discover_paired_fastqs(fastq_dir)

    with open(input_dir / "samples.tsv", "w") as out:
        out.write("sample\tread1\tread2\n")
        for sample, read1, read2 in samples:
            out.write(f"{sample}\t{read1}\t{read2}\n")

    print(f"Found {len(samples)} paired-end samples.", flush=True)

    mapped_bams = []

    for sample, read1, read2 in samples:
        out_prefix = star_bam_dir / f"{sample}_2pass"
        bam_sorted = star_bam_dir / f"{sample}_2passAligned.sortedByCoord.out.bam"

        if bam_sorted.exists() and bam_sorted.stat().st_size > 0:
            print(f"STAR BAM exists, skipping mapping: {bam_sorted}", flush=True)
            mapped_bams.append((sample, bam_sorted))
            continue

        read_files_command = []

        if str(read1).endswith(".gz"):
            read_files_command = ["--readFilesCommand", "zcat"]

        run(
            [
                "STAR",
                "--runMode", "alignReads",
                "--twopassMode", "Basic",
                "--genomeDir", star_index_dir,
                "--runThreadN", args.threads,
                *read_files_command,
                "--readFilesIn", read1, read2,
                "--outSAMstrandField", "intronMotif",
                "--outFileNamePrefix", out_prefix,
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outBAMcompression", "10",
            ],
            log_file=command_log,
            stderr_file=log_dir / f"STAR_{sample}.stderr",
        )

        if not bam_sorted.exists() or bam_sorted.stat().st_size == 0:
            raise SystemExit(f"ERROR: STAR BAM was not created correctly: {bam_sorted}")

        mapped_bams.append((sample, bam_sorted))

    bam_list = input_dir / "star_braker_utr_bams.tsv"

    with open(bam_list, "w") as out:
        out.write("sample\tbam\n")
        for sample, bam in mapped_bams:
            out.write(f"{sample}\t{bam}\n")

    stringtie_gtf_list = stringtie_dir / "stringtie_gtf.list"

    with open(stringtie_gtf_list, "w") as list_out:
        for sample, bam in mapped_bams:
            sample_gtf = stringtie_dir / f"{sample}.gtf"

            run(
                [
                    "stringtie",
                    bam,
                    "-p", args.threads,
                    "-G", pc_gtf,
                    "-o", sample_gtf,
                ],
                log_file=command_log,
                stderr_file=log_dir / f"stringtie_{sample}.stderr",
            )

            list_out.write(str(sample_gtf) + "\n")

    merged_gtf = stringtie_dir / "stringtie_merged.gtf"

    run(
        [
            "stringtie",
            "--merge",
            "-p", args.threads,
            "-G", pc_gtf,
            "-o", merged_gtf,
            stringtie_gtf_list,
        ],
        log_file=command_log,
        stderr_file=log_dir / "stringtie_merge.stderr",
    )

    candidate_gtf = feelnc_dir / "candidate_lncRNA.gtf"

    run(
        [
            "FEELnc_filter.pl",
            "-i", merged_gtf,
            "-a", pc_gtf,
            f"--monoex={args.monoex}",
            "-p", args.threads,
        ],
        log_file=command_log,
        stdout_file=candidate_gtf,
        stderr_file=log_dir / "FEELnc_filter.stderr",
    )

    run(
        [
            "FEELnc_codpot.pl",
            "-i", candidate_gtf,
            "-a", pc_gtf,
            "-g", genome,
            "--mode=shuffle",
            "--outdir", codpot_dir,
            "-o", args.prefix,
        ],
        log_file=command_log,
        stdout_file=log_dir / "FEELnc_codpot.stdout",
        stderr_file=log_dir / "FEELnc_codpot.stderr",
    )

    lnc_gtf = codpot_dir / f"{args.prefix}.lncRNA.gtf"
    mrna_gtf = codpot_dir / f"{args.prefix}.mRNA.gtf"

    if not lnc_gtf.exists() or lnc_gtf.stat().st_size == 0:
        sys.exit(f"ERROR: expected FEELnc lncRNA output missing: {lnc_gtf}")

    classes_txt = feelnc_dir / f"{args.prefix}_classes.txt"

    run(
        [
            "FEELnc_classifier.pl",
            "-i", lnc_gtf,
            "-a", pc_gtf,
        ],
        log_file=command_log,
        stdout_file=classes_txt,
        stderr_file=log_dir / "FEELnc_classifier.stderr",
    )

    lnc_gff3 = final_dir / f"{args.prefix}.lncRNA.gff3"

    if args.make_gff3:
        run(
            [
                "agat_convert_sp_gxf2gxf.pl",
                "-g", lnc_gtf,
                "-o", lnc_gff3,
            ],
            log_file=command_log,
            stdout_file=log_dir / "agat_lncRNA.stdout",
            stderr_file=log_dir / "agat_lncRNA.stderr",
        )

    print("\nDONE", flush=True)

    print(f"Main FEELnc directory:   {out_dir}", flush=True)
    print(f"STAR index:              {star_index_dir}", flush=True)
    print(f"STAR BAM directory:      {star_bam_dir}", flush=True)
    print(f"BAM list:                {bam_list}", flush=True)
    print(f"StringTie merged GTF:    {merged_gtf}", flush=True)
    print(f"FEELnc candidate GTF:    {candidate_gtf}", flush=True)
    print(f"FEELnc lncRNA GTF:       {lnc_gtf}", flush=True)
    print(f"FEELnc mRNA-like GTF:    {mrna_gtf}", flush=True)
    print(f"FEELnc class table:      {classes_txt}", flush=True)

    if args.make_gff3:
        print(f"Converted lncRNA GFF3:   {lnc_gff3}", flush=True)

    print(f"Command log:             {command_log}", flush=True)


if __name__ == "__main__":
    main()

