#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import math

def compute_genomeSAindex_nbases(genome_path):
    """
    Compute the STAR parameter genomeSAindexNbases based on the genome length.
    
    The value is calculated as:
        floor(log2(total_genome_length) / 2) - 1
    and then taking the minimum between 14 and that value.
    
    Parameters:
        genome_path (str): Path to the genome FASTA file.
    
    Returns:
        int: The computed value to be used for --genomeSAindexNbases.
    """
    total_length = 0
    with open(genome_path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                continue
            total_length += len(line.strip())
    # Calculate based on the formula: floor(log2(total_length)/2) - 1
    computed = int(math.log2(total_length) / 2) - 1
    return min(14, computed)

def main():
    # Parse command-line arguments (all parameters are required)
    parser = argparse.ArgumentParser(
        description="STAR mapping script converted from bash. All parameters are required."
    )
    parser.add_argument("--genomepath", required=True,
                        help="Path to the genome FASTA file")
    parser.add_argument("--genomedir", required=True,
                        help="Directory to hold the STAR genome index")
    parser.add_argument("--wd", required=True,
                        help="Working directory containing RNA-seq FASTQ files")
    parser.add_argument("--out_dir", required=True,
                        help="Output directory for RNA mapping")
    parser.add_argument("--threads", required=True, type=int,
                        help="Number of threads to use")
    args = parser.parse_args()

    # Log the full command-line invocation
    print("Command line:", " ".join(sys.argv))

    # Create necessary directories if they don't exist
    os.makedirs(args.genomedir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    # Check if the STAR genome index already exists by looking for 'SAindex'
    star_index_path = os.path.join(args.genomedir, "SAindex")
    if os.path.isdir(args.genomedir) and os.path.exists(star_index_path):
        print("STAR genome index already generated. Skipping.")
    else:
        print("STAR genome index not found. Generating...")
        # Compute genomeSAindexNbases from the genome FASTA file
        genome_saindex_nbases = compute_genomeSAindex_nbases(args.genomepath)
        star_cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeFastaFiles", args.genomepath,
            "--genomeSAindexNbases", str(genome_saindex_nbases),
            "--runThreadN", str(args.threads),
            "--genomeDir", args.genomedir
        ]
        print("Running command:", " ".join(star_cmd))
        subprocess.run(star_cmd, check=True)

    # Get unique basenames by scanning the working directory for files containing '.cleaned.fastq'
    unique_basenames = set()
    for filename in os.listdir(args.wd):
        if ".cleaned.fastq" in filename:
            # Mimic the bash 'cut' commands: take text before first dot, then text before first underscore.
            base = filename.split('.')[0]
            basename = base.split('_')[0]
            unique_basenames.add(basename)
    unique_basenames = sorted(unique_basenames)

    # Process each unique basename
    for basename in unique_basenames:
        paired1 = os.path.join(args.wd, f"{basename}_1.cleaned.fastq")
        paired2 = os.path.join(args.wd, f"{basename}_2.cleaned.fastq")
        single  = os.path.join(args.wd, f"{basename}.cleaned.fastq")

        out_prefix   = os.path.join(args.out_dir, f"{basename}_2pass")
        log_progress = f"{out_prefix}Log.progress.out"
        sj_tab       = f"{out_prefix}SJ.out.tab"
        bam_sorted   = os.path.join(args.out_dir, f"{basename}_2passAligned.sortedByCoord.out.bam")
        bam_flag2    = os.path.join(args.out_dir, f"{basename}_2pass_flag2.bam")

        # For paired-end FASTQ files
        if os.path.isfile(paired1) and os.path.isfile(paired2):
            if os.path.isfile(log_progress) and os.path.isfile(sj_tab):
                print(f"{basename}_1/2.cleaned.fastq has been mapped.")
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                print("Running command:", " ".join(samtools_cmd))
                subprocess.run(samtools_cmd, check=True)
            else:
                star_cmd = [
                    "STAR",
                    "--runMode", "alignReads",
                    "--twopassMode", "Basic",
                    "--genomeDir", args.genomedir,
                    "--runThreadN", str(args.threads),
                    "--readFilesIn", paired1, paired2,
                    "--outSAMstrandField", "intronMotif",
                    "--outFileNamePrefix", out_prefix,
                    "--outSAMtype", "BAM", "SortedByCoordinate"
                ]
                print("Running command:", " ".join(star_cmd))
                subprocess.run(star_cmd, check=True)
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                print("Running command:", " ".join(samtools_cmd))
                subprocess.run(samtools_cmd, check=True)

        # For single-end FASTQ file
        elif os.path.isfile(single):
            if os.path.isfile(log_progress) and os.path.isfile(sj_tab):
                print(f"{basename}.cleaned.fastq has been mapped.")
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                print("Running command:", " ".join(samtools_cmd))
                subprocess.run(samtools_cmd, check=True)
            else:
                star_cmd = [
                    "STAR",
                    "--runMode", "alignReads",
                    "--twopassMode", "Basic",
                    "--genomeDir", args.genomedir,
                    "--runThreadN", str(args.threads),
                    "--readFilesIn", single,
                    "--outSAMstrandField", "intronMotif",
                    "--outFileNamePrefix", out_prefix,
                    "--outSAMtype", "BAM", "SortedByCoordinate"
                ]
                print("Running command:", " ".join(star_cmd))
                subprocess.run(star_cmd, check=True)
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                print("Running command:", " ".join(samtools_cmd))
                subprocess.run(samtools_cmd, check=True)
        else:
            print(f"Error: Could not find paired or single-end FASTQ files for {basename}.")
            sys.exit(1)

if __name__ == '__main__':
    main()
