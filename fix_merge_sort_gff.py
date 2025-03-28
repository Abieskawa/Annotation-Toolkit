#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import re
import shutil

def usage():
    print("Usage: {} -g <input_genome> -r <input_gff> -p <prefix> -a <annotator> [-v]".format(sys.argv[0]))
    sys.exit(1)

def run_command(cmd, shell=False):
    """Run a command and exit if it fails."""
    print("Running: " + " ".join(cmd) if not shell else cmd)
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        print("Error running command: ", cmd, file=sys.stderr)
        sys.exit(result.returncode)

def merge_intervals(intervals):
    """
    Given a list of intervals (start, end), merge overlapping intervals.
    Returns a new list of intervals (tuples) sorted by start.
    """
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def get_hierarchy_level(feature, attributes_str):
    """
    Return an integer level for a given feature type and its attribute string.
    Level 0: gene (or biological_region)
    Level 2: exon or CDS
    Level 1: if the attribute string contains "Parent=" (and not gene or exon/CDS), treat as transcript.
    """
    f = feature.lower()
    if f in ('gene', 'biological_region'):
         return 0
    elif f in ('exon', 'cds'):
         return 2
    elif "Parent=" in attributes_str:
         return 1
    else:
         return 0

def extract_bed_from_gff(gff_file, feature, bed_file, adjust_start=True, verbose=False):
    """
    Process a GFF file by grouping features into a hierarchy:
      Gene -> Transcript -> Exon/Introns

    This function uses get_hierarchy_level() to assign a level for each record.
    
    Grouping logic:
      - Level 0 records (gene) create a new gene group.
      - Level 1 records (transcript) are mapped to their parent gene (using the Parent attribute).
      - Level 2 records (exon/CDS) are stored under the transcript (using the Parent attribute).
        If no transcript record is encountered, they are stored directly under the gene.
      - Intron records (type "intron") are stored similarly as level 1 (if present).
    
    Then:
      - For feature=="gene": all exon intervals (from transcripts and direct gene-level) are merged and the overall gene span is output.
      - For feature=="transcript": one BED line is output per transcript using its merged exon intervals (flattened into one overall interval).
      - For feature=="intron": only intron intervals directly recorded under transcripts are output.
    
    The BED file is written with columns: chrom, start (0-based), end, ID, ., strand.
    
    If verbose is enabled, only gene groups that have a gene record but no exon-level data are printed.
    """
    groups = {}  # key: gene id, value: dict with keys "gene", "exons" (direct), and "transcripts"
    with open(gff_file) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attr = parts
            level = get_hierarchy_level(ftype, attr)
            # Parse attributes into a dictionary.
            attributes = {}
            for a in attr.split(";"):
                a = a.strip()
                if not a:
                    continue
                if "=" in a:
                    key, value = a.split("=", 1)
                    attributes[key] = value
            if level == 0:
                gene_id = attributes.get("ID") or attributes.get("gene_id")
                if gene_id:
                    groups.setdefault(gene_id, {"gene": None, "exons": [], "transcripts": {}})
                    groups[gene_id]["gene"] = (chrom, int(start), int(end), strand, attr)
                    if verbose:
                        print(f"DEBUG: Found gene: {gene_id} on {chrom}:{start}-{end} ({strand})")
            elif level == 1:
                transcript_id = attributes.get("ID")
                parent_gene = attributes.get("Parent") or attributes.get("gene_id")
                if transcript_id and parent_gene:
                    groups.setdefault(parent_gene, {"gene": None, "exons": [], "transcripts": {}})
                    groups[parent_gene]["transcripts"].setdefault(transcript_id, {"transcript": None, "exons": [], "introns": []})
                    groups[parent_gene]["transcripts"][transcript_id]["transcript"] = (chrom, int(start), int(end), strand, attr)
                    if verbose:
                        print(f"DEBUG: Found transcript: {transcript_id} for gene {parent_gene} on {chrom}:{start}-{end} ({strand})")
            elif level == 2:
                parent = attributes.get("Parent") or attributes.get("transcript_id")
                if parent:
                    # Check if parent is a transcript in any gene.
                    found = False
                    for gene_id, group in groups.items():
                        if parent in group["transcripts"]:
                            group["transcripts"][parent]["exons"].append((int(start), int(end)))
                            found = True
                            break
                    if not found:
                        # If no transcript record, store directly at gene level.
                        groups.setdefault(parent, {"gene": None, "exons": [], "transcripts": {}})
                        groups[parent]["exons"].append((int(start), int(end)))
            # For intron entries (we assume introns are level 1, though they might be typed explicitly as "intron")
            elif ftype.lower() == "intron":
                parent = attributes.get("Parent") or attributes.get("transcript_id")
                if parent:
                    for gene_id, group in groups.items():
                        if parent in group["transcripts"]:
                            group["transcripts"][parent].setdefault("introns", []).append((int(start), int(end)))
                            break
    # Print debug info only for gene groups with a gene record but no exon-level data.
    if verbose:
        for gene_id, group in groups.items():
            total_exons = len(group["exons"])
            for t_id, transcript in group["transcripts"].items():
                total_exons += len(transcript["exons"])
            if group["gene"] is not None and total_exons == 0:
                print(f"DEBUG: Gene {gene_id} has a gene record but no exon-level data.")
    
    with open(bed_file, "w") as fout:
        if feature.lower() == "gene":
            # Merge all exon intervals from both direct gene-level exons and from transcripts.
            for gene_id, group in groups.items():
                if group["gene"] is None:
                    continue
                all_exons = list(group["exons"])
                for transcript in group["transcripts"].values():
                    all_exons.extend(transcript["exons"])
                if not all_exons:
                    continue
                merged_exons = merge_intervals(all_exons)
                chrom, gstart, gend, strand, gene_attr = group["gene"]
                bed_start = gstart - 1 if adjust_start else gstart
                line_out = f"{chrom}\t{bed_start}\t{gend}\t{gene_id}\t.\t{strand}"
                fout.write(line_out + "\n")
        elif feature.lower() == "transcript":
            # For each transcript, merge its exon intervals and flatten into one overall interval.
            for gene_id, group in groups.items():
                for transcript_id, transcript in group["transcripts"].items():
                    if transcript["transcript"] is None or len(transcript["exons"]) == 0:
                        continue
                    merged_exons = merge_intervals(transcript["exons"])
                    tx_start = min(interval[0] for interval in merged_exons)
                    tx_end = max(interval[1] for interval in merged_exons)
                    chrom, tstart, tend, strand, t_attr = transcript["transcript"]
                    bed_start = tx_start - 1 if adjust_start else tx_start
                    line_out = f"{chrom}\t{bed_start}\t{tx_end}\t{transcript_id}\t.\t{strand}"
                    fout.write(line_out + "\n")
        elif feature.lower() == "intron":
            # Output only intron entries that were directly recorded under transcripts.
            for gene_id, group in groups.items():
                for transcript_id, transcript in group["transcripts"].items():
                    if transcript.get("introns") and len(transcript["introns"]) > 0:
                        sorted_introns = sorted(transcript["introns"], key=lambda x: x[0])
                        for i, (istart, iend) in enumerate(sorted_introns):
                            chrom = group["gene"][0]
                            bed_start = istart - 1 if adjust_start else istart
                            line_out = f"{chrom}\t{bed_start}\t{iend}\t{transcript_id}_intron{i+1}\t.\t{group['gene'][3]}"
                            fout.write(line_out + "\n")
    return

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline using a GFF input (converted to GTF) to generate BED/FASTA outputs at gene and transcript levels"
    )
    parser.add_argument("-g", "--genome", required=True, help="Input genome file")
    parser.add_argument("-r", "--gff", required=True, help="Input GFF file")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output files")
    parser.add_argument("-a", "--annotator", required=True, help="Annotator name (e.g. braker3)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose debug output")
    args = parser.parse_args()

    input_genome = args.genome
    input_gff = args.gff
    prefix = args.prefix
    annotator = args.annotator
    verbose = args.verbose

    # Define output file names.
    output_intron_lowercase_genome = f"{prefix}_genome_intron-lower.fa"
    intron_bed = f"{prefix}_{annotator}_intron.bed"
    gene_bed = f"{prefix}_{annotator}_gene.bed"
    trans_bed = f"{prefix}_{annotator}_trans.bed"
    gene_fa = f"{prefix}_{annotator}_gene.fa"
    trans_fa = f"{prefix}_{annotator}_trans.fa"
    cds_fa_temp = f"{prefix}_{annotator}_CDS_70.fa"
    pep_fa_temp = f"{prefix}_{annotator}_pep_70.fa"
    cds_fa = f"{prefix}_{annotator}_CDS.fa"
    pep_fa = f"{prefix}_{annotator}_pep.fa"
    aa_fa = f"{prefix}_{annotator}_aa.fa"

    # -----------------------------------------------------
    # Generate a GTF file from the input GFF using gffread.
    gtf_generated = "merged_fix_tran.gtf"
    cmd_gffread_gtf = ["gffread", "-T", input_gff, "-o", gtf_generated]
    run_command(cmd_gffread_gtf)
    print(f"GTF file generated: {gtf_generated}")
    # We'll use the generated GTF for CDS/peptide extraction.
    # For BED extractions, we work with the original GFF.
    # -----------------------------------------------------

    # Generate BED files.
    extract_bed_from_gff(input_gff, "intron", intron_bed, adjust_start=True, verbose=verbose)
    print(f"Intron BED file written: {intron_bed}")
    extract_bed_from_gff(input_gff, "gene", gene_bed, adjust_start=True, verbose=verbose)
    print(f"Gene BED file written: {gene_bed}")
    extract_bed_from_gff(input_gff, "transcript", trans_bed, adjust_start=True, verbose=verbose)
    print(f"Transcript BED file written: {trans_bed}")

    # Use bedtools maskfasta to soft-mask the genome using the intron BED file.
    cmd_mask = [
        "bedtools", "maskfasta", "-soft",
        "-fi", input_genome,
        "-bed", intron_bed,
        "-fo", output_intron_lowercase_genome
    ]
    run_command(cmd_mask)
    print(f"Soft-masked genome written: {output_intron_lowercase_genome}")
    # -----------------------------------------------------

    # Use bedtools getfasta to extract gene and transcript sequences.
    cmd_getfasta_gene = [
        "bedtools", "getfasta", "-fi", output_intron_lowercase_genome,
        "-bed", gene_bed, "-nameOnly", "-s", "-fo", gene_fa
    ]
    run_command(cmd_getfasta_gene)
    print(f"Gene FASTA written: {gene_fa}")
    cmd_getfasta_trans = [
        "bedtools", "getfasta", "-fi", output_intron_lowercase_genome,
        "-bed", trans_bed, "-nameOnly", "-s", "-fo", trans_fa
    ]
    run_command(cmd_getfasta_trans)
    print(f"Transcript FASTA written: {trans_fa}")

    # Remove any "(+)" and "(-)" from the FASTA headers.
    for fasta_file in [gene_fa, trans_fa]:
        with open(fasta_file, "r") as f:
            content = f.read()
        content = re.sub(r'\(\+\)', '', content)
        content = re.sub(r'\(-\)', '', content)
        with open(fasta_file, "w") as f:
            f.write(content)
        print(f"Cleaned {fasta_file} of stranded annotations.")

    # Run gffread (using the generated GTF) to extract CDS and peptide FASTA sequences.
    cmd_gffread = [
        "gffread", "-g", output_intron_lowercase_genome,
        "-x", cds_fa_temp,
        "-y", pep_fa_temp,
        gtf_generated
    ]
    run_command(cmd_gffread)
    print(f"gffread produced {cds_fa_temp} and {pep_fa_temp}")

    # Process CDS and peptide FASTA files (simulate awk/sort/awk pipeline).
    for temp_file in [cds_fa_temp, pep_fa_temp]:
        output_file = temp_file.replace("_70", "")
        records = []
        with open(temp_file, "r") as fin:
            header = None
            seq_lines = []
            for line in fin:
                if line.startswith(">"):
                    if header:
                        records.append((header, "".join(seq_lines)))
                    header = line.strip()
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
            if header:
                records.append((header, "".join(seq_lines)))
        records.sort(key=lambda x: x[0])
        with open(output_file, "w") as fout:
            for rec in records:
                fout.write(f"{rec[0]}\n{rec[1]}\n")
        print(f"Processed FASTA file written: {output_file}")

    # Rename peptide file: replace '.t' with '.p' in header if needed.
    with open(pep_fa_temp.replace("_70", ""), "r") as fin:
        content = fin.read()
    content = re.sub(r'\.t', '.p', content)
    with open(pep_fa, "w") as fout:
        fout.write(content)
    print(f"Peptide file renamed to: {pep_fa}")

    # Clean up temporary files.
    for f in [intron_bed, gene_bed, trans_bed, cds_fa_temp, pep_fa_temp]:
        if os.path.exists(f):
            os.remove(f)
    print("Temporary BED and intermediate FASTA files removed.")

    # Move output files into MOLAS input directory.
    output_dir = f"./{prefix}_MOLAS_input/FASTA/"
    os.makedirs(output_dir, exist_ok=True)
    files_to_move = [
        output_intron_lowercase_genome,
        output_intron_lowercase_genome + ".fai",
        gene_fa, trans_fa, cds_fa, pep_fa, aa_fa
    ]
    for f in files_to_move:
        if os.path.exists(f):
            shutil.move(f, os.path.join(output_dir, os.path.basename(f)))
    print(f"Moved output files to {output_dir}")

    # Remove the temporary GTF file.
    if os.path.exists(gtf_generated):
        os.remove(gtf_generated)
        print(f"Removed temporary GTF file: {gtf_generated}")

if __name__ == "__main__":
    main()
