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

def extract_bed_from_gff(gff_file, feature, bed_file, adjust_start=True, verbose=False):
    """
    Build a hierarchy of gene/transcript/exon/intron from the GFF file and output BED entries.
    This function is used only for transcript-level BED generation since it requires
    assembling exons into transcripts.
    
    Hierarchy details:
      - Gene-level records (types "gene" or "ncrna_gene") are saved.
      - Exon-level records (exon, CDS, UTR) are collected and linked via Parent.
      - Intron-level records are also collected.
      - Any record (not gene, exon, CDS, UTR, intron) that has a Parent attribute is treated as transcript.
      
    Only transcripts whose parent gene exists are kept. BED entries for transcripts are produced
    by merging their exon intervals.
    """
    genes = {}         # key: gene id, value: (chrom, start, end, strand, attr)
    transcripts = {}   # key: transcript id, value: dict with keys:
                       #   "transcript": (chrom, start, end, strand, attr),
                       #   "parent_gene": gene id,
                       #   "exons": list of (start, end),
                       #   "introns": list of (start, end)
    gene_exons = {}    # For exon-level records directly attached to a gene

    with open(gff_file) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attr = parts
            start = int(start)
            end = int(end)
            # Parse attributes (assumes key=value; format)
            attributes = {}
            for a in attr.split(";"):
                a = a.strip()
                if not a:
                    continue
                if "=" in a:
                    key, value = a.split("=", 1)
                    attributes[key] = value

            ftype_low = ftype.lower()
            # Gene-level records.
            if ftype_low in ["gene", "ncrna_gene"]:
                gene_id = attributes.get("ID") or attributes.get("gene_id")
                if gene_id:
                    genes[gene_id] = (chrom, start, end, strand, attr)
            # Exon-level records.
            elif ftype_low in ["exon", "cds", "utr"]:
                parent = attributes.get("Parent")
                if parent:
                    if parent in transcripts:
                        transcripts[parent]["exons"].append((start, end))
                    elif parent in genes:
                        gene_exons.setdefault(parent, []).append((start, end))
                    else:
                        transcripts.setdefault(parent, {"transcript": None, "exons": [], "introns": []})
                        transcripts[parent]["exons"].append((start, end))
            # Intron-level records.
            elif ftype_low == "intron":
                parent = attributes.get("Parent")
                if parent:
                    transcripts.setdefault(parent, {"transcript": None, "exons": [], "introns": []})
                    transcripts[parent]["introns"].append((start, end))
            # Transcript-level records: no limitation on type
            elif ("Parent" in attributes or "gene_id" in attributes) and attributes.get("ID"):
                transcript_id = attributes.get("ID")
                parent_gene = attributes.get("Parent") or attributes.get("gene_id")
                transcripts.setdefault(transcript_id, {"transcript": None, "exons": [], "introns": []})
                transcripts[transcript_id]["transcript"] = (chrom, start, end, strand, attr)
                transcripts[transcript_id]["parent_gene"] = parent_gene

    # Filter transcripts: keep only those whose parent gene exists in genes.
    valid_transcripts = {tid: t for tid, t in transcripts.items() 
                         if ("parent_gene" in t and t["parent_gene"] in genes) or (tid in genes)}

    if verbose:
        for gene_id, gene_data in genes.items():
            total_exons = len(gene_exons.get(gene_id, []))
            for t in valid_transcripts.values():
                if t.get("parent_gene") == gene_id:
                    total_exons += len(t["exons"])
            if total_exons == 0:
                print(f"DEBUG: Gene {gene_id} has a gene record but no exon-level data.")

    with open(bed_file, "w") as fout:
        if feature.lower() == "transcript":
            for transcript_id, t in valid_transcripts.items():
                if t["transcript"] is None or len(t["exons"]) == 0:
                    continue
                merged_exons = merge_intervals(t["exons"])
                tx_start = min(interval[0] for interval in merged_exons)
                tx_end = max(interval[1] for interval in merged_exons)
                chrom, tstart, tend, strand, t_attr = t["transcript"]
                bed_start = tx_start - 1 if adjust_start else tx_start
                line_out = f"{chrom}\t{bed_start}\t{tx_end}\t{transcript_id}\t.\t{strand}"
                fout.write(line_out + "\n")
        elif feature.lower() == "intron":
            # Although intron extraction can be done directly, here we use the hierarchy if available.
            for transcript_id, t in valid_transcripts.items():
                if t.get("introns") and len(t["introns"]) > 0:
                    sorted_introns = sorted(t["introns"], key=lambda x: x[0])
                    for i, (istart, iend) in enumerate(sorted_introns):
                        if t["transcript"]:
                            chrom = t["transcript"][0]
                            strand = t["transcript"][3]
                        else:
                            parent_gene = t.get("parent_gene")
                            chrom, _, _, strand, _ = genes[parent_gene]
                        bed_start = istart - 1 if adjust_start else istart
                        line_out = f"{chrom}\t{bed_start}\t{iend}\t{transcript_id}_intron{i+1}\t.\t{strand}"
                        fout.write(line_out + "\n")
        elif feature.lower() == "gene":
            for gene_id, gene_data in genes.items():
                combined_exons = list(gene_exons.get(gene_id, []))
                for t in valid_transcripts.values():
                    if t.get("parent_gene") == gene_id:
                        combined_exons.extend(t["exons"])
                if not combined_exons:
                    continue
                merged_exons = merge_intervals(combined_exons)
                chrom, gstart, gend, strand, gene_attr = gene_data
                bed_start = gstart - 1 if adjust_start else gstart
                line_out = f"{chrom}\t{bed_start}\t{gend}\t{gene_id}\t.\t{strand}"
                fout.write(line_out + "\n")
    return

def extract_bed_direct(gff_file, feature, bed_file, adjust_start=True):
    """
    Directly extract BED lines for features that do not require hierarchy building.
    For example, for gene and intron, simply scan the GFF file and output the corresponding lines.
    """
    with open(gff_file) as fin, open(bed_file, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attr = parts
            ftype_low = ftype.lower()
            if feature.lower() == "gene":
                if ftype_low not in ["gene", "ncrna_gene"]:
                    continue
            elif feature.lower() == "intron":
                if ftype_low != "intron":
                    continue
            else:
                continue  # only handle gene and intron here
            start = int(start)
            end = int(end)
            bed_start = start - 1 if adjust_start else start
            # Parse attributes for an ID.
            attributes = {}
            for a in attr.split(";"):
                a = a.strip()
                if not a:
                    continue
                if "=" in a:
                    key, value = a.split("=", 1)
                    attributes[key] = value
            feature_id = attributes.get("ID", ".")
            fout.write(f"{chrom}\t{bed_start}\t{end}\t{feature_id}\t.\t{strand}\n")
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
    # For gene and intron features, we extract directly to save time.
    extract_bed_direct(input_gff, "intron", intron_bed, adjust_start=True)
    print(f"Intron BED file written: {intron_bed}")
    extract_bed_direct(input_gff, "gene", gene_bed, adjust_start=True)
    print(f"Gene BED file written: {gene_bed}")
    # For transcript, we need the heavy hierarchy build.
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

    # Process CDS and peptide FASTA files (simulate the awk/sort/awk pipeline).
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
