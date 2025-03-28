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

def extract_bed_from_gff(gff_file, features, bed_files, adjust_start=True, verbose=False):
    """
    Build a hierarchy of gene/transcript/exon/intron from the GFF file and output BED entries
    for the specified features in one pass.
    
    Features can include: "gene", "intron", "transcript".
    
    For transcript-level BED, the transcript boundaries are taken directly from the transcript record,
    without merging exon intervals.
    """
    # Ensure features is a list
    if isinstance(features, str):
        features = [features]
        
    genes = {}         # key: gene id, value: (chrom, start, end, strand, attr)
    transcripts = {}   # key: transcript id, value: dict with keys:
                       #   "transcript": (chrom, start, end, strand, attr),
                       #   "parent_gene": gene id,
                       #   "exons": list of (start, end),
                       #   "introns": list of (start, end)
    # Removed gene_exons since exon-level records directly attached to genes are unlikely.
    warned_genes = set()  # Track genes for which a warning has already been printed

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
            if ftype_low in ["gene"]:
                gene_id = attributes.get("ID") or attributes.get("gene_id")
                if gene_id:
                    genes[gene_id] = (chrom, start, end, strand, attr)
            # Exon-level records (including CDS, UTR, start_codon, stop_codon).
            elif ftype_low in ["exon", "cds", "utr", "start_codon", "stop_codon"]:
                parent = attributes.get("Parent")
                if parent:
                    if parent in transcripts:
                        transcripts[parent]["exons"].append((start, end))
                    elif parent in genes:
                        # Warn once per gene and ignore storing these exons.
                        if parent not in warned_genes:
                            print(f"WARNING: Exon-level record directly attached to gene {parent} encountered. Ignoring.", file=sys.stderr)
                            warned_genes.add(parent)
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

    # Filter transcripts: only keep those with a valid parent gene.
    valid_transcripts = {}
    for tid, t in transcripts.items():
        if "parent_gene" in t:
            if t["parent_gene"] in genes:
                valid_transcripts[tid] = t
            else:
                print(f"WARNING: Transcript {tid} has parent_gene {t['parent_gene']} that is not found in genes. Skipping.", file=sys.stderr)
        else:
            print(f"WARNING: Transcript {tid} does not have a parent_gene attribute. Skipping.", file=sys.stderr)

    if verbose:
        for gene_id, gene_data in genes.items():
            # Only count exons from transcript-level records.
            total_exons = 0
            for t in valid_transcripts.values():
                if t.get("parent_gene") == gene_id:
                    total_exons += len(t["exons"])
            if total_exons == 0:
                print(f"DEBUG: Gene {gene_id} has a gene record but no exon-level data.")

    # For each requested feature, write to its corresponding BED file.
    for feat in features:
        feat_lower = feat.lower()
        with open(bed_files[feat_lower], "w") as fout:
            if feat_lower == "transcript":
                for transcript_id, t in valid_transcripts.items():
                    if t["transcript"] is None:
                        continue
                    # Use transcript record boundaries directly without merging exons.
                    chrom, tstart, tend, strand, t_attr = t["transcript"]
                    bed_start = tstart - 1 if adjust_start else tstart
                    line_out = f"{chrom}\t{bed_start}\t{tend}\t{transcript_id}\t.\t{strand}"
                    fout.write(line_out + "\n")
            elif feat_lower == "intron":
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
            elif feat_lower == "gene":
                for gene_id, gene_data in genes.items():
                    combined_exons = []  # No direct gene_exon storage used.
                    for t in valid_transcripts.values():
                        if t.get("parent_gene") == gene_id:
                            combined_exons.extend(t["exons"])
                    if not combined_exons:
                        continue
                    
                    chrom, gstart, gend, strand, gene_attr = gene_data
                    bed_start = gstart - 1 if adjust_start else gstart
                    line_out = f"{chrom}\t{bed_start}\t{gend}\t{gene_id}\t.\t{strand}"
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

    # Generate BED files for intron, gene, and transcript in one pass.
    bed_files = {
        "intron": intron_bed,
        "gene": gene_bed,
        "transcript": trans_bed
    }
    features = ["intron", "gene", "transcript"]
    extract_bed_from_gff(input_gff, features, bed_files, adjust_start=True, verbose=verbose)
    print(f"Intron BED file written: {intron_bed}")
    print(f"Gene BED file written: {gene_bed}")
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
