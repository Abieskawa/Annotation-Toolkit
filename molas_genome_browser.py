#!/usr/bin/env python3
import argparse
import subprocess
import os
import shutil
import csv
import sys

# Increase CSV field size limit to avoid errors.
csv.field_size_limit(sys.maxsize)

def run_cmd(cmd, shell=False):
    # Print the command as a space-separated string.
    if isinstance(cmd, list):
        print("Running command:", " ".join(cmd))
    else:
        print("Running command:", cmd)
    subprocess.run(cmd, shell=shell, check=True)

def get_tool_path(tool, tools_dir):
    """
    If tools_dir is provided, check if the tool exists there.
    Otherwise, try to locate the tool using shutil.which().
    If not found, return the tool name (hoping it is in the PATH).
    """
    if tools_dir:
        candidate = os.path.join(tools_dir, tool)
        if os.path.exists(candidate):
            return candidate
        else:
            print(f"Warning: {candidate} not found. Using '{tool}' from PATH.")
            return tool
    else:
        candidate = shutil.which(tool)
        if candidate:
            return candidate
        else:
            print(f"Warning: Could not locate {tool} in PATH; using '{tool}'.")
            return tool

def get_file_path(filename, tools_dir):
    """
    If tools_dir is provided, check if the file exists there.
    Otherwise, try to locate the file in the directory of bedToBigBed.
    If found, return its absolute path. Otherwise, return the filename.
    """
    if tools_dir:
        fp = os.path.join(tools_dir, filename)
        if os.path.exists(fp):
            return fp
        else:
            print(f"Warning: {fp} not found. Using '{filename}' from current directory.")
            return filename
    else:
        # Try to find bedToBigBed in PATH and then use its directory
        bedToBigBed = shutil.which("bedToBigBed")
        if bedToBigBed:
            dir_path = os.path.dirname(bedToBigBed)
            fp = os.path.join(dir_path, filename)
            if os.path.exists(fp):
                return fp
            else:
                print(f"Warning: {fp} not found. Using '{filename}' from current directory.")
        return filename

def get_gene_biotypes(gff_file):
    """
    Parse the input GFF3 file and return a mapping of gene ID to gene_biotype.
    If no gene_biotype is found, default to "protein coding".
    If the gene_biotype is "protein_coding", convert it to "protein coding";
    otherwise, replace underscores with spaces.
    """
    biotypes = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature == "gene":
                attr_field = fields[8]
                attrs = {}
                for attr in attr_field.split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    if "=" in attr:
                        key, value = attr.split("=", 1)
                        attrs[key] = value
                gene_id = attrs.get("ID", None)
                biotype = attrs.get("gene_biotype", "protein coding")
                if biotype == "protein_coding":
                    biotype = "protein coding"
                else:
                    biotype = biotype.replace("_", " ")
                if gene_id:
                    biotypes[gene_id] = biotype
    return biotypes

def convert_genome_and_gff(input_genome, input_gff, prefix, tools_dir):
    # Convert genome to 2bit.
    output_twobit = f"{prefix}_genome.2bit"
    faToTwoBit = get_tool_path("faToTwoBit", tools_dir)
    run_cmd([faToTwoBit, input_genome, output_twobit])
    
    # Convert GFF3 to genePred file.
    genePred_file = f"{prefix}_genePred.gp"
    gff3ToGenePred = get_tool_path("gff3ToGenePred", tools_dir)
    run_cmd([gff3ToGenePred, "-useName", input_gff, genePred_file])
    
    # Convert genePred to bigGenePred, sort, and cut first 16 columns.
    bgp_bed = f"{prefix}_bgp.bed"
    genePredToBigGenePred = get_tool_path("genePredToBigGenePred", tools_dir)
    cmd = (
        f"{genePredToBigGenePred} -known {genePred_file} stdout | "
        "sort -k1,1 -k2,2n | cut -f 1-16 > " + bgp_bed
    )
    run_cmd(cmd, shell=True)
    
    return output_twobit, genePred_file, bgp_bed

def add_nr_annotation(nr_ann, bgp_bed, prefix, gene_biotypes):
    input_file = bgp_bed
    output_file = f"{prefix}_bgp-NR.bed"
    ann = {}
    
    # First pass: read the NR annotation file.
    with open(nr_ann, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] == "qseqid":
                continue
            key = row[0].split('.')[0]
            annotation = row[1]
            if ">" in annotation:
                annotation = annotation.split('>')[0]
            ann[key] = annotation
    
    # Second pass: process bgp.bed file and add NR annotation and transcript type.
    with open(input_file, 'r') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            # Column 13 (index 12) is assumed to be gene ID.
            gene = row[12]
            matched_ann = ann.get(gene, "_")
            # Get gene_biotype from GFF3 mapping; default to "protein coding" if not found.
            transcript_type = gene_biotypes.get(gene, "protein coding")
            row.extend([matched_ann, transcript_type])
            writer.writerow(row)
    return output_file

def extract_genePred(genePred_file, prefix):
    extracted_file = f"{prefix}_extracted-gp.tsv"
    with open(genePred_file, 'r') as fin, open(extracted_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if len(row) < 15:
                continue
            # Extract columns: 1, 13, 14, 15 (indices 0, 12, 13, 14)
            writer.writerow([row[0], row[12], row[13], row[14]])
    return extracted_file

def join_extracted_with_nr(extracted_file, nr_bed, prefix):
    output_file = f"{prefix}_bgp-NR-CDS.bed"
    mapping = {}
    with open(extracted_file, 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        for row in reader:
            if len(row) < 4:
                continue
            mapping[row[0]] = [row[1], row[2], row[3]]
    
    with open(nr_bed, 'r') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if len(row) < 16:
                row += [''] * (16 - len(row))
            key = row[3]  # 4th column as key.
            if key in mapping:
                row[13] = mapping[key][0]
                row[14] = mapping[key][1]
                row[15] = mapping[key][2]
            writer.writerow(row)
    return output_file

def sed_replace(file_path):
    # Replace "cmpl" with "complete" and "incmpl" with "incomplete".
    with open(file_path, 'r') as f:
        content = f.read()
    content = content.replace("cmpl", "complete").replace("incmpl", "incomplete")
    with open(file_path, 'w') as f:
        f.write(content)

def reorder_columns(bed_file):
    """
    Reorder columns as:
    new order: 1,2,3,4,5,6,7,8,9,10,11,12,17,13,14,15,18,16.
    """
    tmp_file = bed_file + ".tmp"
    with open(bed_file, 'r') as fin, open(tmp_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if len(row) < 18:
                row += [''] * (18 - len(row))
            # Reorder indices: 0,1,2,3,4,5,6,7,8,9,10,11,16,12,13,14,17,15.
            new_order = [row[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11,16,12,13,14,17,15]]
            writer.writerow(new_order)
    shutil.move(tmp_file, bed_file)

def create_bigbed(prefix, tools_dir):
    bed_file = f"{prefix}_bgp-NR-CDS.bed"
    genome_2bit = f"{prefix}_genome.2bit"
    output_bb = f"{prefix}_bgp-NR-CDS.bb"
    as_file = get_file_path("bigGenePred.as", tools_dir)
    bedToBigBed = get_tool_path("bedToBigBed", tools_dir)
    cmd = [
        bedToBigBed,
        "-tab",
        f"-as={as_file}",
        "-type=bed12+6",
        "-sizesIs2Bit",
        bed_file,
        genome_2bit,
        output_bb
    ]
    run_cmd(cmd)
    return output_bb

def organize_files(prefix, outdir):
    # Determine base output directory. If outdir is provided, use it; otherwise, use current directory.
    if outdir:
        base_dir = os.path.abspath(outdir)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir, exist_ok=True)
    else:
        base_dir = os.getcwd()
    
    # Create the directory structure under base_dir.
    browser_input_dir = os.path.join(base_dir, f"{prefix}_MOLAS_input", "genome_browser")
    intermediate_dir = os.path.join(browser_input_dir, "intermediated_file")
    os.makedirs(intermediate_dir, exist_ok=True)
    
    # Remove temporary files.
    for f in [f"{prefix}_bgp-NR.bed", f"{prefix}_extracted-gp.tsv"]:
        if os.path.exists(f):
            os.remove(f)
    # Move intermediate files.
    for f in [f"{prefix}_genome.2bit", f"{prefix}_genePred.gp", f"{prefix}_bgp.bed", f"{prefix}_bgp-NR-CDS.bed"]:
        if os.path.exists(f):
            shutil.move(f, intermediate_dir)
    # Move the BigBed file.
    bb_file = f"{prefix}_bgp-NR-CDS.bb"
    if os.path.exists(bb_file):
        shutil.move(bb_file, browser_input_dir)

def main():
    parser = argparse.ArgumentParser(
        description="Generate Genome BigBed File for MOLAS Genome Browser using GFF3 input."
    )
    parser.add_argument("-g", required=True, metavar="GENOME", 
                        help="Input genome file in FASTA format (e.g., genome.fasta)")
    parser.add_argument("-r", required=True, metavar="GFF3", 
                        help="Input GFF3 file (e.g., annotation.gff3)")
    parser.add_argument("-p", required=True, metavar="PREFIX", 
                        help="Prefix for output files (e.g., FT2.1)")
    parser.add_argument("-n", required=True, metavar="NR_ANNOT", 
                        help="NR annotation file in TSV format")
    parser.add_argument("--tools_dir", default=None, metavar="TOOLS_DIR",
                        help="(Optional) Directory containing UCSC tools and bigGenePred.as. "
                             "If not provided, tools are assumed to be in the PATH and bigGenePred.as is searched in the bedToBigBed directory.")
    parser.add_argument("-o", "--outdir", default=None, metavar="OUTPUT_DIR",
                        help="(Optional) Uppermost output directory. If provided, output files will be written below this directory.")
    args = parser.parse_args()
    
    prefix = args.p
    tools_dir = args.tools_dir
    outdir = args.outdir
    
    # Parse gene biotypes from the GFF3 file.
    global gene_biotypes
    gene_biotypes = get_gene_biotypes(args.r)
    
    # Step 1: Convert genome and GFF3 to intermediate formats.
    genome_2bit, genePred_file, bgp_bed = convert_genome_and_gff(args.g, args.r, prefix, tools_dir)
    
    # Step 2: Add NR annotation to bgp.bed file.
    nr_bed = add_nr_annotation(args.n, bgp_bed, prefix, gene_biotypes)
    
    # Step 3: Extract columns from genePred file.
    extracted_file = extract_genePred(genePred_file, prefix)
    
    # Step 4: Join extracted info into NR BED to create NR-CDS BED.
    cds_bed = join_extracted_with_nr(extracted_file, nr_bed, prefix)
    
    # Step 5: Replace "cmpl"/"incmpl" with "complete"/"incomplete".
    sed_replace(cds_bed)
    
    # Step 6: Reorder columns as specified.
    reorder_columns(cds_bed)
    
    # Step 7: Create BigBed file using bedToBigBed.
    create_bigbed(prefix, tools_dir)
    
    # Step 8: Organize output and intermediate files.
    organize_files(prefix, outdir)
    
    print("BigBed file generation complete.")

if __name__ == "__main__":
    main()
