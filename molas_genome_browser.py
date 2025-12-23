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
        bedToBigBed = shutil.which("bedToBigBed")
        if bedToBigBed:
            dir_path = os.path.dirname(bedToBigBed)
            fp = os.path.join(dir_path, filename)
            if os.path.exists(fp):
                return fp
            else:
                print(f"Warning: {fp} not found. Using '{filename}' from current directory.")
        return filename

def _normalize_transcript_type(s):
    """
    Minimal normalizer to keep your existing 'protein coding' style, while accepting NCBI values.
    """
    s = (s or "").strip()
    if not s:
        return "unknown"
    if s == "protein_coding":
        return "protein coding"
    # many NCBI noncoding types are already like rRNA/tRNA/misc_RNA; keep but unify underscores
    return s.replace("_", " ")

def get_gene_biotypes(gff_file):
    """
    MOD (minimal): NCBI gene lines have multiple identifiers (ID/Name/gene).
    Store biotype under all of them so transcriptType lookup works for NCBI GFF.
    ALSO: if gene_biotype is missing, fall back to gbkey on the gene feature.
    """
    biotypes = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue

            attrs = {}
            for attr in fields[8].split(";"):
                attr = attr.strip()
                if not attr or "=" not in attr:
                    continue
                key, value = attr.split("=", 1)
                attrs[key] = value

            # Key change: if no gene_biotype, use gbkey (NCBI often has gbkey=Gene/rRNA/etc.)
            biotype_raw = attrs.get("gene_biotype") or attrs.get("gbkey") or "unknown"
            biotype = _normalize_transcript_type(biotype_raw)

            for k in (attrs.get("ID"), attrs.get("Name"), attrs.get("gene")):
                if k:
                    biotypes[k] = biotype

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

    # Convert genePred to bigGenePred, sort, and cut first 16 columns (keep original behavior).
    bgp_bed = f"{prefix}_bgp.bed"
    genePredToBigGenePred = get_tool_path("genePredToBigGenePred", tools_dir)
    cmd = (
        f"{genePredToBigGenePred} -known {genePred_file} stdout | "
        "sort -k1,1 -k2,2n | cut -f 1-16 > " + bgp_bed
    )
    run_cmd(cmd, shell=True)

    return output_twobit, genePred_file, bgp_bed

def add_nr_annotation(nr_ann, bgp_bed, prefix, gene_biotypes):
    """
    Keep original behavior:
      - input bgp.bed has 16 columns
      - append 2 columns:
          col17: NR annotation (string)
          col18: transcriptType (string)
    """
    input_file = bgp_bed
    output_file = f"{prefix}_bgp-NR.bed"
    ann = {}
    with open(nr_ann, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] == "qseqid":
                continue
            key = row[0].split('.')[0]
            annotation = row[1] if len(row) > 1 else "_"
            if ">" in annotation:
                annotation = annotation.split('>')[0]
            ann[key] = annotation

    with open(input_file, 'r') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            gene = row[12] if len(row) > 12 else ""
            matched_ann = ann.get(gene, "_")
            transcript_type = gene_biotypes.get(gene, "protein coding")
            row.extend([matched_ann, transcript_type])
            writer.writerow(row)

    return output_file

def extract_genePred(genePred_file, prefix):
    """
    MOD (minimal but critical): Fix genePred columns used.
    genePred indices (0-based):
      0  name
      11 name2
      12 cdsStartStat
      13 cdsEndStat
      14 exonFrames
    """
    extracted_file = f"{prefix}_extracted-gp.tsv"
    with open(genePred_file, 'r') as fin, open(extracted_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if len(row) < 15:
                continue
            writer.writerow([row[0], row[11], row[12], row[13], row[14]])
    return extracted_file

def join_extracted_with_nr(extracted_file, nr_bed, prefix):
    """
    MOD: Put geneName/cdsStartStat/cdsEndStat/exonFrames back into columns 13-16
    (pre-reorder) so reorder_columns() yields your custom bigGenePred.as layout.
    """
    output_file = f"{prefix}_bgp-NR-CDS.bed"
    mapping = {}

    with open(extracted_file, 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        for row in reader:
            if len(row) < 5:
                continue
            mapping[row[0]] = [row[1], row[2], row[3], row[4]]  # geneName, cdsStartStat, cdsEndStat, exonFrames

    with open(nr_bed, 'r') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if len(row) < 18:
                row += [''] * (18 - len(row))
            tx = row[3]
            if tx in mapping:
                geneName, cdsStartStat, cdsEndStat, exonFrames = mapping[tx]
                row[12] = geneName
                row[13] = cdsStartStat
                row[14] = cdsEndStat
                row[15] = exonFrames
            writer.writerow(row)

    return output_file

def sed_replace(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
    content = content.replace("cmpl", "complete").replace("incmpl", "incomplete")
    with open(file_path, 'w') as f:
        f.write(content)

def reorder_columns(bed_file):
    """
    KEEP (unchanged): reorder into your custom bigGenePred.as column order.
    new order: 1..12, 17, 13, 14, 15, 18, 16
    """
    tmp_file = bed_file + ".tmp"
    with open(bed_file, 'r') as fin, open(tmp_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')
        for row in reader:
            if len(row) < 18:
                row += [''] * (18 - len(row))
            new_order = [row[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11,16,12,13,14,17,15]]
            writer.writerow(new_order)
    shutil.move(tmp_file, bed_file)

# --- Added (minimal): sanitizer to avoid bedToBigBed abort on NCBI edge cases ---

def _parse_list_field(s):
    s = (s or "").strip()
    if s.endswith(","):
        s = s[:-1]
    if s == "":
        return []
    return [x for x in s.split(",") if x != ""]

def _format_list_field(items):
    return ",".join(str(x) for x in items) + ","

def sanitize_lists_for_bigbed(bed_file):
    """
    Ensure list fields match blockCount:
      col10 blockCount
      col11 blockSizes
      col12 chromStarts
      col18 exonFrames

    Fix by truncating/padding; exonFrames pads with -1 ("no coding frame").
    """
    tmp_file = bed_file + ".san.tmp"
    fixed = 0

    with open(bed_file, 'r') as fin, open(tmp_file, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')

        for row in reader:
            if len(row) < 18:
                row += [''] * (18 - len(row))

            try:
                bc = int(row[9])
            except Exception:
                bc = 1
                row[9] = "1"

            bs = _parse_list_field(row[10])
            if len(bs) != bc:
                fixed += 1
                bs = (bs[:bc] + ["0"] * max(0, bc - len(bs)))
                row[10] = _format_list_field(bs)

            cs = _parse_list_field(row[11])
            if len(cs) != bc:
                fixed += 1
                cs = (cs[:bc] + ["0"] * max(0, bc - len(cs)))
                row[11] = _format_list_field(cs)

            ef = _parse_list_field(row[17])
            if len(ef) != bc:
                fixed += 1
                ef = (ef[:bc] + ["-1"] * max(0, bc - len(ef)))
                row[17] = _format_list_field(ef)

            writer.writerow(row)

    shutil.move(tmp_file, bed_file)
    print(f"sanitize_lists_for_bigbed: fixed {fixed} list-length issues in {bed_file}")

# --- end sanitizer ---

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
    """
    MOD (minimal): overwrite existing files if they already exist in destination.
    """
    if outdir:
        base_dir = os.path.abspath(outdir)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir, exist_ok=True)
    else:
        base_dir = os.getcwd()

    genome_browser_dir = os.path.join(base_dir, "genome_browser")
    intermediated_file_dir = os.path.join(genome_browser_dir, "intermediated_file")
    os.makedirs(genome_browser_dir, exist_ok=True)
    os.makedirs(intermediated_file_dir, exist_ok=True)

    # Remove temporary files.
    for f in [f"{prefix}_bgp-NR.bed", f"{prefix}_extracted-gp.tsv"]:
        if os.path.exists(f):
            os.remove(f)

    # Move intermediate files with overwrite.
    for f in [f"{prefix}_genome.2bit", f"{prefix}_genePred.gp", f"{prefix}_bgp.bed", f"{prefix}_bgp-NR-CDS.bed"]:
        if os.path.exists(f):
            dst = os.path.join(intermediated_file_dir, os.path.basename(f))
            if os.path.exists(dst):
                os.remove(dst)
            shutil.move(f, intermediated_file_dir)

    # Move the BigBed file with overwrite.
    bb_file = f"{prefix}_bgp-NR-CDS.bb"
    if os.path.exists(bb_file):
        dst = os.path.join(genome_browser_dir, os.path.basename(bb_file))
        if os.path.exists(dst):
            os.remove(dst)
        shutil.move(bb_file, genome_browser_dir)

def main():
    parser = argparse.ArgumentParser(
        description="Generate Genome BigBed File for MOLAS Genome Browser using GFF3 input."
    )
    parser.add_argument("-g", required=True, metavar="GENOME",
                        help="Input genome file in FASTA format (e.g., genome.fasta)")
    parser.add_argument("-r", required=True, metavar="GFF3",
                        help="Input GFF3 file (e.g., annotation.gff3)")
    parser.add_argument("-p", required=True, metavar="PREFIX",
                        help="Prefix for output files (e.g., beltfish_v2)")
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

    gene_biotypes = get_gene_biotypes(args.r)

    genome_2bit, genePred_file, bgp_bed = convert_genome_and_gff(args.g, args.r, prefix, tools_dir)
    nr_bed = add_nr_annotation(args.n, bgp_bed, prefix, gene_biotypes)
    extracted_file = extract_genePred(genePred_file, prefix)
    cds_bed = join_extracted_with_nr(extracted_file, nr_bed, prefix)
    sed_replace(cds_bed)

    # Keep your original reorder (required for your custom bigGenePred.as)
    reorder_columns(cds_bed)

    # Added: ensure lists satisfy bedToBigBed strict checks
    sanitize_lists_for_bigbed(cds_bed)

    create_bigbed(prefix, tools_dir)
    organize_files(prefix, outdir)

    print("BigBed file generation complete.")

if __name__ == "__main__":
    main()
