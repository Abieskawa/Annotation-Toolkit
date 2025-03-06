#!/usr/bin/env python3
"""
Combined GFF processing tool with integrated functionalities:

1) Fix: An integrated GFF fixer that processes Infernal and tRNAscan-SE lines,
   rewriting them into a gene/child/exon hierarchy, removing "evalue=..." attributes,
   and ensuring the file begins with '##gff-version 3'.

2) Tblout-to-GFF Conversion & Filtering:  
   When the input file is an Infernal tblout file (specified via -I infernal or --tblout),
   the file is first filtered for overlapping hits. For each overlapping region (same query/chromosome and strand),
   only the hit with the lowest E-value (or highest bit score if tied) is retained.
   Optionally, the filtered tblout hits are saved (in GFF format) to a separate output file (via --filtered-out).
   Then the filtered entries are converted to GFF and fixed.

3) Sorting: The fixed GFF is concatenated and sorted by scaffold, start, and hierarchy.
   A warning is issued if duplicate lines are detected and removed.

Supported input formats (use -I FORMAT INPUT_FILE):
  - infernal   : Infernal tblout (conversion & filtering applied)
  - trna-scan  : GFF from tRNAscan-SE
  - braker3    : Braker3 GFF

Usage examples:
  - For Infernal tblout:
      -I infernal my_infernal.tblout --csv mydata.csv --filtered-out filtered.gff
  - For tRNAscan-SE GFF:
      -I trna-scan my_trnascan.gff --csv mydata.csv
  - For Braker3 GFF:
      -I braker3 my_braker3.gff --csv mydata.csv

Other options allow you to run only the fix step, only the sort step, or the full merge (fix then sort) pipeline.

gff3 ref: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#Date: 18 August 2020   #Version: 1.26

Note: Only fix mode uses GFF reference info; the braker3 format remains unchanged.
"""

import sys
import argparse
import csv
import logging
import re
from collections import defaultdict
import pandas as pd

# ========================
# Helper Functions
# ========================
def load_csv(csv_path):
    d = {}
    with open(csv_path, 'r', newline='') as infile:
        rdr = csv.DictReader(infile)
        for row in rdr:
            col_id = row['ID'].strip()
            acc = row['Accession'].strip()
            tlist = [x.strip() for x in row['Type'].split(';')]
            desc = row['Description'].strip()
            d[col_id] = {
                'Accession': acc,
                'Type': tlist,
                'Description': desc
            }
    return d

def remove_evalue(attr_str):
    if not attr_str or attr_str == '.':
        return ''
    parts = [p.strip() for p in attr_str.split(';')]
    filtered = [p for p in parts if not p.lower().startswith('evalue=')]
    return ';'.join(filtered)

def parse_attrs(a):
    d = {}
    if not a or a == '.':
        return d
    for piece in a.split(';'):
        piece = piece.strip()
        if '=' in piece:
            k, v = piece.split('=', 1)
            d[k.strip()] = v.strip()
    return d

def reconst_attrs(d):
    if not d:
        return '.'
    return ';'.join(f"{k}={v}" for k, v in d.items()) + ';'

def pick_first_non_gene(types_list):
    for t in types_list:
        if t.lower() != 'gene':
            return t
    return 'gene'

def normalize_type(s):
    return re.sub(r"[_.\-'/]", '', s.lower())

def generate_unique_id(seqid, ftype, counter, suffix=""):
    return f"{seqid}.{ftype}{counter}{suffix}"

# ========================
# Functions for GFF Fixing
# ========================
def fix_gff_lines(in_lines, csv_path):
    id_to_info = load_csv(csv_path)
    logging.info(f"Loaded CSV with {len(id_to_info)} entries.")
    final_lines = []
    counters = defaultdict(int)
    tRNA_counters = defaultdict(int)
    line_num = 0
    for raw_line in in_lines:
        line_num += 1
        line_str = raw_line.rstrip('\n')
        # Preserve version and comment lines
        if line_str.startswith('##gff-version 3'):
            final_lines.append(line_str)
            continue
        if not line_str or line_str.startswith('#'):
            final_lines.append(line_str)
            continue
        parts = line_str.split('\t')
        if len(parts) != 9:
            logging.warning(f"Skipping malformed line {line_num}: {line_str}")
            continue
        seqid, source, ftype_in, start, end, score, strand, phase, attr_in = parts
        cleaned = remove_evalue(attr_in)
        # Case A: Infernal lines (tblout conversion) => create gene, child, exon lines
        if source.lower() == 'infernal' and (ftype_in in id_to_info):
            info = id_to_info[ftype_in]
            if info['Type'] and info['Type'][0].lower() == 'gene':
                first_non_gene = pick_first_non_gene(info['Type'])
                gene_key = (seqid, first_non_gene + "_g")
                counters[gene_key] += 1
                gcount = counters[gene_key]
                gene_id = generate_unique_id(seqid, normalize_type(first_non_gene), gcount, "_g1")
                gene_attrs = {
                    'ID': gene_id,
                    'Name': gene_id,
                    'gene_biotype': first_non_gene
                }
                gene_line = "\t".join([seqid, source, "gene", start, end, score, strand, phase, reconst_attrs(gene_attrs)])
                final_lines.append(gene_line)
                child_key = (seqid, first_non_gene)
                counters[child_key] += 1
                cval2 = counters[child_key]
                child_id = generate_unique_id(seqid, normalize_type(first_non_gene), cval2, "")
                child_attrs = {
                    'ID': child_id,
                    'Name': child_id,
                    'Parent': gene_id
                }
                child_line = "\t".join([seqid, source, normalize_type(first_non_gene), start, end, score, strand, phase, reconst_attrs(child_attrs)])
                final_lines.append(child_line)
                exon_id = f"{seqid}.exon1"
                exon_attrs = {
                    'ID': exon_id,
                    'Name': exon_id,
                    'Parent': child_id
                }
                exon_line = "\t".join([seqid, source, "exon", start, end, score, strand, phase, reconst_attrs(exon_attrs)])
                final_lines.append(exon_line)
            continue
        # Case B: tRNAscan-SE GFF lines => create gene and rewrite tRNA line
        elif source.lower() == 'trnascan-se' and ftype_in.lower() == 'trna':
            tRNA_counters[(seqid, "child_tRNA")] += 1
            nval = tRNA_counters[(seqid, "child_tRNA")]
            tRNA_counters[(seqid, "gene_tRNA")] += 1
            mval = tRNA_counters[(seqid, "gene_tRNA")]
            if ftype_in in id_to_info:
                info = id_to_info[ftype_in]
                first_non_gene = pick_first_non_gene(info['Type'])
            else:
                first_non_gene = "tRNA"
            gene_id = f"{seqid}.trna{nval}_g{mval}"
            gene_attrs = {
                'ID': gene_id,
                'Name': gene_id,
                'gene_biotype': first_non_gene
            }
            gene_line = "\t".join([seqid, source, "gene", start, end, score, strand, phase, reconst_attrs(gene_attrs)])
            final_lines.append(gene_line)
            attr_d = parse_attrs(cleaned)
            new_tRNA_id = f"{seqid}.trna{nval}"
            attr_d['ID'] = new_tRNA_id
            attr_d['Parent'] = gene_id
            updated_col9 = reconst_attrs(attr_d)
            tRNA_line = "\t".join([seqid, source, ftype_in, start, end, score, strand, phase, updated_col9])
            final_lines.append(tRNA_line)
        else:
            c9 = cleaned if cleaned.strip() else '.'
            updated_line = "\t".join([seqid, source, ftype_in, start, end, score, strand, phase, c9])
            final_lines.append(updated_line)
    if not any(ln.startswith('##gff-version 3') for ln in final_lines):
        final_lines.insert(0, "##gff-version 3")
    return final_lines

# ========================
# Functions for tblout-to-GFF Conversion (Integrated)
# ========================
def convert_tblout_to_gff_lines(in_lines, args):
    logging.info("Converting tblout to GFF format...")
    gff_lines = []
    source = "cmscan" if args.cmscan else "cmsearch"
    if args.version:
        source += "-" + args.version
    if args.source:
        source = args.source
    for line in in_lines:
        if line.startswith("#"):
            continue
        line = line.rstrip("\n")
        fields = line.split()
        # Use fmt2 if specified (expected for Infernal tblout)
        if args.fmt2:
            if len(fields) < 27:
                logging.error(f"Expected at least 27 fields for fmt2 but got {len(fields)} in line: {line}")
                continue
            idx = fields[0]
            if args.cmscan:
                seqname = fields[3]
                seqaccn = fields[4]
                mdlname = fields[1]
                mdlaccn = fields[2]
            else:
                seqname = fields[1]
                seqaccn = fields[2]
                mdlname = fields[3]
                mdlaccn = fields[4]
            clan    = fields[5]
            mdl     = fields[6]
            mdlfrom = fields[7]
            mdlto   = fields[8]
            seqfrom = fields[9]
            seqto   = fields[10]
            strand  = fields[11]
            trunc   = fields[12]
            pass_val = fields[13]
            gc      = fields[14]
            bias    = fields[15]
            score   = fields[16]
            evalue  = fields[17]
            inc     = fields[18]
            olp     = fields[19]
            anyidx  = fields[20]
            anyfrct1= fields[21]
            anyfrct2= fields[22]
            winidx  = fields[23]
            winfrct1= fields[24]
            winfrct2= fields[25]
            desc    = fields[26]
            if len(fields) > 27:
                extra_desc = "_".join(fields[27:])
                desc = desc + "_" + extra_desc
        else:
            if len(fields) < 18:
                logging.error(f"Expected at least 18 fields for default format but got {len(fields)} in line: {line}")
                continue
            if args.cmscan:
                seqname = fields[2]
                seqaccn = fields[3]
                mdlname = fields[0]
                mdlaccn = fields[1]
            else:
                seqname = fields[0]
                seqaccn = fields[1]
                mdlname = fields[2]
                mdlaccn = fields[3]
            mdl     = fields[4]
            mdlfrom = fields[5]
            mdlto   = fields[6]
            seqfrom = fields[7]
            seqto   = fields[8]
            strand  = fields[9]
            trunc   = fields[10]
            pass_val = fields[11]
            gc      = fields[12]
            bias    = fields[13]
            score   = fields[14]
            evalue  = fields[15]
            inc     = fields[16]
            desc    = fields[17]
            if len(fields) > 18:
                extra_desc = "_".join(fields[18:])
                desc = desc + "_" + extra_desc
        if strand not in ["+", "-"]:
            logging.error(f"Invalid strand '{strand}' in line: {line}")
            continue
        try:
            score_val = float(score)
            evalue_val = float(evalue)
        except ValueError:
            logging.error(f"Score or evalue not a valid float in line: {line}")
            continue
        if args.min_bit is not None and score_val < args.min_bit:
            continue
        if args.max_evalue is not None and evalue_val > args.max_evalue:
            continue
        if args.ignore_trna and mdlaccn in ["RF00005", "RF01852"]:
            continue
        attributes = f"evalue={evalue}"
        if args.all:
            desc_prefix = "" if args.hidedesc else "desc="
            if args.fmt2:
                attributes = (f"idx={idx};seqaccn={seqaccn};mdlaccn={mdlaccn};clan={clan};mdl={mdl};"
                              f"mdlfrom={mdlfrom};mdlto={mdlto};trunc={trunc};pass={pass_val};"
                              f"gc={gc};bias={bias};inc={inc};olp={olp};anyidx={anyidx};anyfrct1={anyfrct1};"
                              f"anyfrct2={anyfrct2};winidx={winidx};winfrct1={winfrct1};winfrct2={winfrct2};"
                              f"{desc_prefix}{desc}")
            else:
                attributes = (f"seqaccn={seqaccn};mdlaccn={mdlaccn};mdl={mdl};"
                              f"mdlfrom={mdlfrom};mdlto={mdlto};trunc={trunc};pass={pass_val};"
                              f"gc={gc};bias={bias};inc={inc};{desc_prefix}{desc}")
        elif args.none:
            attributes = "-"
        elif args.desc:
            attributes = desc
        if args.extra:
            if attributes == "-":
                attributes = ""
            elif not attributes.endswith(";"):
                attributes += ";"
            attributes += args.extra + ";"
        try:
            if strand == "+":
                start = int(seqfrom)
                end = int(seqto)
            else:
                start = int(seqto)
                end = int(seqfrom)
        except ValueError:
            logging.error(f"Invalid seqfrom/seqto in line: {line}")
            continue
        if start > end:
            start, end = end, start
        score_formatted = f"{score_val:.1f}"
        gff_line = f"{seqname}\t{source}\t{mdlname}\t{start}\t{end}\t{score_formatted}\t{strand}\t.\t{attributes}"
        gff_lines.append(gff_line)
    return gff_lines

# ========================
# Overlapping Hit Filtering for Infernal tblout
# ========================
def filter_overlapping_hits(tblout_lines, args):
    logging.info("Filtering overlapping hits from Infernal tblout...")
    hits = []
    for line in tblout_lines:
        if line.startswith("#"):
            continue
        fields = line.split()
        try:
            if args.fmt2:
                if args.cmscan:
                    query = fields[3]
                    seqfrom = int(fields[9])
                    seqto = int(fields[10])
                    strand = fields[11]
                    bit = float(fields[16])
                    evalue = float(fields[17])
                else:
                    query = fields[1]
                    seqfrom = int(fields[9])
                    seqto = int(fields[10])
                    strand = fields[11]
                    bit = float(fields[16])
                    evalue = float(fields[17])
            else:
                if args.cmscan:
                    query = fields[2]
                    seqfrom = int(fields[7])
                    seqto = int(fields[8])
                    strand = fields[9]
                    bit = float(fields[14])
                    evalue = float(fields[15])
                else:
                    query = fields[0]
                    seqfrom = int(fields[7])
                    seqto = int(fields[8])
                    strand = fields[9]
                    bit = float(fields[14])
                    evalue = float(fields[15])
        except Exception as ex:
            logging.error(f"Error parsing line: {line}. Exception: {ex}")
            continue
        # Compute actual start/end based on strand.
        if strand == "+":
            start = seqfrom
            end = seqto
        else:
            start = seqto
            end = seqfrom
        hit = {
            "query": query,
            "strand": strand,
            "start": start,
            "end": end,
            "evalue": evalue,
            "bit": bit,
            "line": line
        }
        hits.append(hit)
    
    # Group hits by (query, strand)
    groups = defaultdict(list)
    for hit in hits:
        key = (hit["query"], hit["strand"])
        groups[key].append(hit)
    
    filtered_hits = []
    # For each group, sort by start and cluster overlapping intervals.
    for key, group_hits in groups.items():
        group_hits.sort(key=lambda h: h["start"])
        cluster = []
        current_end = None
        for hit in group_hits:
            if not cluster:
                cluster.append(hit)
                current_end = hit["end"]
            else:
                if hit["start"] <= current_end:
                    cluster.append(hit)
                    current_end = max(current_end, hit["end"])
                else:
                    best = min(cluster, key=lambda h: (h["evalue"], -h["bit"]))
                    filtered_hits.append(best)
                    cluster = [hit]
                    current_end = hit["end"]
        if cluster:
            best = min(cluster, key=lambda h: (h["evalue"], -h["bit"]))
            filtered_hits.append(best)
    
    filtered_lines = [hit["line"] for hit in filtered_hits]
    logging.info(f"Filtering complete. {len(filtered_lines)} hits retained after filtering.")
    return filtered_lines

# ========================
# Main Processing in Fix Module (with Step Control and New -I Input Format)
# ========================
def process_fix(args):
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    
    # Determine input source based on -I flag if provided.
    if args.I is not None:
        input_format = args.I[0].lower()
        input_file = args.I[1]
        with open(input_file, 'r') as f:
            input_lines = f.readlines()
        logging.info(f"Reading input from {input_file} with format '{input_format}'")
        if input_format == "infernal":
            tblout_flag = True
        elif input_format in ["trna-scan", "trnascan", "braker3"]:
            tblout_flag = False
        else:
            logging.error("Unsupported input format. Supported formats: infernal, trna-scan, braker3.")
            sys.exit(1)
    elif args.input:
        with open(args.input, 'r') as f:
            input_lines = f.readlines()
        logging.info(f"Reading input from {args.input}")
        tblout_flag = args.tblout
        input_format = "infernal" if tblout_flag else "gff"
    else:
        input_lines = sys.stdin.readlines()
        logging.info("Reading input from stdin.")
        tblout_flag = args.tblout
        input_format = "infernal" if tblout_flag else "gff"
    
    # If input is Infernal tblout, filter overlapping hits first.
    if tblout_flag and ((args.I is not None and input_format == "infernal") or (args.I is None)):
        filtered_lines = filter_overlapping_hits(input_lines, args)
        if args.filtered_out:
            filtered_gff = convert_tblout_to_gff_lines(filtered_lines, args)
            with open(args.filtered_out, 'w') as outf:
                for ln in filtered_gff:
                    outf.write(ln + "\n")
            logging.info(f"Filtered tblout entries saved to {args.filtered_out}")
        gff_lines = convert_tblout_to_gff_lines(filtered_lines, args)
    elif tblout_flag:
        gff_lines = convert_tblout_to_gff_lines(input_lines, args)
    else:
        gff_lines = input_lines

    # Process steps based on --step option.
    if args.step == "fix":
        fixed_lines = fix_gff_lines(gff_lines, args.csv)
        if not any(ln.startswith('##gff-version 3') for ln in fixed_lines):
            fixed_lines.insert(0, "##gff-version 3")
        if args.output:
            with open(args.output, 'w') as outf:
                for ln in fixed_lines:
                    outf.write(ln + "\n")
            logging.info(f"Fixed GFF file saved to {args.output}")
        else:
            for ln in fixed_lines:
                print(ln)
            logging.info("Fixed GFF output printed to stdout.")
    elif args.step == "sort":
        # Sorting only (assumes input is already fixed GFF).
        fixed_lines = gff_lines
        if not any(ln.startswith('##gff-version 3') for ln in fixed_lines):
            fixed_lines.insert(0, "##gff-version 3")
        header_lines = []
        data_lines = []
        for line in fixed_lines:
            if line.startswith("#"):
                header_lines.append(line.rstrip('\n'))
            else:
                parts = line.strip().split('\t')
                if len(parts) == 9:
                    data_lines.append(parts)
                else:
                    logging.warning(f"Skipping invalid line during sort: {line}")
        if data_lines:
            df = pd.DataFrame(data_lines, columns=['scaffold', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            df['start'] = pd.to_numeric(df['start'])
            df['end'] = pd.to_numeric(df['end'])
        else:
            df = pd.DataFrame(columns=['scaffold', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        parent_dict = {}
        for _, row in df.iterrows():
            if 'ID=' in row['attributes']:
                feature_id = row['attributes'].split('ID=')[1].split(';')[0]
                parent_dict[feature_id] = 0
        for _, row in df.iterrows():
            if 'Parent=' in row['attributes']:
                child_id = row['attributes'].split('ID=')[1].split(';')[0]
                parent_id = row['attributes'].split('Parent=')[1].split(';')[0]
                parent_dict[child_id] = parent_dict.get(parent_id, 0) + 1
        df['hierarchy'] = df.apply(lambda row: parent_dict.get(row['attributes'].split('Parent=')[1].split(';')[0], 0) + 1 if 'Parent=' in row['attributes'] else 0, axis=1)
        original_count = len(df)
        df_no_dup = df.drop_duplicates()
        if len(df_no_dup) < original_count:
            logging.warning("Duplicates found in GFF sorting and removed.")
        sorted_df = df_no_dup.sort_values(by=['scaffold', 'start', 'hierarchy'], ascending=[True, True, True])
        output_lines = header_lines[:]
        for _, row in sorted_df.iterrows():
            out_line = "\t".join([str(row[col]) for col in ['scaffold','source','feature','start','end','score','strand','phase','attributes']])
            output_lines.append(out_line)
        if args.output:
            with open(args.output, 'w') as outf:
                for ln in output_lines:
                    outf.write(ln + "\n")
            logging.info(f"Sorted GFF file saved to {args.output}")
        else:
            for ln in output_lines:
                print(ln)
            logging.info("Sorted GFF output printed to stdout.")
    elif args.step == "merge":
        fixed_lines = fix_gff_lines(gff_lines, args.csv)
        if not any(ln.startswith('##gff-version 3') for ln in fixed_lines):
            fixed_lines.insert(0, "##gff-version 3")
        header_lines = []
        data_lines = []
        for line in fixed_lines:
            if line.startswith("#"):
                header_lines.append(line.rstrip('\n'))
            else:
                parts = line.strip().split('\t')
                if len(parts) == 9:
                    data_lines.append(parts)
                else:
                    logging.warning(f"Skipping invalid line during sort: {line}")
        if data_lines:
            df = pd.DataFrame(data_lines, columns=['scaffold', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            df['start'] = pd.to_numeric(df['start'])
            df['end'] = pd.to_numeric(df['end'])
        else:
            df = pd.DataFrame(columns=['scaffold', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        parent_dict = {}
        for _, row in df.iterrows():
            if 'ID=' in row['attributes']:
                feature_id = row['attributes'].split('ID=')[1].split(';')[0]
                parent_dict[feature_id] = 0
        for _, row in df.iterrows():
            if 'Parent=' in row['attributes']:
                child_id = row['attributes'].split('ID=')[1].split(';')[0]
                parent_id = row['attributes'].split('Parent=')[1].split(';')[0]
                parent_dict[child_id] = parent_dict.get(parent_id, 0) + 1
        df['hierarchy'] = df.apply(lambda row: parent_dict.get(row['attributes'].split('Parent=')[1].split(';')[0], 0) + 1 if 'Parent=' in row['attributes'] else 0, axis=1)
        original_count = len(df)
        df_no_dup = df.drop_duplicates()
        if len(df_no_dup) < original_count:
            logging.warning("Duplicates found in GFF sorting and removed.")
        sorted_df = df_no_dup.sort_values(by=['scaffold', 'start', 'hierarchy'], ascending=[True, True, True])
        output_lines = header_lines[:]
        for _, row in sorted_df.iterrows():
            out_line = "\t".join([str(row[col]) for col in ['scaffold','source','feature','start','end','score','strand','phase','attributes']])
            output_lines.append(out_line)
        if args.output:
            with open(args.output, 'w') as outf:
                for ln in output_lines:
                    outf.write(ln + "\n")
            logging.info(f"Final sorted GFF file saved to {args.output}")
        else:
            for ln in output_lines:
                print(ln)
            logging.info("Final sorted GFF output printed to stdout.")

# ========================
# Main Entry Point
# ========================
def main():
    parser = argparse.ArgumentParser(
        description="Combined GFF processing tool that fixes GFF format, optionally converts Infernal tblout (with overlapping hit filtering) to GFF, "
                    "and allows you to test separate steps: fix, sort, or full merge (fix then sort).\n\n"
                    "Supported input formats (use -I FORMAT INPUT_FILE):\n"
                    "  infernal   : Infernal tblout (conversion & filtering applied)\n"
                    "  trna-scan  : GFF from tRNAscan-SE\n"
                    "  braker3    : Braker3 GFF\n\n"
                    "If -I is not provided, -i is used for input and --tblout indicates tblout format."
    )
    subparsers = parser.add_subparsers(dest='command', required=True, help='Sub-command to run')

    fix_parser = subparsers.add_parser('fix', help='Run GFF processing. Use --step to select processing stage.')
    fix_parser.add_argument("-I", nargs=2, metavar=("FORMAT", "INPUT_FILE"),
                            help="Specify the input format and file. Supported formats: infernal, trna-scan, braker3. For example: -I infernal my_infernal.tblout")
    fix_parser.add_argument("-i", "--input", type=str, default=None,
                            help="Input GFF or tblout file (if -I is not used).")
    # Set default CSV to Rfam_15_0.csv
    fix_parser.add_argument("-c", "--csv", type=str, default="Rfam_15_0.csv",
                            help="CSV with columns: Accession,ID,Type,Description (Rfam or similar). Default: Rfam_15_0.csv")
    fix_parser.add_argument("-o", "--output", type=str, default=None,
                            help="Output file for final result (or stdout).")
    fix_parser.add_argument("--tblout", action="store_true",
                            help="Specify that the input file is in tblout format (if -I is not used).")
    fix_parser.add_argument("--step", choices=["fix", "sort", "merge"], default="merge",
                            help="Select processing step: 'fix' to only fix, 'sort' to only sort (input must be fixed GFF), 'merge' (default) to run full pipeline.")
    fix_parser.add_argument("--filtered-out", type=str, default=None,
                            help="Output file for storing the filtered tblout entries (in GFF format) after overlapping filtering (only for infernal tblout).")
    # Tblout conversion options with improved flag names.
    fix_parser.add_argument("--min-bit", "-T", type=float, default=None,
                            help="Minimum bit score to include (tblout conversion).")
    fix_parser.add_argument("--max-evalue", "-E", type=float, default=None,
                            help="Maximum E-value to include (tblout conversion).")
    fix_parser.add_argument("--cmscan", action="store_true",
                            help="Indicate that the tblout file was created by cmscan (tblout conversion).")
    fix_parser.add_argument("--fmt2", action="store_true",
                            help="Indicate that the tblout file was created with --fmt2 option (expects at least 27 fields) (tblout conversion).")
    fix_parser.add_argument("--all", action="store_true",
                            help="Output all info in 'attributes' column for tblout conversion.")
    fix_parser.add_argument("--none", action="store_true",
                            help="Output no info in 'attributes' column for tblout conversion.")
    fix_parser.add_argument("--desc", action="store_true",
                            help="Output description field in 'attributes' column instead of E-value for tblout conversion.")
    fix_parser.add_argument("--version", type=str, default=None,
                            help="Append '-<version>' to the source field for tblout conversion.")
    fix_parser.add_argument("--extra", type=str, default=None,
                            help="Append '<extra>;' to the attributes column for tblout conversion.")
    fix_parser.add_argument("--hidedesc", action="store_true",
                            help="Do not include 'desc:' prior to description value in attributes column for tblout conversion.")
    fix_parser.add_argument("--source", type=str, default=None,
                            help="Specify source field to use instead of default (cmscan/cmsearch) for tblout conversion.")
    fix_parser.add_argument("--ignore-trna", action="store_true",
                            help="Ignore Infernal/Rfam predictions for tRNAs in tblout conversion.")

    args = parser.parse_args()

    # Validate tblout conversion options if needed.
    if (args.tblout or (args.I and args.I[0].lower() == "infernal")):
        if args.min_bit is not None and args.max_evalue is not None:
            parser.error("ERROR: --min-bit and --max-evalue cannot be used in combination. Pick one.")
        if (args.all and args.none) or (args.all and args.desc) or (args.none and args.desc):
            parser.error("ERROR: --all, --none, and --desc options are mutually exclusive.")
        if args.fmt2 and not args.cmscan:
            parser.error("ERROR: --fmt2 only makes sense with --cmscan.")
        if args.source and args.version:
            parser.error("ERROR: --source and --version are incompatible.")

    if args.command == 'fix':
        process_fix(args)

if __name__ == "__main__":
    main()
