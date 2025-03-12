#!/usr/bin/env python3
"""
Combined GFF tool:
  - Fixes GFF from various sources (infernal, trnascan‑se, braker3, mitos) to create
    gene/child/exon hierarchies and ensures only 9 columns with a single "##gff-version 3".
  - "fix" mode: processes each file separately => <basename>_fix.gff.
  - "all" mode: processes all files => for each file also writes an <basename>_fix.gff,
    then merges them into one merged_fix.gff.

Assumption:
  - The exon number of ncRNA genes predicted with cmscan is only one (remove tRNA entries in this script)

References:
1. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
2. https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
3. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/
"""

import sys
import argparse
import csv
import logging
import re
import os
from collections import defaultdict
import pandas as pd

def load_csv(fp):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not os.path.isabs(fp):
        fp_in_script = os.path.join(script_dir, fp)
        if os.path.exists(fp_in_script):
            fp = fp_in_script
    if not os.path.exists(fp):
        fp_in_cwd = os.path.join(os.getcwd(), fp)
        if os.path.exists(fp_in_cwd):
            fp = fp_in_cwd
    if not os.path.exists(fp):
        raise FileNotFoundError(f"CSV file not found: {fp}")
    d = {}
    with open(fp, 'r', newline='') as f:
        for row in csv.DictReader(f):
            key = row['ID'].strip()
            d[key] = {
                'Accession': row['Accession'].strip(),
                'Type': [t.strip() for t in row['Type'].split(';')],
                'Description': row['Description'].strip()
            }
    return d

def parse_attrs(a):
    d = {}
    if not a or a == '.':
        return d
    for p in a.split(';'):
        p = p.strip()
        if '=' in p:
            k, v = p.split('=', 1)
            d[k.strip()] = v.strip()
    return d

def reconst_attrs(d):
    if not d:
        return '.'
    return ';'.join(f"{k}={v}" for k, v in d.items()) + ';'

# Modified: Allow period (".") to remain in the string.
def normalize_type(s):
    return re.sub(r"[\-'/]", '', s.lower())

def gen_uid(seq, typ, cnt, basename=None):
    if basename:
        return f"{basename}_{normalize_type(typ)}{cnt}"
    else:
        return f"{seq}_{normalize_type(typ)}{cnt}"

# ---------------- Mitos-specific ----------------

def fix_mitos_attr_field(attr, seq, args):
    new_attrs = []
    for chunk in attr.split(';'):
        chunk = chunk.strip()
        if not chunk:
            continue
        if '=' in chunk:
            key, value = chunk.split('=', 1)
            if key in ("ID", "Parent", "gene_id", "Name"):
                if args.basename:
                    if not value.startswith(f"{args.basename}_"):
                        value = f"{args.basename}_{value}"
                else:
                    if not value.startswith(f"{seq}."):
                        value = f"{seq}.{value}"
            new_attrs.append(f"{key}={value}")
        else:
            new_attrs.append(chunk)
    out = ';'.join(new_attrs)
    if out and not out.endswith(';'):
        out += ';'
    return out

def fix_mitos_lines(lines, args):
    fixed = []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#"):
            fixed.append(line)
            continue
        parts = line.split('\t')
        if len(parts) != 9:
            fixed.append(line)
            continue
        if parts[2] == "ncRNA_gene":
            parts[2] = "gene"
        parts[8] = fix_mitos_attr_field(parts[8], parts[0], args)
        fixed.append("\t".join(parts))
    return fixed

# ---------------- Infernal Tblout → GFF ----------------

def convert_tblout_to_gff(lines, args):
    accepted = []
    dropped = []
    if args.source:
        src = args.source
    else:
        src = "cmscan"
    if args.version:
        src += "-" + args.version
    for line in lines:
        if line.startswith("#"):
            continue
        fs = line.rstrip("\n").split()
        if args.fmt2:         
            if len(fs) < 27:
                dropped.append((line, "Insufficient fields"))
                continue
            seq = fs[3]
            s_from, s_to = fs[9], fs[10]
            strand, score, ev = fs[11], fs[16], fs[17]
            feature = fs[1]
        else:
            if len(fs) < 18:
                dropped.append((line, "Insufficient fields"))
                continue
            seq = fs[2]
            s_from, s_to = fs[7], fs[8]
            strand, score, ev = fs[9], fs[14], fs[15]
            feature = fs[2]
        # Filter out tRNA entries: if "trna" appears anywhere in the feature name.
        if args.ignore_trna and "trna" in feature.lower():
            dropped.append((line, "Ignored tRNA (due to --ignore-trna)"))
            continue
        try:
            if strand == "+":
                s, e = int(s_from), int(s_to)
            else:
                s, e = int(s_to), int(s_from)
            score_f = ev
        except Exception as ex:
            dropped.append((line, f"Coord parse error: {ex}"))
            continue
        if args.min_bit is not None and float(score) < args.min_bit:
            dropped.append((line, "Bit score below threshold"))
            continue
        if args.max_evalue is not None and float(ev) > args.max_evalue:
            dropped.append((line, "E-value above threshold"))
            continue
        accepted.append("\t".join([seq, src, feature, str(s), str(e), score_f, strand, ".", "."]))
    return accepted, dropped

# ---------------- Braker3 fix ----------------

def fix_braker3(lines, args):
    new_lines = []
    for line in lines:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith("#"):
            new_lines.append(line)
            continue
        parts = line.split('\t')
        if len(parts) != 9:
            new_lines.append(line)
            continue
        seq, src, feat, start, end, score, strand, phase, attr = parts
        if feat.lower() == "gene":
            if attr == ".":
                attr = "gene_biotype=protein_coding;"
            elif "gene_biotype=" not in attr:
                if not attr.endswith(";"):
                    attr += ";"
                attr += "gene_biotype=protein_coding;"
        new_lines.append("\t".join([seq, src, feat, start, end, score, strand, phase, attr]))
    return new_lines

# ---------------- GFF fix for infernal / tRNAscan‑se ----------------

def fix_gff_lines(lines, csv_fp, in_fmt, args):
    csvdata = load_csv(csv_fp)
    out = []
    cnts, tcnts = defaultdict(int), defaultdict(int)
    for line in lines:
        if line.startswith("#") or not line.strip():
            out.append(line.rstrip("\n"))
            continue
        parts = line.rstrip("\n").split('\t')
        if len(parts) != 9:
            out.append(line.rstrip("\n"))
            continue
        seq, src, typ, start, end, score, strand, phase, attr = parts
        if in_fmt == "infernal":
            if typ in csvdata:
                info = csvdata[typ]
                toks = info['Type']
                if any(x.lower() == "gene" for x in toks):
                    # Use the second token if available; otherwise, use the CSV key (ID) as proxy.
                    child_tok = toks[1] if len(toks) > 1 else typ
                    cnts[(seq, child_tok.lower()+"_g")] += 1
                    gid = gen_uid(seq, child_tok+"_g", cnts[(seq, child_tok.lower()+"_g")], basename=args.basename)
                    pattrs = {'ID': gid, 'gene_biotype': child_tok, 'description': info['Description']}
                    out.append("\t".join([seq, src, "gene", start, end, score, strand, phase, reconst_attrs(pattrs)]))
                    cnts[(seq, child_tok.lower())] += 1
                    cid = gen_uid(seq, child_tok, cnts[(seq, child_tok.lower())], basename=args.basename)
                    cattrs = {'ID': cid, 'Parent': gid}
                    out.append("\t".join([seq, src, child_tok, start, end, score, strand, phase, reconst_attrs(cattrs)]))
                    exon_id = f"{cid}.exon1"
                    exon_attrs = {'ID': exon_id, 'Parent': cid}
                    out.append("\t".join([seq, src, "exon", start, end, score, strand, phase, reconst_attrs(exon_attrs)]))
                else:
                    cnts[(seq, "bio_reg")] += 1
                    bid = gen_uid(seq, "biological_region", cnts[(seq, "bio_reg")], basename=args.basename)
                    battrs = {'ID': bid, 'description': info['Description']}
                    out.append("\t".join([seq, src, "biological_region", start, end, score, strand, phase, reconst_attrs(battrs)]))
            else:
                out.append("\t".join(parts[:8] + [attr if attr else "."]))
        elif in_fmt == "trnascan-se":
            if typ in csvdata:
                info = csvdata[typ]
                toks = info['Type']
                if any(x.lower() == "gene" for x in toks):
                    child_tok = toks[1] if len(toks) > 1 else typ
                    cnts[(seq, child_tok.lower()+"_g")] += 1
                    gid = gen_uid(seq, child_tok+"_g", cnts[(seq, child_tok.lower()+"_g")], basename=args.basename)
                    pattrs = {'ID': gid, 'gene_biotype': child_tok, 'description': info['Description']}
                    out.append("\t".join([seq, src, "gene", start, end, score, strand, phase, reconst_attrs(pattrs)]))
                    d = parse_attrs(attr)
                    d['Parent'] = gid
                    if args.basename:
                        if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                            d['ID'] = f"{args.basename}.{d['ID']}"
                    out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
                else:
                    d = parse_attrs(attr)
                    if args.basename:
                        if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                            d['ID'] = f"{args.basename}.{d['ID']}"
                        if 'Parent' in d and not d['Parent'].startswith(f"{args.basename}"):
                            d['Parent'] = f"{args.basename}.{d['Parent']}"
                    out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
            else:
                d = parse_attrs(attr)
                if args.basename:
                    if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                        d['ID'] = f"{args.basename}.{d['ID']}"
                    if 'Parent' in d and not d['Parent'].startswith(f"{args.basename}"):
                        d['Parent'] = f"{args.basename}.{d['Parent']}"
                out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
        else:
            out.append("\t".join(parts[:8] + [attr if attr else "."]))
    final = []
    for l in out:
        if l.startswith("##gff-version"):
            continue
        final.append(l)
    final.insert(0, "##gff-version 3")
    return final

def fix_gff_lines_main(lines, csv_fp, in_fmt, args):
    if in_fmt == "mitos":
        return fix_mitos_lines(lines, args)
    elif in_fmt == "braker3":
        return fix_braker3(lines, args)
    elif in_fmt in ("infernal", "trnascan-se"):
        return fix_gff_lines(lines, csv_fp, in_fmt, args)
    else:
        return lines

# ---------------- Cleanup ----------------

def cleanup_for_fix(lines):
    out = []
    version_found = False
    for line in lines:
        if line.startswith("#!gff-spec-version"):
            continue
        if line.startswith("##gff-version"):
            if not version_found:
                out.append("##gff-version 3")
                version_found = True
            continue
        out.append(line)
    final = []
    for line in out:
        if line.startswith("#"):
            final.append(line)
            continue
        parts = line.split("\t")
        if len(parts) >= 9:
            final.append("\t".join(parts[:9]))
        else:
            final.append(line)
    return final

# ---------------- Hierarchy & Sorting ----------------

def get_hierarchy_level(feature, attributes):
    f = feature.lower()
    if f in ('gene', 'biological_region'):
         return 0
    elif f in ('exon', 'cds'):
         return 2
    elif "Parent=" in attributes:
         return 1
    else:
         return 0

def sort_merge_gff_lines(lines):
    headers = []
    data = []
    for l in lines:
        if l.startswith("#"):
            if l.startswith("##gff-version"):
                continue
            headers.append(l)
        else:
            data.append(l)
    rows = []
    for l in data:
        parts = l.split("\t")
        if len(parts) < 9:
            continue
        rows.append(parts[:9])
    if not rows:
        return ["##gff-version 3"] + headers
    df = pd.DataFrame(rows, columns=['scaffold','source','feature','start','end','score','strand','phase','attributes'])
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df = df.dropna(subset=['start'])
    df['hierarchy'] = df.apply(lambda row: get_hierarchy_level(row['feature'], row['attributes']), axis=1)
    df = df.sort_values(by=['scaffold','start','hierarchy'])
    df = df.drop(columns=['hierarchy'])
    merged = ["\t".join(map(str, row)) for row in df.values]
    unique_hdr = []
    for h in headers:
        if h not in unique_hdr:
            unique_hdr.append(h)
    final = ["##gff-version 3"] + unique_hdr + merged
    return final

# ---------------- File Processing ----------------

def process_file(fmt, fname, args):
    with open(fname, 'r') as f:
        lines = f.readlines()
    input_count = len(lines)
    dropped = []
    fmt_lower = fmt.lower()
    conv_count = None
    if fmt_lower == "infernal":
        converted, drops2 = convert_tblout_to_gff(lines, args)
        conv_count = len(converted)
        lines = converted
        dropped = drops2
    elif fmt_lower in ("braker3", "mitos", "trnascan-se"):
        pass
    fixed = fix_gff_lines_main(lines, args.csv, fmt_lower, args)
    return fixed, dropped, input_count, conv_count

# ---------------- fix mode ----------------

def process_fix_cmd(args):
    outputs = {}
    drop_logs = {}
    metrics = {}
    for (fmt, fname) in args.I:
        fixed, dropped, input_count, conv_count = process_file(fmt, fname, args)
        cleaned = cleanup_for_fix(fixed)
        outputs[(fmt, fname)] = cleaned
        drop_logs[(fmt, fname)] = dropped
        metrics[(fmt, fname)] = (input_count, conv_count)
    return outputs, drop_logs, metrics

# ---------------- all mode ----------------

def process_all(args):
    non_mitos_lines = []
    mitos_lines = []
    all_drops = {}
    metrics = {}
    for (fmt, fname) in args.I:
        fixed, dropped, input_count, conv_count = process_file(fmt, fname, args)
        cleaned = cleanup_for_fix(fixed)
        base = os.path.basename(fname)
        root, _ = os.path.splitext(base)
        out_fname = os.path.join(args.outdir, f"{root}_{fmt}_fix.gff")
        with open(out_fname, 'w') as outf:
            outf.write("\n".join(cleaned) + "\n")
        print(f"----\nWrote intermediate all-fix file: {out_fname} with {len(cleaned)} lines.")
        if fmt.lower() == "mitos":
            mitos_lines.extend(cleaned)
        else:
            non_mitos_lines.extend(cleaned)
        all_drops[fname] = dropped
        metrics[(fmt, fname)] = (input_count, conv_count)
    merged_non_mitos = sort_merge_gff_lines(non_mitos_lines)
    mitos_clean = [line for line in mitos_lines if not line.startswith("##gff-version")]
    merged = merged_non_mitos + mitos_clean
    return merged, all_drops, metrics

# ---------------- main ----------------

def main():
    parser = argparse.ArgumentParser(
        description="GFF combine & fix tool. In 'all' mode, also produce intermediate _fix.gff files."
    )
    sub = parser.add_subparsers(dest='command', required=True)
    
    fixp = sub.add_parser('fix', help="Process each file separately => <root>_fix.gff.")
    fixp.add_argument("-I", nargs=2, metavar=("FORMAT","FILE"), action="append", required=True,
                      help="(FORMAT in {infernal, trnascan-se, braker3, mitos})")
    fixp.add_argument("--csv", type=str, default="Rfam_15_0.csv",
                      help="CSV for lookup (default=Rfam_15_0.csv).")
    fixp.add_argument("--outdir", type=str, required=True,
                      help="Output directory for separate <root>_fix.gff.")
    fixp.add_argument("--save-filtered-hits", type=str, default=None,
                      help="Prefix for dropped-hits log (infernal only).")
    fixp.add_argument("--min-bit", "-T", type=float, default=None,
                      help="Min bit score for infernal.")
    fixp.add_argument("--max-evalue", "-E", type=float, default=None,
                      help="Max evalue for infernal.")
    fixp.add_argument("--fmt2", action="store_true", help="Infernal tblout with --fmt2 (>=27 fields).")
    fixp.add_argument("--version", type=str, default=None, help="Append version to source.")
    fixp.add_argument("--ignore-trna", action="store_true", help="Skip tRNA hits in infernal.")
    fixp.add_argument("--basename", type=str, default=None,
                      help="Prepend to generated gene IDs, e.g. 'M_tai' => 'M_tai_rrna.g1'.")
    fixp.add_argument("--source", type=str, default=None,
                      help="Specify source field for infernal-to-GFF conversion. If provided, this value is used in the second column; otherwise, default is 'cmscan'.")
    
    allp = sub.add_parser('all', help="Fix all => <root>_fix.gff, then merge => merged_fix.gff.")
    allp.add_argument("-I", nargs=2, metavar=("FORMAT","FILE"), action="append", required=True,
                      help="(FORMAT in {infernal, trnascan-se, braker3, mitos})")
    allp.add_argument("--csv", type=str, default="Rfam_15_0.csv",
                      help="CSV for lookup (default=Rfam_15_0.csv).")
    allp.add_argument("--outdir", type=str, required=True,
                      help="Output directory. Also writes merged_fix.gff.")
    allp.add_argument("--save-filtered-hits", type=str, default=None,
                      help="Prefix for dropped-hits log (infernal only).")
    allp.add_argument("--min-bit", "-T", type=float, default=None,
                      help="Min bit score for infernal.")
    allp.add_argument("--max-evalue", "-E", type=float, default=None,
                      help="Max evalue for infernal.")
    allp.add_argument("--fmt2", action="store_true", help="Infernal tblout with --fmt2 (>=27 fields).")
    allp.add_argument("--version", type=str, default=None, help="Append version to source.")
    allp.add_argument("--ignore-trna", action="store_true", help="Skip tRNA hits in infernal.")
    allp.add_argument("--basename", type=str, default=None,
                      help="Prepend to generated gene IDs, e.g. 'M_tai' => 'M_tai_rrna.g1'.")
    allp.add_argument("--source", type=str, default=None,
                      help="Specify source field for infernal-to-GFF conversion. If provided, this value is used in the second column; otherwise, default is 'cmscan'.")
    
    args = parser.parse_args()
    csv_data = load_csv(args.csv)
    csv_count = len(csv_data)
    
    if args.command == "fix":
        outputs, drop_logs, metrics = process_fix_cmd(args)
        summary_by_type = defaultdict(list)
        for (fmt, in_fname), lines in outputs.items():
            input_count, conv_count = metrics[(fmt, in_fname)]
            base = os.path.basename(in_fname)
            root, _ = os.path.splitext(base)
            out_fname = os.path.join(args.outdir, f"{root}_{fmt}_fix.gff")
            with open(out_fname, 'w') as outF:
                outF.write("\n".join(lines) + "\n")
            output_count = len(lines)
            dropped_count = len(drop_logs[(fmt, in_fname)])
            drop_reason = ""
            if fmt.lower() == "infernal" and args.ignore_trna:
                drop_reason = "tRNA entries ignored due to --ignore-trna"
            elif dropped_count > 0:
                ds = defaultdict(int)
                for ln, reason in drop_logs[(fmt, in_fname)]:
                    ds[reason] += 1
                drop_reason = ", ".join(f"{v} {k}" for k, v in ds.items())
            if fmt.lower() == "infernal":
                parent_count = sum(1 for line in lines if not line.startswith("#") and line.split("\t")[2].lower() == "gene")
                exon_count = sum(1 for line in lines if not line.startswith("#") and line.split("\t")[2].lower() == "exon")
                msg = (f"Infernal Processing\n\n"
                       f"Input File: {in_fname}\n"
                       f"Total Lines: {input_count}\n\n"
                       f"Processing Steps:\n"
                       f"  - Converted to GFF: {conv_count} lines\n"
                       f"  - Dropped Entries: {dropped_count} ({drop_reason})\n"
                       f"  - Loaded Rfam CSV: {csv_count} rows\n\n"
                       f"Output File: {os.path.basename(out_fname)}\n"
                       f"Total Lines: {output_count}\n"
                       f"Drop Log: {args.save_filtered_hits if args.save_filtered_hits else f'{root}_dropout.log'}\n\n"
                       f"Added {parent_count} parent lines, {exon_count} exon lines\n"
                       f"------------------------------------------------------------")
                summary_by_type["infernal"].append(msg)
            elif fmt.lower() == "trnascan-se":
                msg = (f"tRNAscan-SE Processing\n\n"
                       f"Input File: {in_fname}\n"
                       f"Total Lines: {input_count}\n\n"
                       f"Processing Steps:\n"
                       f"  - Preserved original naming logic and added basename (and parent gene if applicable)\n"
                       f"  - Dropped Entries: {dropped_count} ({drop_reason})\n\n"
                       f"Output File: {os.path.basename(out_fname)}\n"
                       f"Total Lines: {output_count}\n"
                       f"------------------------------------------------------------")
                summary_by_type["trnascan-se"].append(msg)
            else:
                msg = (f"{fmt.upper()} Processing\n\n"
                       f"Input File: {in_fname}\n"
                       f"Total Lines: {input_count}\n\n"
                       f"Output File: {os.path.basename(out_fname)}\n"
                       f"Total Lines: {output_count}\n"
                       f"------------------------------------------------------------")
                summary_by_type[fmt.lower()].append(msg)
        for file_type, summaries in summary_by_type.items():
            print(f"\n==== Summary for {file_type.upper()} Files ====\n")
            for msg in summaries:
                print(msg)
    
    elif args.command == "all":
        merged, all_drops, metrics = process_all(args)
        merged_fname = os.path.join(args.outdir, "merged_fix.gff")
        with open(merged_fname, 'w') as outF:
            outF.write("\n".join(merged) + "\n")
        merged_count = len(merged)
        print(f"\nFinal Merging\n\nMerged File Created: {os.path.basename(merged_fname)}\nTotal Lines in Merged File: {merged_count}\n------------------------------------------------------------")
        summary_by_type = defaultdict(list)
        for (fmt, fn) in args.I:
            if fmt.lower() in ("infernal", "trnascan-se"):
                input_count, conv_count = metrics[(fmt, fn)]
                base = os.path.basename(fn)
                root, _ = os.path.splitext(base)
                intermediate_fname = os.path.join(args.outdir, f"{root}_{fmt}_fix.gff")
                with open(intermediate_fname, 'r') as inf:
                    intermediate_lines = [line.strip() for line in inf.readlines()]
                output_count = len(intermediate_lines)
                drops = all_drops[fn]
                dropped_count = len(drops)
                drop_reason = ""
                if fmt.lower() == "infernal" and args.ignore_trna:
                    drop_reason = "tRNA entries ignored due to --ignore-trna"
                elif dropped_count > 0:
                    ds = defaultdict(int)
                    for ln, reason in drops:
                        ds[reason] += 1
                    drop_reason = ", ".join(f"{v} {k}" for k, v in ds.items())
                if fmt.lower() == "infernal":
                    parent_count = sum(1 for line in intermediate_lines if not line.startswith("#") and line.split("\t")[2].lower() == "gene")
                    exon_count = sum(1 for line in intermediate_lines if not line.startswith("#") and line.split("\t")[2].lower() == "exon")
                    msg = (f"Infernal Processing\n\n"
                           f"Input File: {fn}\n"
                           f"Total Lines: {input_count}\n\n"
                           f"Processing Steps:\n"
                           f"  - Converted to GFF: {conv_count} lines\n"
                           f"  - Dropped Entries: {dropped_count} ({drop_reason})\n"
                           f"  - Loaded Rfam CSV: {csv_count} rows\n\n"
                           f"Output File: {os.path.basename(intermediate_fname)}\n"
                           f"Total Lines: {output_count}\n"
                           f"Drop Log: {args.save_filtered_hits if args.save_filtered_hits else f'{root}_dropout.log'}\n\n"
                           f"Added {parent_count} parent lines, {exon_count} exon lines\n"
                           f"------------------------------------------------------------")
                    summary_by_type["infernal"].append(msg)
                elif fmt.lower() == "trnascan-se":
                    msg = (f"tRNAscan-SE Processing\n\n"
                           f"Input File: {fn}\n"
                           f"Total Lines: {input_count}\n\n"
                           f"Processing Steps:\n"
                           f"  - Preserved original naming logic and added basename (and parent gene if applicable)\n"
                           f"  - Dropped Entries: {dropped_count} ({drop_reason})\n\n"
                           f"Output File: {os.path.basename(intermediate_fname)}\n"
                           f"Total Lines: {output_count}\n"
                           f"------------------------------------------------------------")
                    summary_by_type["trnascan-se"].append(msg)
        for file_type, summaries in summary_by_type.items():
            print(f"\n==== Summary for {file_type.upper()} Files ====\n")
            for msg in summaries:
                print(msg)
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
