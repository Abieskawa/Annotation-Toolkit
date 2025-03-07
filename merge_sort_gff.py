#!/usr/bin/env python3
"""
Combined GFF tool:
  - Fixes GFF from various sources (infernal, trnascan‑se, braker3, mitos) to create
    gene/child/exon hierarchies and ensures only 9 columns with a single "##gff-version 3".
  - "fix" mode: processes each file separately => <basename>_fix.gff.
  - "all" mode: processes all files => for each file also writes an <basename>_allfix.gff,
    then merges them into one merged_fix.gff.
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

def remove_evalue(s):
    if not s or s == '.':
        return ''
    parts = []
    for p in s.split(';'):
        p = p.strip()
        if not p.lower().startswith('evalue='):
            parts.append(p)
    return ';'.join(parts)

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

def normalize_type(s):
    return re.sub(r"[.\-'/]", '', s.lower())

def gen_uid(seq, typ, cnt, suf="", basename=None):
    if basename:
        if not seq.startswith(f"{basename}_"):
            new_seq = f"{basename}_{seq}"
        else:
            new_seq = seq
        return f"{new_seq}.{normalize_type(typ)}{cnt}{suf}"
    else:
        return f"{seq}_{normalize_type(typ)}{cnt}{suf}"

# ---------------- Mitos-specific ----------------

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
        
        # Replace "ncRNA_gene" with "gene"
        if parts[2] == "ncRNA_gene":
            parts[2] = "gene"
        
        # remove parentheses from attributes
        parts[8] = fix_mitos_attr_field(parts[8], parts[0], args)
        fixed.append("\t".join(parts))
    return fixed

def fix_mitos_attr_field(attr, seq, args):
    attr = re.sub(r"\(.*?\)", "", attr).strip()
    new_attrs = []
    for chunk in attr.split(';'):
        chunk = chunk.strip()
        if not chunk:
            continue
        if '=' in chunk:
            key, value = chunk.split('=', 1)
            if key in ("ID", "Parent", "gene_id", "Name"):
                if not value.startswith(f"{seq}."):
                    value = f"{seq}.{value}"
            new_attrs.append(f"{key}={value}")
        else:
            new_attrs.append(chunk)
    out = ';'.join(new_attrs)
    if out and not out.endswith(';'):
        out += ';'
    return out

# ---------------- Infernal Overlap Filtering ----------------

def filter_overlapping_hits(lines, args):
    logging.info("Filtering overlapping Infernal hits...")
    hits = []
    for line in lines:
        if line.startswith("#"):
            continue
        fs = line.split()
        try:
            if args.fmt2:
                query = fs[3] if args.cmscan else fs[1]
                s, e = int(fs[9]), int(fs[10])
                strand, bit, ev = fs[11], float(fs[16]), float(fs[17])
            else:
                query = fs[2] if args.cmscan else fs[0]
                s, e = int(fs[7]), int(fs[8])
                strand, bit, ev = fs[9], float(fs[14]), float(fs[15])
        except Exception as ex:
            logging.error(f"Parse error: {line} => {ex}")
            continue
        start, end = (s, e) if strand == "+" else (e, s)
        hits.append({
            "query": query,
            "strand": strand,
            "start": start,
            "end": end,
            "bit": bit,
            "evalue": ev,
            "line": line
        })
    from collections import defaultdict
    groups = defaultdict(list)
    for h in hits:
        groups[(h["query"], h["strand"])].append(h)
    
    filt = []
    dropped = []
    for (key, hs) in groups.items():
        hs.sort(key=lambda x: x["start"])
        cluster = []
        cur_end = None
        for h in hs:
            if not cluster:
                cluster.append(h)
                cur_end = h["end"]
            elif h["start"] <= cur_end:
                cluster.append(h)
                cur_end = max(cur_end, h["end"])
            else:
                chosen = min(cluster, key=lambda x: (x["evalue"], -x["bit"]))
                for h_drop in cluster:
                    if h_drop is not chosen:
                        dropped.append(
                            (h_drop["line"], f"Overlap with better-scoring: {chosen['line'].strip()}")
                        )
                filt.append(chosen)
                cluster = [h]
                cur_end = h["end"]
        if cluster:
            chosen = min(cluster, key=lambda x: (x["evalue"], -x["bit"]))
            for h_drop in cluster:
                if h_drop is not chosen:
                    dropped.append(
                        (h_drop["line"], f"Overlap with better-scoring: {chosen['line'].strip()}")
                    )
            filt.append(chosen)
    
    logging.info(f"Retained {len(filt)} hits after overlap filtering, dropped {len(dropped)}.")
    return [h["line"] for h in filt], dropped

# ---------------- Infernal Tblout → GFF ----------------

def convert_tblout_to_gff(lines, args):
    logging.info("Converting Infernal tblout to GFF...")
    accepted = []
    dropped = []
    src = "cmscan" if args.cmscan else "cmsearch"
    if args.version:
        src += "-" + args.version
    if args.source:
        src = args.source
    
    for line in lines:
        if line.startswith("#"):
            continue
        fs = line.rstrip("\n").split()
        if args.fmt2:
            if len(fs) < 27:
                dropped.append((line, "Insufficient fields (>=27 needed for --fmt2)"))
                continue
            seq = fs[3] if args.cmscan else fs[1]
            s_from, s_to = fs[9], fs[10]
            strand, score, ev = fs[11], fs[16], fs[17]
            feature = fs[1] if args.cmscan else fs[2]
        else:
            if len(fs) < 18:
                dropped.append((line, "Insufficient fields (>=18 needed)"))
                continue
            seq = fs[2] if args.cmscan else fs[0]
            s_from, s_to = fs[7], fs[8]
            strand, score, ev = fs[9], fs[14], fs[15]
            feature = fs[2] if args.cmscan else fs[2]
        
        if args.ignore_trna and feature.lower() == "trna":
            dropped.append((line, "Ignored tRNA (due to --ignore-trna)"))
            continue
        
        try:
            s, e = (int(s_from), int(s_to)) if strand == "+" else (int(s_to), int(s_from))
            score_f = f"{float(score):.1f}"
        except Exception as ex:
            dropped.append((line, f"coord parse error: {ex}"))
            continue
        
        if args.min_bit is not None and float(score) < args.min_bit:
            dropped.append((line, f"Bit score {score} < threshold {args.min_bit}"))
            continue
        if args.max_evalue is not None and float(ev) > args.max_evalue:
            dropped.append((line, f"E-value {ev} > threshold {args.max_evalue}"))
            continue
        
        attr_field = "." if args.none else f"evalue={ev}"
        
        accepted.append("\t".join([
            seq, src, feature, str(s), str(e), score_f, strand, ".", attr_field
        ]))
    
    logging.info(f"Converted {len(accepted)} lines to GFF; dropped {len(dropped)}.")
    return accepted, dropped

# ---------------- Braker3 fix ----------------

def fix_braker3(lines, args):
    logging.info("Fixing Braker3 lines: add gene_biotype=protein_coding to gene feats.")
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

# ---------------- GFF fix for infernal / tRNAscan-se ----------------

def fix_gff_lines(lines, csv_fp, in_fmt, args):
    csvdata = load_csv(csv_fp)
    logging.info(f"Loaded CSV with {len(csvdata)} rows.")
    
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
            # if 'typ' in csv => treat as Rfam
            if typ in csvdata:
                info = csvdata[typ]
                toks = info['Type']
                if any(x.lower() == "gene" for x in toks):
                    child_tok = toks[1] if len(toks) > 1 else toks[0]
                    cnts[(seq, child_tok.lower()+"_g")] += 1
                    gid = gen_uid(seq, child_tok+"_g", cnts[(seq, child_tok.lower()+"_g")], basename=args.basename)
                    pattrs = {
                        'ID': gid,
                        'gene_biotype': child_tok,
                        'description': info['Description']
                    }
                    out.append("\t".join([
                        seq, src, "gene", start, end, score, strand, phase,
                        "." if args.none else reconst_attrs(pattrs)
                    ]))
                    
                    cnts[(seq, child_tok.lower())] += 1
                    cid = gen_uid(seq, child_tok, cnts[(seq, child_tok.lower())], basename=args.basename)
                    cattrs = {'ID': cid, 'Parent': gid}
                    out.append("\t".join([
                        seq, src, child_tok, start, end, score, strand, phase,
                        "." if args.none else reconst_attrs(cattrs)
                    ]))
                    
                    exon_id = f"{cid}.exon1"
                    exon_attrs = {'ID': exon_id, 'Parent': cid}
                    out.append("\t".join([
                        seq, src, "exon", start, end, score, strand, phase,
                        "." if args.none else reconst_attrs(exon_attrs)
                    ]))
                else:
                    # treat as "biological_region"
                    cnts[(seq, "bio_reg")] += 1
                    bid = gen_uid(seq, "biological_region", cnts[(seq, "bio_reg")], basename=args.basename)
                    battrs = {
                        'ID': bid,
                        'description': info['Description']
                    }
                    out.append("\t".join([
                        seq, src, "biological_region", start, end, score, strand, phase,
                        "." if args.none else reconst_attrs(battrs)
                    ]))
            else:
                # not in CSV => keep, remove evalue
                out.append("\t".join(parts[:8] + [
                    "." if args.none else remove_evalue(attr) or "."
                ]))
        
        elif in_fmt == "trnascan-se":
            if typ.lower() in ("trna", "pseudogene"):
                tcnts[(seq, "gene_tRNA")] += 1
                gnum = tcnts[(seq, "gene_tRNA")]
                gid = gen_uid(seq, "trna_g", gnum, basename=args.basename)
                gattrs = {'ID': gid, 'gene_biotype': "tRNA"}
                out.append("\t".join([
                    seq, src, "gene", start, end, score, strand, phase,
                    "." if args.none else reconst_attrs(gattrs)
                ]))
                
                tcnts[(seq, "child_tRNA")] += 1
                cnum = tcnts[(seq, "child_tRNA")]
                cid = gen_uid(seq, "trna", cnum, basename=args.basename)
                d = parse_attrs(remove_evalue(attr))
                d['ID'] = cid
                d['Parent'] = gid
                out.append("\t".join([
                    seq, src, "tRNA", start, end, score, strand, phase,
                    "." if args.none else reconst_attrs(d)
                ]))
            else:
                out.append("\t".join(parts[:8] + [
                    "." if args.none else remove_evalue(attr) or "."
                ]))
        else:
            out.append("\t".join(parts[:8] + [
                "." if args.none else remove_evalue(attr) or "."
            ]))
    # single version line at top
    final = []
    version_found = False
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

# ---------------- Cleanup: ensure 9 columns & single version line ----------------

def cleanup_for_fix(lines):
    out = []
    version_found = False
    for line in lines:
        if line.startswith("##gff-version"):
            if not version_found:
                out.append("##gff-version 3")
                version_found = True
            # skip additional version lines
            continue
        out.append(line)
    
    final = []
    for line in out:
        if line.startswith("#"):
            final.append(line)
            continue
        parts = line.split("\t")
        if len(parts) >= 9:
            # discard columns beyond the 9th
            final.append("\t".join(parts[:9]))
        else:
            # if < 9, keep as is
            final.append(line)
    return final

# ---------------- Sorting & Merging for "all" ----------------

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
        parts = parts[:9]
        rows.append(parts)
    
    if not rows:
        return ["##gff-version 3"] + headers
    
    df = pd.DataFrame(rows, columns=[
        'scaffold','source','feature','start','end','score','strand','phase','attributes'
    ])
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df = df.dropna(subset=['start'])
    df['hierarchy'] = df['attributes'].apply(lambda x: 1 if 'Parent=' in x else 0)
    df = df.sort_values(by=['scaffold','start','hierarchy'])
    df = df[['scaffold','source','feature','start','end','score','strand','phase','attributes']]
    merged = ["\t".join(map(str, row)) for row in df.values]
    
    # single version line + unique comments
    unique_hdr = []
    for h in headers:
        if h not in unique_hdr:
            unique_hdr.append(h)
    final = ["##gff-version 3"] + unique_hdr + merged
    return final

# ---------------- Single-file processing ----------------

def process_file(fmt, fname, args):
    with open(fname, 'r') as f:
        lines = f.readlines()
    logging.info(f"Processing file: {fname}, format={fmt}, {len(lines)} lines.")
    
    dropped = []
    fmt_lower = fmt.lower()
    if fmt_lower == "infernal":
        filtered, drops1 = filter_overlapping_hits(lines, args)
        converted, drops2 = convert_tblout_to_gff(filtered, args)
        lines = converted
        dropped = drops1 + drops2
    elif fmt_lower == "braker3":
        pass
    elif fmt_lower == "mitos":
        pass
    elif fmt_lower == "trnascan-se":
        pass
    
    fixed = fix_gff_lines_main(lines, args.csv, fmt_lower, args)
    return fixed, dropped

# ---------------- fix mode ----------------

def process_fix_cmd(args):
    outputs = {}
    drop_logs = {}
    for (fmt, fname) in args.I:
        gff_lines, dropped = process_file(fmt, fname, args)
        cleaned = cleanup_for_fix(gff_lines)
        outputs[fname] = cleaned
        drop_logs[fname] = dropped
    return outputs, drop_logs

# ---------------- all mode ----------------

def process_all(args):
    """
    In 'all' mode:
     - We fix each file, save an intermediate <root>_allfix.gff for the user
     - Also gather all lines into a list to produce merged_fix.gff
    """
    all_data = []
    all_drops = {}
    
    for (fmt, fname) in args.I:
        gff_lines, dropped = process_file(fmt, fname, args)
        cleaned = cleanup_for_fix(gff_lines)
        
        # Write the intermediate file for user inspection
        base = os.path.basename(fname)
        root, _ = os.path.splitext(base)
        out_fname = os.path.join(args.outdir, f"{root}_allfix.gff")
        with open(out_fname, 'w') as outf:
            outf.write("\n".join(cleaned) + "\n")
        logging.info(f"Wrote intermediate all-fix: {out_fname}")
        
        all_data.extend(cleaned)
        all_drops[fname] = dropped
    
    # Separate Mitos vs non-Mitos
    non_mitos, mitos = [], []
    for line in all_data:
        cols = line.strip().split("\t")
        if len(cols) < 2:
            non_mitos.append(line)
            continue
        if cols[1].lower() == "mitos":
            mitos.append(line)
        else:
            non_mitos.append(line)
    
    merged_non_mitos = sort_merge_gff_lines(non_mitos)
    # remove leading version lines from mitos
    mitos_clean = [m for m in mitos if not m.startswith("##gff-version")]
    merged = merged_non_mitos + mitos_clean
    return merged, all_drops

# ---------------- main ----------------

def main():
    parser = argparse.ArgumentParser(
        description="GFF combine & fix tool. In 'all' mode, also produce intermediate _allfix.gff files."
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
    fixp.add_argument("--cmscan", action="store_true", help="Infernal input from cmscan tblout.")
    fixp.add_argument("--fmt2", action="store_true", help="Infernal tblout with --fmt2 (>=27 fields).")
    grp = fixp.add_mutually_exclusive_group()
    grp.add_argument("--all-attrs", action="store_true", help="(Unused) keep all attrs.")
    grp.add_argument("--none", action="store_true", help="Use '.' for attributes.")
    grp.add_argument("--desc", action="store_true", help="(Unused) store description only.")
    fixp.add_argument("--version", type=str, default=None, help="Append version to source.")
    fixp.add_argument("--extra", type=str, default=None, help="(Unused) extra info in attributes.")
    fixp.add_argument("--hidedesc", action="store_true", help="(Unused) hide 'desc:'.")
    fixp.add_argument("--source", type=str, default=None, help="Override source field (infernal).")
    fixp.add_argument("--ignore-trna", action="store_true", help="Skip tRNA hits in infernal.")
    fixp.add_argument("--basename", type=str, default=None,
                      help="Prepend to generated gene IDs, e.g. 'M_tai' => 'M_tai_contigXY.gene1'.")
    
    allp = sub.add_parser('all', help="Fix all => <root>_allfix.gff, then merge => merged_fix.gff.")
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
    allp.add_argument("--cmscan", action="store_true", help="Infernal input from cmscan tblout.")
    allp.add_argument("--fmt2", action="store_true", help="Infernal tblout with --fmt2 (>=27 fields).")
    grp2 = allp.add_mutually_exclusive_group()
    grp2.add_argument("--all-attrs", action="store_true", help="(Unused) keep all attrs.")
    grp2.add_argument("--none", action="store_true", help="Use '.' for attributes.")
    grp2.add_argument("--desc", action="store_true", help="(Unused) store description only.")
    allp.add_argument("--version", type=str, default=None, help="Append version to source.")
    allp.add_argument("--extra", type=str, default=None, help="(Unused) extra info in attributes.")
    allp.add_argument("--hidedesc", action="store_true", help="(Unused) hide 'desc:'.")
    allp.add_argument("--source", type=str, default=None, help="Override source field (infernal).")
    allp.add_argument("--ignore-trna", action="store_true", help="Skip tRNA hits in infernal.")
    allp.add_argument("--basename", type=str, default=None,
                      help="Prepend to generated gene IDs, e.g. 'M_tai' => 'M_tai_contigXY.gene1'.")
    
    args = parser.parse_args()
    
    if args.command == "fix":
        outputs, drop_logs = process_fix_cmd(args)
        for in_fname, lines in outputs.items():
            base = os.path.basename(in_fname)
            root, _ = os.path.splitext(base)
            out_fname = os.path.join(args.outdir, f"{root}_fix.gff")
            with open(out_fname, 'w') as outF:
                outF.write("\n".join(lines) + "\n")
            logging.info(f"Wrote: {out_fname}")
            
            # If infernal => drop log
            for (fmt, fn) in args.I:
                if fmt.lower() == "infernal" and fn == in_fname:
                    drop_prefix = args.save_filtered_hits if args.save_filtered_hits else f"{root}_dropout"
                    drop_fname = os.path.join(args.outdir, f"{drop_prefix}.log")
                    with open(drop_fname, 'w') as dF:
                        if drop_logs[in_fname]:
                            for ln, reason in drop_logs[in_fname]:
                                dF.write(f"{ln.rstrip()}  # Dropped: {reason}\n")
                        else:
                            dF.write("No dropped hits.\n")
                    logging.info(f"Wrote drop log: {drop_fname}")
    
    elif args.command == "all":
        if len(args.I) < 2:
            parser.error("The 'all' command requires at least two input files.")
        merged, all_drops = process_all(args)
        
        # Write merged
        merged_fname = os.path.join(args.outdir, "merged_fix.gff")
        with open(merged_fname, 'w') as outF:
            outF.write("\n".join(merged) + "\n")
        logging.info(f"Wrote merged file: {merged_fname}")
        
        # Write drop logs
        for (fmt, fn) in args.I:
            if fmt.lower() == "infernal":
                base = os.path.basename(fn)
                root, _ = os.path.splitext(base)
                drop_prefix = args.save_filtered_hits if args.save_filtered_hits else f"{root}_dropout"
                drop_fname = os.path.join(args.outdir, f"{drop_prefix}.log")
                drops = all_drops[fn]
                with open(drop_fname, 'w') as dF:
                    if drops:
                        for ln, reason in drops:
                            dF.write(f"{ln.rstrip()}  # Dropped: {reason}\n")
                    else:
                        dF.write("No dropped hits.\n")
                logging.info(f"Wrote drop log: {drop_fname}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
