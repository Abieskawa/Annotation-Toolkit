#!/usr/bin/env python3
"""
Combined GFF tool:
  - Fixes Infernal/tRNAscan-SE GFF to create gene/child/exon hierarchies.
  - For Infernal tblout input, filters overlapping hits then converts to GFF.
    * First, overlapping hits are reduced (with dropped ones logged).
    * Then, if --min-bit or --max-evalue are provided, hits with a bit score below
      the --min-bit threshold or an E-value above the --max-evalue threshold are dropped.
      Dropped hits (with reasons) are logged if --save-filtered-hits is specified
      (this option applies only in infernal mode).
  - For tRNA hits:
    * If --ignore-trna is provided, tRNA hits are dropped (this option applies only to Infernal input).
    * Otherwise, a parent gene feature and a child tRNA feature are created.
  - For Braker3, gene features get gene_biotype=protein_coding in the attribute column.
  
Note:
  For Infernal output, we assume there is only one exon per hit.
  
Default CSV: Rfam_15_0.csv
Usage: -I FORMAT FILE (supported: infernal, trnascan-se, braker3); --step {fix, sort, merge}
"""

import sys, argparse, csv, logging, re, os
from collections import defaultdict
import pandas as pd

# --- Utility functions ---
def load_csv(fp):
    # Get the directory where the script is stored.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # If fp is not an absolute path, assume it is relative to the script's directory.
    if not os.path.isabs(fp):
        fp_in_script = os.path.join(script_dir, fp)
        if os.path.exists(fp_in_script):
            fp = fp_in_script
    # Otherwise, try the current working directory.
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
            d[key] = {'Accession': row['Accession'].strip(),
                      'Type': [t.strip() for t in row['Type'].split(';')],
                      'Description': row['Description'].strip()}
    return d

def remove_evalue(s):
    return ';'.join(p.strip() for p in s.split(';') if not p.lower().startswith('evalue=')) if s and s!='.' else ''

def parse_attrs(a):
    d = {}
    if not a or a=='.': 
        return d
    for p in a.split(';'):
        if '=' in p:
            k, v = p.split('=', 1)
            d[k.strip()] = v.strip()
    return d

def reconst_attrs(d):
    return '.' if not d else ';'.join(f"{k}={v}" for k,v in d.items()) + ';'

def normalize_type(s):
    # Remove dots, dashes, slashes, and apostrophes but preserve underscores.
    return re.sub(r"[.\-'/]", '', s.lower())

def gen_uid(seq, typ, cnt, suf="", basename=None):
    # If a basename is provided, check whether the seq already starts with it.
    if basename:
        if not seq.startswith(f"{basename}_"):
            new_seq = f"{basename}_{seq}"
        else:
            new_seq = seq
        return f"{new_seq}.{normalize_type(typ)}{cnt}{suf}"
    else:
        return f"{seq}_{normalize_type(typ)}{cnt}{suf}"

# --- Overlap filtering for Infernal tblout ---
def filter_overlapping_hits(lines, args):
    logging.info("Filtering overlapping hits...")
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
            logging.error(f"Parse error: {line} - {ex}")
            continue
        start, end = (s, e) if strand == "+" else (e, s)
        hits.append({"query": query, "strand": strand, "start": start, "end": end,
                     "evalue": ev, "bit": bit, "line": line})
    groups = defaultdict(list)
    for h in hits:
        groups[(h["query"], h["strand"])].append(h)
    filt = []
    dropped = []
    for key, hs in groups.items():
        hs.sort(key=lambda h: h["start"])
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
                chosen = min(cluster, key=lambda h: (h["evalue"], -h["bit"]))
                for h_drop in cluster:
                    if h_drop is not chosen:
                        dropped.append((h_drop["line"], f"Overlaps with higher scoring hit selected: {chosen['line'].strip()}"))
                filt.append(chosen)
                cluster = [h]
                cur_end = h["end"]
        if cluster:
            chosen = min(cluster, key=lambda h: (h["evalue"], -h["bit"]))
            for h_drop in cluster:
                if h_drop is not chosen:
                    dropped.append((h_drop["line"], f"Overlaps with higher scoring hit selected: {chosen['line'].strip()}"))
            filt.append(chosen)
    logging.info(f"Retained {len(filt)} hits after overlapping filtering.")
    logging.info(f"Dropped {len(dropped)} overlapping hits due to overlaps.")
    if dropped:
        for dline, reason in dropped:
            logging.debug(f"Dropped hit: {dline.strip()}  # Reason: {reason}")
    else:
        logging.debug("No overlapping hits were dropped.")
    return [h["line"] for h in filt], dropped

# --- Tblout conversion with filtering and logging dropped lines ---
def convert_tblout_to_gff(lines, args):
    logging.info("Converting tblout to GFF...")
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
                logging.error(f"Insufficient fields: {line}")
                dropped.append((line, "Insufficient fields (expected >=27 for --fmt2)"))
                continue
            seq = fs[3] if args.cmscan else fs[1]
            s_from, s_to = fs[9], fs[10]
            strand, score, ev = fs[11], fs[16], fs[17]
            feature = fs[1] if args.cmscan else fs[2]
        else:
            if len(fs) < 18:
                logging.error(f"Insufficient fields: {line}")
                dropped.append((line, "Insufficient fields (expected >=18)"))
                continue
            seq = fs[2] if args.cmscan else fs[0]
            s_from, s_to = fs[7], fs[8]
            strand, score, ev = fs[9], fs[14], fs[15]
            feature = fs[2] if args.cmscan else fs[2]
        # For Infernal input, if --ignore-trna is set, drop tRNA hits.
        if args.ignore_trna and args.I[0].lower() == "infernal" and feature.lower() == "trna":
            dropped.append((line, "Ignored tRNA hit due to --ignore-trna"))
            continue
        try:
            s, e = (int(s_from), int(s_to)) if strand == "+" else (int(s_to), int(s_from))
            score_f = f"{float(score):.1f}"
        except Exception as ex:
            logging.error(f"Conversion error: {line} - {ex}")
            dropped.append((line, f"Conversion error: {ex}"))
            continue
        if args.min_bit is not None and float(score) < args.min_bit:
            dropped.append((line, f"Bit score {score} below threshold {args.min_bit}"))
            continue
        if args.max_evalue is not None and float(ev) > args.max_evalue:
            dropped.append((line, f"E-value {ev} above threshold {args.max_evalue}"))
            continue
        attr_field = "." if args.none else f"evalue={ev}"
        accepted.append("\t".join([seq, src, feature, str(s), str(e), score_f, strand, ".", attr_field]))
    logging.info(f"Converted {len(accepted)} tblout lines to GFF. Dropped {len(dropped)} lines due to filtering.")
    return accepted, dropped

# --- GFF Fixing for Infernal and trnascan‑se ---
def fix_gff_lines(lines, csv_fp, in_fmt, args):
    csvdata = load_csv(csv_fp)
    # Log detailed CSV info to help the user understand what was loaded.
    sample_keys = list(csvdata.keys())[:5]
    logging.info(f"CSV file '{csv_fp}' loaded with {len(csvdata)} entries. Example entry IDs: {sample_keys}")
    if in_fmt == "infernal":
        logging.info("Assumption: Each Infernal hit is assumed to have only one exon.")
    out = []
    cnts, tcnts = defaultdict(int), defaultdict(int)
    for line in lines:
        if line.startswith("#") or not line.strip():
            out.append(line.rstrip("\n"))
            continue
        parts = line.rstrip("\n").split('\t')
        if len(parts) != 9:
            logging.warning(f"Malformed: {line}")
            continue
        seq, src, typ, start, end, score, strand, phase, attr = parts

        # Case 1. Infernal input: build hierarchy using CSV info.
        if in_fmt == "infernal" and typ in csvdata:
            info = csvdata[typ]
            toks = [t.strip() for t in info['Type']]
            if any(t.lower() == "gene" for t in toks):
                # Note: We assume there is only one exon per Infernal hit.
                parent_tok = toks[0]
                child_tok = toks[1] if len(toks) > 1 else toks[0]
                cnts[(seq, child_tok.lower()+"_g")] += 1
                gid = gen_uid(seq, child_tok + "_g", cnts[(seq, child_tok.lower()+"_g")], basename=args.basename)
                pattrs = {'ID': gid, 'gene_biotype': child_tok, 'description': info['Description']}
                out.append("\t".join([seq, src, "gene", start, end, score, strand, phase,
                                       ("." if args.none else reconst_attrs(pattrs))]))
                cnts[(seq, child_tok.lower())] += 1
                cid = gen_uid(seq, child_tok, cnts[(seq, child_tok.lower())], basename=args.basename)
                cattrs = {'ID': cid, 'Parent': gid}
                out.append("\t".join([seq, src, child_tok, start, end, score, strand, phase,
                                       ("." if args.none else reconst_attrs(cattrs))]))
                exon_id = f"{cid}.exon1"
                exon_attrs = {'ID': exon_id, 'Parent': cid}
                out.append("\t".join([seq, src, "exon", start, end, score, strand, phase,
                                       ("." if args.none else reconst_attrs(exon_attrs))]))
            else:
                cnts[(seq, "bio_reg")] += 1
                bid = gen_uid(seq, "biological_region", cnts[(seq, "bio_reg")], basename=args.basename)
                pattrs = {'ID': bid, 'description': info['Description']}
                out.append("\t".join([seq, src, "biological_region", start, end, score, strand, phase,
                                       ("." if args.none else reconst_attrs(pattrs))]))
        # Case 2. trnascan-se input: process tRNAscan-SE hits (for types "trna" or "pseudogene")
        elif in_fmt == "trnascan-se" and typ.lower() in ("trna", "pseudogene"):
            tcnts[(seq, "child_tRNA")] += 1
            child_num = tcnts[(seq, "child_tRNA")]
            tcnts[(seq, "gene_tRNA")] += 1
            gene_num = tcnts[(seq, "gene_tRNA")]
            gid = gen_uid(seq, "trna_g", gene_num, basename=args.basename)
            gattrs = {'ID': gid, 'gene_biotype': "tRNA"}
            out.append("\t".join([seq, src, "gene", start, end, score, strand, phase,
                                   ("." if args.none else reconst_attrs(gattrs))]))
            child_id = gen_uid(seq, "trna", child_num, basename=args.basename)
            d = parse_attrs(remove_evalue(attr))
            d['ID'] = child_id
            d['Parent'] = gid
            out.append("\t".join([seq, src, "trna", start, end, score, strand, phase,
                                   ("." if args.none else reconst_attrs(d))]))
        # Case 3. Fallback: for any feature with type "trna" (non‑trnascan‑se)
        elif typ.lower() == "trna":
            if in_fmt == "infernal" and args.ignore_trna:
                continue
            tcnts[(seq, "child_tRNA")] += 1
            child_num = tcnts[(seq, "child_tRNA")]
            tcnts[(seq, "gene_tRNA")] += 1
            gene_num = tcnts[(seq, "gene_tRNA")]
            gid = gen_uid(seq, "trna_g", gene_num, basename=args.basename)
            gattrs = {'ID': gid, 'gene_biotype': "tRNA"}
            out.append("\t".join([seq, src, "gene", start, end, score, strand, phase,
                                   ("." if args.none else reconst_attrs(gattrs))]))
            child_id = gen_uid(seq, "trna", child_num, basename=args.basename)
            d = parse_attrs(remove_evalue(attr))
            d['ID'] = child_id
            d['Parent'] = gid
            out.append("\t".join([seq, src, typ, start, end, score, strand, phase,
                                   ("." if args.none else reconst_attrs(d))]))
        else:
            # For all other cases, simply remove any evalue from the attributes (unless --none is set).
            out.append("\t".join(parts[:-1] + [("." if args.none else remove_evalue(attr) or ".")]))
    if not any(l.startswith("##gff-version 3") for l in out):
        out.insert(0, "##gff-version 3")
    logging.info(f"Fixed GFF entries: {len(out)} lines.")
    return out

# --- GFF Fixing for Braker3 ---
def fix_braker3(lines, args):
    logging.info("Processing Braker3 input: adding gene_biotype=protein_coding for gene features and removing empty lines.")
    new_lines = []
    for line in lines:
        # Skip blank lines or lines with only whitespace.
        if not line.strip():
            continue
        if line.startswith("#"):
            new_lines.append(line.rstrip("\n"))
            continue
        parts = line.rstrip("\n").split('\t')
        if len(parts) != 9:
            new_lines.append(line.rstrip("\n"))
            continue
        seq, src, feat, start, end, score, strand, phase, attr = parts
        if feat.lower() == "gene":
            if attr == ".":
                attr = "gene_biotype=protein_coding;"
            else:
                if "gene_biotype=" not in attr:
                    if attr.endswith(";"):
                        attr = attr + "gene_biotype=protein_coding;"
                    else:
                        attr = attr + ";gene_biotype=protein_coding;"
        new_line = "\t".join([seq, src, feat, start, end, score, strand, phase, attr])
        new_lines.append(new_line)
    # Remove any possible empty lines
    new_lines = [line for line in new_lines if line.strip() != ""]
    return new_lines

# --- Main Processing ---
def process_fix(args):
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    if args.I:
        fmt = args.I[0].lower()
        with open(args.I[1], 'r') as f:
            lines = f.readlines()
        total_input_lines = len(lines)
        in_fmt = fmt
        logging.info(f"Input format: {fmt} | {total_input_lines} lines read.")
    else:
        lines = sys.stdin.readlines()
        total_input_lines = len(lines)
        in_fmt = "gff"
        logging.info(f"STDIN read: {total_input_lines} lines; assuming GFF format.")

    if in_fmt == "infernal":
        filtered_lines, overlapped_dropped = filter_overlapping_hits(lines, args)
        logging.info(f"After overlapping filtering: {len(filtered_lines)} lines remain (from original {total_input_lines}).")
        accepted, conv_dropped = convert_tblout_to_gff(filtered_lines, args)
        total_dropped = overlapped_dropped + conv_dropped
        logging.info(f"After score filtering: {len(accepted)} lines accepted; {len(total_dropped)} lines dropped.")
        if args.save_filtered_hits:
            with open(args.save_filtered_hits, 'w') as f:
                f.write(f"Total input lines: {total_input_lines}\n")
                f.write(f"Lines after overlapping filtering: {len(filtered_lines)}\n")
                for ln, reason in total_dropped:
                    f.write(f"{ln.rstrip()}  # Dropped: {reason}\n")
                f.write(f"\nAccepted lines: {len(accepted)}\nDropped lines: {len(total_dropped)}\n")
            logging.info(f"Filtered tblout saved to {args.save_filtered_hits}")
        gff_lines = accepted
    elif in_fmt == "braker3":
        # For braker3, process with fix_braker3
        gff_lines = fix_braker3(lines, args)
    else:
        # For non-Infernal inputs (including trnascan-se), assume the file is already in GFF.
        gff_lines = lines
    logging.info(f"GFF lines count after conversion: {len(gff_lines)}")
    
    if args.step == "fix":
        if in_fmt == "braker3":
            fixed = gff_lines
        else:
            fixed = fix_gff_lines(gff_lines, args.csv, in_fmt, args)
        # Remove any stray empty lines.
        fixed = [line for line in fixed if line.strip() != ""]
        logging.info(f"Fixed GFF contains {len(fixed)} lines.")
        if args.output:
            with open(args.output, 'w') as f:
                f.write("\n".join(fixed))
        else:
            print("\n".join(fixed))
    elif args.step == "sort":
        fixed = gff_lines
        if not any(l.startswith("##gff-version 3") for l in fixed):
            fixed.insert(0, "##gff-version 3")
        hdr = [l for l in fixed if l.startswith("#")]
        data = [l for l in fixed if not l.startswith("#")]
        rows = [l.split('\t') for l in data if len(l.split('\t')) == 9]
        df = pd.DataFrame(rows, columns=['scaffold','source','feature','start','end','score','strand','phase','attributes'])
        df['start'] = pd.to_numeric(df['start'])
        df['hierarchy'] = df['attributes'].apply(lambda a: 1 if "Parent=" in a else 0)
        sorted_df = df.sort_values(by=['scaffold','start','hierarchy'])
        outl = hdr + ["\t".join(map(str, row)) for row in sorted_df.values]
        logging.info(f"Sorted GFF contains {len(outl)} lines.")
        if args.output:
            with open(args.output, 'w') as f:
                f.write("\n".join(outl))
        else:
            print("\n".join(outl))
    elif args.step == "merge":
        fixed = gff_lines if in_fmt == "braker3" else fix_gff_lines(gff_lines, args.csv, in_fmt, args)
        if not any(l.startswith("##gff-version 3") for l in fixed):
            fixed.insert(0, "##gff-version 3")
        hdr = [l for l in fixed if l.startswith("#")]
        data = [l for l in fixed if not l.startswith("#")]
        rows = [l.split('\t') for l in data if len(l.split('\t')) == 9]
        df = pd.DataFrame(rows, columns=['scaffold','source','feature','start','end','score','strand','phase','attributes'])
        df['start'] = pd.to_numeric(df['start'])
        df['hierarchy'] = df['attributes'].apply(lambda a: 1 if "Parent=" in a else 0)
        sorted_df = df.sort_values(by=['scaffold','start','hierarchy'])
        merged = hdr + ["\t".join(map(str, row)) for row in sorted_df.values]
        logging.info(f"Merged GFF contains {len(merged)} lines.")
        if args.output:
            with open(args.output, 'w') as f:
                f.write("\n".join(merged))
        else:
            print("\n".join(merged))

def main():
    parser = argparse.ArgumentParser(
        description="GFF tool: fixes GFF (with tblout filtering/conversion for infernal and trnascan-se), adds gene_biotype for Braker3, and sorts output. "
                    "Input is specified using -I FORMAT FILE. "
                    "For tblout (infernal) inputs, overlapping hits are first filtered then additional filtering based on --min-bit "
                    "and/or --max-evalue is applied (dropped hits with reasons are logged if --save-filtered-hits is provided [applies only in infernal mode]). "
                    "For tRNA hits, if --ignore-trna is set they are dropped (this applies only to Infernal input); otherwise, a parent gene and child tRNA feature are created."
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    fixp = subparsers.add_parser('fix', help="Run processing; --step: fix, sort, or merge.")
    fixp.add_argument("-I", nargs=2, metavar=("FORMAT","FILE"),
                      help="Input format and file. Supported: infernal, trnascan-se, braker3.")
    fixp.add_argument("-c", "--csv", type=str, default="Rfam_15_0.csv", 
                      help="CSV file (default: Rfam_15_0.csv). If a relative path is provided, the file is searched for in the script's directory by default.")
    fixp.add_argument("-o", "--output", type=str, default=None, help="Output file (or stdout).")
    fixp.add_argument("--step", choices=["fix", "sort", "merge"], default="merge",
                      help="Step: fix only, sort only, or full merge (fix then sort).")
    fixp.add_argument("--save-filtered-hits", type=str, default=None,
                      help="File to save dropped tblout hits with reasons (applies only in infernal mode).")
    fixp.add_argument("--min-bit", "-T", type=float, default=None, 
                      help="Minimum bit score threshold for tblout hits. Hits with a bit score below this threshold will be dropped after overlapping filtering. "
                           "Dropped hits (with reasons) are logged if --save-filtered-hits is specified.")
    fixp.add_argument("--max-evalue", "-E", type=float, default=None, 
                      help="Maximum E-value threshold for tblout hits. Hits with an E-value above this threshold will be dropped after overlapping filtering. "
                           "Dropped hits (with reasons) are logged if --save-filtered-hits is specified.")
    fixp.add_argument("--cmscan", action="store_true", help="Tblout from cmscan.")
    fixp.add_argument("--fmt2", action="store_true", help="Tblout with --fmt2 (27+ fields).")
    grp = fixp.add_mutually_exclusive_group()
    grp.add_argument("--all", action="store_true", help="Output all attributes (tblout).")
    grp.add_argument("--none", action="store_true", help="Output no attributes (tblout): attribute field will be '.'")
    grp.add_argument("--desc", action="store_true", help="Output description instead of E-value (tblout).")
    fixp.add_argument("--version", type=str, default=None, help="Append version to source (tblout).")
    fixp.add_argument("--extra", type=str, default=None, help="Append extra info to attributes (tblout).")
    fixp.add_argument("--hidedesc", action="store_true", help="Do not prefix description with 'desc:' (tblout).")
    fixp.add_argument("--source", type=str, default=None, help="Specify source field (tblout).")
    fixp.add_argument("--ignore-trna", action="store_true", help="Ignore (drop) tRNA hits (applies only to Infernal input).")
    fixp.add_argument("--basename", type=str, default=None,
                      help="A basename to prepend to all generated IDs. For example, if provided as 'ABC', then gene IDs will be formatted as 'ABC_<scaffold>.<type>_g<counter>', e.g. ABC_chrom1.rrna_g2. If the input sequence already begins with the basename, the prefix is not added again.")
    args = parser.parse_args()
    if args.command == "fix":
        if args.I and args.I[0].lower() == "infernal":
            if args.min_bit is not None and args.max_evalue is not None:
                parser.error("ERROR: --min-bit and --max-evalue cannot be used together.")
            if args.fmt2 and not args.cmscan:
                parser.error("ERROR: --fmt2 only makes sense with --cmscan.")
            if args.source and args.version:
                parser.error("ERROR: --source and --version are incompatible.")
        process_fix(args)

if __name__ == "__main__":
    main()
