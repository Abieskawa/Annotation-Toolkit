#!/usr/bin/env python3
"""
Combined GFF tool:  
  - Fixes GFF from various sources (infernal, trnascan‑se, braker3, mitos, gb/gbf (mitohifi, mitoz)) to create 
    gene/child/exon hierarchies and ensures only 9 columns with a single "##gff-version 3".
  - Processes all input files (provided via -I) and writes intermediate <basename>_fix.gff files, 
    then merges them into one merged_fix.gff.
  - infernal2gff: we use bit score as the score value in gff file.
  
Assumption:
  - The exon number of ncRNA genes predicted with cmscan is only one (remove tRNA entries in this script)
  
References:
1. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
2. https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
3. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/
tblout2gff reference:
1. https://github.com/nawrockie/jiffy-infernal-hmmer-scripts/blob/master/infernal-tblout2gff.pl
"""

import sys, argparse, csv, logging, re, os, subprocess, tempfile
from collections import defaultdict
import pandas as pd

# --- Helper: Normalize IDs by stripping content in parentheses ---
def normalize_id(id_str):
    return re.sub(r'\(.*?\)', '', id_str)

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
    return ';'.join(f"{k}={v}" for k, v in d.items())

def normalize_type(s):
    return re.sub(r"[\-'/]", '', s.lower())

def gen_uid(seq, typ, cnt, basename=None):
    if basename:
        return f"{basename}_{seq}_{normalize_type(typ)}{cnt}"
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
            if key == "Is_circular":
                value = "true"
            if key in ("ID", "Parent", "gene_id"):
                if args.basename and not value.startswith(f"{args.basename}_"):
                    value = f"{args.basename}_{value}"
            new_attrs.append(f"{key}={value}")
        else:
            new_attrs.append(chunk)
    out = ';'.join(new_attrs)
    if out and not out.endswith(';'):
        out += ';'
    return out

def fix_mitos_lines(lines, args):
    header_lines = [l for l in lines if l.startswith("#")]
    data_lines = [l for l in lines if not l.startswith("#")]
    try:
        data_lines.sort(key=lambda l: int(l.split('\t')[3]) if len(l.split('\t')) >= 9 and l.split('\t')[3].isdigit() else 0)
    except Exception as e:
        logging.warning("Sorting mitos input failed: " + str(e))
    lines = header_lines + data_lines
    local_dropped = []
    inserted_mrna_ids = set()
    gene_info = {} 
    child_biotype = {}
    transcript_map = {}
    output_lines = []
    cds_counter = {}
    exon_counter = {}
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#"):
            output_lines.append(line)
            continue
        parts = line.split('\t')
        if len(parts) != 9:
            local_dropped.append((line, "[mitos] Incorrect number of columns"))
            continue
        parts[8] = fix_mitos_attr_field(parts[8], parts[0], args)
        feature = parts[2]
        if feature in ("gene", "ncRNA_gene"):
            orig = feature
            if feature == "ncRNA_gene":
                feature = "gene"
                parts[2] = "gene"
            attrs = parse_attrs(parts[8])
            gene_id = attrs.get("ID", "")
            if gene_id:
                gene_info[gene_id] = {
                    "seq": parts[0],
                    "source": parts[1],
                    "start": parts[3],
                    "end": parts[4],
                    "score": parts[5],
                    "strand": parts[6],
                    "phase": parts[7],
                    "attributes": parts[8],
                    "original": orig,
                    "final_id": gene_id
                }
            output_lines.append("\t".join(parts))
        elif feature == "exon":
            attrs = parse_attrs(parts[8])
            mrna_id = attrs.get("Parent", None)
            if mrna_id and (("transcript_oh" in mrna_id.lower()) or ("transcript_ol" in mrna_id.lower())):
                local_dropped.append((line, "[mitos] Parent transcript indicates origin-of-replication error"))
                continue
            if mrna_id:
                norm_mrna = normalize_id(mrna_id)
                if norm_mrna in transcript_map:
                    mrna_id = transcript_map[norm_mrna]
                if "transcript" in mrna_id:
                    computed_gene = mrna_id.replace("transcript", "gene", 1)
                else:
                    computed_gene = mrna_id
                norm_computed = normalize_id(computed_gene)
                for gid, rec in gene_info.items():
                    if normalize_id(rec["final_id"]) != norm_computed and normalize_id(gid) == norm_computed:
                        rec["final_id"] = computed_gene
                        break
                for gid, rec in gene_info.items():
                    if normalize_id(rec["final_id"]) == norm_computed and rec["original"] == "gene":
                        if mrna_id not in inserted_mrna_ids:
                            inserted_mrna_ids.add(mrna_id)
                            mrna_line = "\t".join([
                                rec["seq"],
                                rec["source"],
                                "mRNA",
                                rec["start"],
                                rec["end"],
                                rec["score"],
                                rec["strand"],
                                rec["phase"],
                                reconst_attrs({'ID': mrna_id, 'Parent': rec["final_id"]})
                            ])
                            output_lines.append(mrna_line)
                        break
                if mrna_id in inserted_mrna_ids:
                    if mrna_id not in cds_counter:
                        cds_counter[mrna_id] = 0
                    cds_counter[mrna_id] += 1
                    cds_id = f"{mrna_id}.cds{cds_counter[mrna_id]}"
                    cds_line = "\t".join([
                        parts[0],
                        parts[1],
                        "CDS",
                        parts[3],
                        parts[4],
                        ".",
                        parts[6],
                        "0",
                        reconst_attrs({'ID': cds_id, 'Parent': mrna_id})
                    ])
                    output_lines.append(cds_line)
                if parts[1] == "mitos":
                    if mrna_id not in exon_counter:
                        exon_counter[mrna_id] = 0
                    exon_counter[mrna_id] += 1
                    exon_id = f"{mrna_id}.exon{exon_counter[mrna_id]}"
                else:
                    name_attr = attrs.get("Name", "").strip()
                    if name_attr:
                        norm_name = re.sub(r'\s+', '_', name_attr)
                        exon_id = f"{mrna_id}.{norm_name}"
                    else:
                        exon_id = f"{mrna_id}.exon1"
                attrs["ID"] = exon_id
                attrs["Parent"] = mrna_id
                parts[8] = reconst_attrs(attrs)
            output_lines.append("\t".join(parts))
        else:
            attrs = parse_attrs(parts[8])
            parent_val = attrs.get("Parent", None)
            if parent_val:
                if "transcript" in parent_val:
                    computed_gene = parent_val.replace("transcript", "gene", 1)
                else:
                    computed_gene = parent_val
                norm_computed = normalize_id(computed_gene)
                for gid, rec in gene_info.items():
                    if normalize_id(rec["final_id"]) == norm_computed and rec["original"] == "ncRNA_gene":
                        if rec["final_id"] != computed_gene:
                            rec["final_id"] = computed_gene
                        if rec["final_id"] not in child_biotype:
                            child_biotype[rec["final_id"]] = feature
                        break
            if parts[2].lower() in ("trna", "rrna"):
                transcript_id = attrs.get("ID", None)
                if transcript_id:
                    transcript_map[normalize_id(transcript_id)] = transcript_id
            output_lines.append("\t".join(parts))
    updated_lines = []
    for line in output_lines:
        if line.startswith("#"):
            updated_lines.append(line)
            continue
        parts = line.split('\t')
        if len(parts) < 9:
            updated_lines.append(line)
            continue
        if parts[2] == "gene":
            attrs = parse_attrs(parts[8])
            orig_id = attrs.get("ID", "")
            for gid, rec in gene_info.items():
                if normalize_id(orig_id) == normalize_id(rec["final_id"]):
                    attrs["ID"] = rec["final_id"]
                    if "gene_biotype" not in attrs:
                        if rec["original"] == "gene":
                            attrs["gene_biotype"] = "protein_coding"
                        else:
                            attrs["gene_biotype"] = child_biotype.get(rec["final_id"], "ncRNA")
                    if "gene_id" in attrs:
                        attrs["gene_id"] = rec["final_id"]
                    break
            parts[8] = reconst_attrs(attrs)
            updated_lines.append("\t".join(parts))
        else:
            updated_lines.append(line)
    sorted_lines = sort_merge_gff_lines(updated_lines)
    return sorted_lines, local_dropped

def convert_tblout_to_gff(lines, args):
    accepted = []
    dropped = []
    if args.source:
        src = args.source
    else:
        print("Error: You must specify --source", file=sys.stderr)
        sys.exit(1)
    
    # Load CSV data for filtering
    csvdata = load_csv(args.csv)
    
    # Normalize ignore list - strip whitespace and convert to lowercase
    ignore_list = []
    if args.ignore:
        ignore_list = [item.strip().lower() for item in args.ignore.split(',')]
    
    for line in lines:
        if line.startswith("#"):
            continue
        fs = line.rstrip("\n").split()
        if args.fmt2:
            if len(fs) < 27:
                dropped.append((line, "[infernal] Insufficient fields"))
                continue
            seq = fs[3]
            s_from, s_to = fs[9], fs[10]
            strand, score, ev = fs[11], fs[16], fs[17]
            feature = fs[1]
        else:
            if len(fs) < 18:
                dropped.append((line, "[infernal] Insufficient fields"))
                continue
            seq = fs[2]
            s_from, s_to = fs[7], fs[8]
            strand, score, ev = fs[9], fs[14], fs[15]
            feature = fs[2]
        
        # No coordinate check is done – simply convert s_from and s_to.
        s = int(s_from)
        e = int(s_to)
        
        # CSV-based filtering logic
        should_ignore = False
        if feature in csvdata:
            # Get the Type information from CSV
            type_info = csvdata[feature]['Type']
            # Check if any ignore pattern matches the Type
            for ignore_pattern in ignore_list:
                for type_entry in type_info:
                    if ignore_pattern in type_entry.lower():
                        should_ignore = True
                        dropped.append((line, f"[infernal] Ignored based on CSV Type '{type_entry}' matching --ignore '{ignore_pattern}': {feature}"))
                        break
                if should_ignore:
                    break
        else:
            # If not in CSV, fall back to simple name matching (for backwards compatibility)
            feature_normalized = feature.strip().lower()
            for ignore_pattern in ignore_list:
                if ignore_pattern in feature_normalized:
                    should_ignore = True
                    dropped.append((line, f"[infernal] Ignored by name matching --ignore '{ignore_pattern}': {feature}"))
                    break
        
        if not should_ignore:
            accepted.append("\t".join([seq, src, feature, str(s), str(e), score, strand, ".", "."]))
    
    return accepted, dropped

def fix_braker3(lines, args):
    new_lines = []
    prefix = f"{args.basename}_" if args.basename else ""

    for raw in lines:
        line = raw.rstrip("\n")
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
        attrs = parse_attrs(attr) if attr and attr != "." else {}

        # 1) add default gene_biotype
        if feat.lower() == "gene" and "gene_biotype" not in attrs:
            attrs["gene_biotype"] = "protein_coding"

        # 2) add basename prefix where needed
        if prefix:
            for key in ("ID", "Parent", "gene_id", "transcript_id"):
                if key in attrs and not attrs[key].startswith(prefix):
                    attrs[key] = f"{prefix}{attrs[key]}"

        new_attr = reconst_attrs(attrs) if attrs else "."
        new_lines.append("\t".join(
            [seq, src, feat, start, end, score, strand, phase, new_attr]
        ))

    return new_lines

def fix_gff_lines(lines, csv_fp, in_fmt, args):
    csvdata = load_csv(csv_fp)
    out = []
    cnts, tcnts = defaultdict(int), defaultdict(int)
    
    if in_fmt == "trnascan-se":
        # First pass: build a set of pseudogene IDs (normalized to lower case)
        pseudogene_ids = set()
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split('\t')
            if len(parts) != 9:
                continue
            typ = parts[2].strip().lower()
            if typ == "pseudogene":
                d = parse_attrs(parts[8])
                pgid = d.get("ID", "").strip().lower()
                if args.basename and not pgid.startswith(f"{args.basename}"):
                    pgid = f"{args.basename}_{pgid}"
                pseudogene_ids.add(pgid)
        gene_type_map = {}
        for line in lines:
            if line.startswith("#") or not line.strip():
                out.append(line.rstrip("\n"))
                continue
            parts = line.rstrip("\n").split('\t')
            if len(parts) != 9:
                out.append(line.rstrip("\n"))
                continue
            seq, src, typ, start, end, score, strand, phase, attr = parts
            # If the line is a pseudogene, output it normally.
            if typ.lower() == "pseudogene":
                d = parse_attrs(attr)
                if args.basename:
                    if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                        d['ID'] = f"{args.basename}_{d['ID']}"
                    if 'Parent' in d and not d['Parent'].startswith(f"{args.basename}"):
                        d['Parent'] = f"{args.basename}_{d['Parent']}"
                gene_id = d.get("ID", "").strip().lower()
                gene_type_map[gene_id] = d.get("gene_biotype", typ).lower()
                out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
            # Else if the type exists in CSV (regular gene record)
            elif typ in csvdata:
                info = csvdata[typ]
                toks = info['Type']
                if any(x.lower() == "gene" for x in toks):
                    child_tok = toks[1] if len(toks) > 1 else typ
                    cnts[(seq, child_tok.lower()+"_g")] += 1
                    if args.basename:
                        gid = f"{args.basename}_{seq}_{normalize_type(child_tok + '_g')}{cnts[(seq, child_tok.lower() + '_g')]}"
                    else:
                        gid = gen_uid(seq, child_tok+"_g", cnts[(seq, child_tok.lower()+"_g")])
                    pattrs = {'ID': gid, 'gene_biotype': child_tok, 'description': info['Description']}
                    out.append("\t".join([seq, src, "gene", start, end, score, strand, phase, reconst_attrs(pattrs)]))
                    gene_type_map[gid] = child_tok.lower()
                    d = parse_attrs(attr)
                    d['Parent'] = gid
                    if args.basename:
                        if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                            d['ID'] = f"{args.basename}_{d['ID']}"
                    out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
                else:
                    d = parse_attrs(attr)
                    if args.basename:
                        if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                            d['ID'] = f"{args.basename}_{d['ID']}"
                        if 'Parent' in d and not d['Parent'].startswith(f"{args.basename}"):
                            d['Parent'] = f"{args.basename}_{d['Parent']}"
                    gene_id = d.get("ID", "").strip().lower()
                    gene_type_map[gene_id] = d.get("gene_biotype", typ).lower()
                    if typ.lower() == "pseudogene":
                        pseudogene_ids.add(gene_id)
                    out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
            else:
                if typ.lower() == "exon":
                    d = parse_attrs(attr)
                    parent_id = d.get("Parent", "").strip().lower()
                    if args.basename and not parent_id.startswith(f"{args.basename}"):
                        parent_id = f"{args.basename}_{parent_id}"
                    if parent_id in pseudogene_ids:
                        continue
                    else:
                        if args.basename:
                            if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                                d['ID'] = f"{args.basename}_{d['ID']}"
                            if 'Parent' in d and not d['Parent'].startswith(f"{args.basename}"):
                                d['Parent'] = f"{args.basename}_{d['Parent']}"
                        out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
                else:
                    d = parse_attrs(attr)
                    if args.basename:
                        if 'ID' in d and not d['ID'].startswith(f"{args.basename}"):
                            d['ID'] = f"{args.basename}_{d['ID']}"
                        if 'Parent' in d and not d['Parent'].startswith(f"{args.basename}"):
                            d['Parent'] = f"{args.basename}_{d['Parent']}"
                    out.append("\t".join([seq, src, typ, start, end, score, strand, phase, reconst_attrs(d)]))
    elif in_fmt == "infernal":
        for line in lines:
            if line.startswith("#") or not line.strip():
                out.append(line.rstrip("\n"))
                continue
            parts = line.rstrip("\n").split('\t')
            if len(parts) != 9:
                out.append(line.rstrip("\n"))
                continue
            seq, src, typ, start, end, score, strand, phase, attr = parts
            if typ in csvdata:
                info = csvdata[typ]
                toks = info['Type']
                if any(x.lower() == "gene" for x in toks):
                    child_tok = toks[1] if len(toks) > 1 else typ
                    cnts[(seq, child_tok.lower()+"_g")] += 1
                    if args.basename:
                        gid = f"{args.basename}_{seq}_{normalize_type(child_tok + '_g')}{cnts[(seq, child_tok.lower() + '_g')]}"
                    else:
                        gid = gen_uid(seq, child_tok+"_g", cnts[(seq, child_tok.lower()+"_g")])
                    pattrs = {'ID': gid, 'gene_biotype': child_tok, 'description': info['Description']}
                    out.append("\t".join([seq, src, "gene", start, end, score, strand, phase, reconst_attrs(pattrs)]))
                    cnts[(seq, child_tok.lower())] += 1
                    if args.basename:
                        cid = f"{args.basename}_{seq}_{normalize_type(child_tok)}{cnts[(seq, child_tok.lower())]}"
                    else:
                        cid = gen_uid(seq, child_tok, cnts[(seq, child_tok.lower())])
                    cattrs = {'ID': cid, 'Parent': gid}
                    out.append("\t".join([seq, src, child_tok, start, end, score, strand, phase, reconst_attrs(cattrs)]))
                    exon_id = f"{cid}.exon1"
                    exon_attrs = {'ID': exon_id, 'Parent': cid}
                    out.append("\t".join([seq, src, "exon", start, end, score, strand, phase, reconst_attrs(exon_attrs)]))
                else:
                    cnts[(seq, "bio_reg")] += 1
                    if args.basename:
                        bid = f"{args.basename}_{seq}_{normalize_type('biological_region')}{cnts[(seq, 'bio_reg')]}"
                    else:
                        bid = f"{seq}_{normalize_type('biological_region')}{cnts[(seq, 'bio_reg')]}"
                    battrs = {'ID': bid, 'description': info['Description']}
                    out.append("\t".join([seq, src, "biological_region", start, end, score, strand, phase, reconst_attrs(battrs)]))
            else:
                out.append("\t".join(parts[:8] + [attr if attr else "."]))
    else:
        for line in lines:
            if line.startswith("#") or not line.strip():
                out.append(line.rstrip("\n"))
            else:
                parts = line.rstrip("\n").split('\t')
                if len(parts) != 9:
                    out.append(line.rstrip("\n"))
                else:
                    seq, src, typ, start, end, score, strand, phase, attr = parts
                    out.append("\t".join(parts[:8] + [attr if attr else "."]))
    final = []
    for l in out:
        if l.startswith("##gff-version"):
            continue
        final.append(l)
    final.insert(0, "##gff-version 3")
    return final

def parse_genbank_strands(gbf_path):
    """Parse a GenBank file and return a dictionary of gene names to strand info."""
    gene_strands = {}
    
    with open(gbf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if '/gene=' in line:
                # Extract gene name
                gene_match = re.search(r'/gene="([^"]+)"', line)
                if gene_match:
                    gene_name = gene_match.group(1)
                    # Look back in previous line for strand information
                    prev_line = prev_lines[-1] if 'prev_lines' in locals() else ""
                    strand = '-' if 'complement(' in prev_line else '+'
                    gene_strands[gene_name] = strand
            
            # Keep track of previous lines for context
            if 'prev_lines' not in locals():
                prev_lines = [line]
            else:
                prev_lines.append(line)
                if len(prev_lines) > 2:
                    prev_lines.pop(0)
    
    return gene_strands

def fix_gb_origin_lines(lines, args):
    """
    Rewrite gb_origin GFF so that:
      • All original header lines are preserved (plus '##gff-version 3').
      • region/source features are copied with optional seq renaming.
      • Records form gene hierarchies with IDs {basename}_{seq}_{childfeature}_{geneName}.
      • Transcript lines use child feature type (tRNA or mRNA for CDS).
      • Only one CDS per protein-coding gene, and every child (tRNA/rRNA) also has an exon.
    """
    # Load strand information from GenBank file if provided
    gene_strands = {}
    if args.gbf:
        gene_strands = parse_genbank_strands(args.gbf)
        logging.info(f"Loaded strand information for {len(gene_strands)} genes from GenBank file")
    
    # Collect headers
    orig_headers = [l.rstrip('') for l in lines if l.startswith("##")]
    out_hdr = ["##gff-version 3"] + [h for h in orig_headers if not h.startswith("##gff-version")]
    # Data lines
    data_lines = [l.rstrip('') for l in lines if not l.startswith("##")]
    
    # Define source value to use
    source_value = args.gb_source if args.gb_source else "."

    # Passthrough region/source with seq override
    passthrough = []
    for l in data_lines:
        parts = l.split('\t')
        if len(parts) >= 3 and parts[2].lower() in ("region", "source"):
            parts[0] = args.gb_seq_name or parts[0]
            parts[1] = source_value  # Use the specified source value
            passthrough.append('\t'.join(parts))

    # Group genes and children by start/end coordinates
    gene_groups = {}
    for l in data_lines:
        parts = l.split('\t')
        if len(parts) < 9:
            continue
        seq = args.gb_seq_name or parts[0]
        source = source_value  # Use the specified source value
        feature = parts[2]
        start, end, strand = parts[3], parts[4], parts[6]
        attrs = parse_attrs(parts[8])
        key = (start, end)
        if feature.lower() == 'gene':
            gene_groups[key] = {'gene': {'seq': seq, 'source': source, 'start': start, 'end': end, 'strand': strand, 'attrs': attrs}}
        elif key in gene_groups:
            gene_groups[key]['child'] = {'seq': seq, 'source': source, 'feature': feature, 'start': start, 'end': end, 'strand': strand, 'attrs': attrs}

    output = out_hdr + passthrough
    for key, grp in gene_groups.items():
        grec = grp['gene']
        seq, source, attrs = grec['seq'], grec['source'], grec['attrs']
        gene_name = attrs.get('Name') or attrs.get('gene') or ''
        
        # Use strand from GenBank file if available
        if gene_name in gene_strands:
            grec['strand'] = gene_strands[gene_name]
            logging.info(f"Setting strand for gene {gene_name} to {grec['strand']} from GenBank file")
        
        child = grp.get('child')
        child_feat = child['feature'].lower() if child else 'gene'
        # Build gene ID
        gene_id = f"{args.basename}_{seq}_{gene_name}" if args.basename else f"{seq}_{gene_name}"
        # Determine biotype
        if 'gene_biotype' in attrs:
            biotype = attrs['gene_biotype']
        else:
            biotype = 'protein_coding' if child_feat == 'cds' else child_feat
        # Gene line
        g_attrs = {'ID': gene_id, 'Name': gene_name, 'gene_id': gene_id, 'gene_biotype': biotype}
        output.append('\t'.join([seq, source, 'gene', grec['start'], grec['end'], '.', grec['strand'], '.', reconst_attrs(g_attrs)]))

        # Child line
        # In most of cases, the mitochondria protein-coding gene only has CDS entries, so the child only count CDS entries here.
        # As for ncRNA gene, child is transcript level
        if child:
            # Child features inherit the strand from parent gene
            child['strand'] = grec['strand']
            
            child_count = grp.get('child_count', 0) + 1
            grp['child_count'] = child_count
            if child_feat == 'cds':
                tid = f"{args.basename}_{seq}_{gene_name}.t{child_count}" if args.basename else f"{seq}_{gene_name}.t{child_count}"
            else:
                tid = f"{args.basename}_{seq}_{gene_name}.{child_feat}{child_count}" if args.basename else f"{seq}_{gene_name}.{child_feat}{child_count}"
            t_feat = child['feature'] if child_feat != 'cds' else 'mRNA'
            output.append('\t'.join([child['seq'], child['source'], t_feat, child['start'], child['end'], '.', child['strand'], '.', reconst_attrs({'ID': tid, 'Parent': gene_id})]))
            # CDS + exon for coding
            if biotype == 'protein_coding':
                # CDS
                count = grp.get('cds_count', 0) + 1
                grp['cds_count'] = count
                cds_id = f"{tid}.cds{count}"
                output.append('\t'.join([child['seq'], child['source'], 'CDS', child['start'], child['end'], '.', child['strand'], '0', reconst_attrs({'ID': cds_id, 'Parent': tid})]))
            # Exon for all child types
            ex_count = grp.get('exon_count', 0) + 1
            grp['exon_count'] = ex_count
            exon_id = f"{tid}.exon{ex_count}"
            output.append('\t'.join([child['seq'], child['source'], 'exon', child['start'], child['end'], '.', child['strand'], '.', reconst_attrs({'ID': exon_id, 'Parent': tid})]))
        else:
            # synthetic transcript/exon when no child
            child_count = grp.get('child_count', 0) + 1
            grp['child_count'] = child_count
            tid = f"{args.basename}_{seq}_{gene_name}.t{child_count}" if args.basename else f"{seq}_transcript_{gene_name}.t{child_count}"
            t_feat = 'mRNA' if biotype == 'protein_coding' else 'gene'
            output.append('\t'.join([seq, source, t_feat, grec['start'], grec['end'], '.', grec['strand'], '0', reconst_attrs({'ID': tid, 'Parent': gene_id})]))
            if biotype == 'protein_coding':
                count = grp.get('cds_count', 0) + 1
                grp['cds_count'] = count
                output.append('\t'.join([seq, source, 'CDS', grec['start'], grec['end'], '.', grec['strand'], '0', reconst_attrs({'ID': f"{tid}.cds{count}", 'Parent': tid})]))
            ex_count = grp.get('exon_count', 0) + 1
            grp['exon_count'] = ex_count
            output.append('\t'.join([seq, source, 'exon', grec['start'], grec['end'], '.', grec['strand'], '.', reconst_attrs({'ID': f"{tid}.exon{ex_count}", 'Parent': tid})]))
    return output
    
def fix_gff_lines_main(lines, csv_fp, in_fmt, args):
    if in_fmt == "mitos":
        return fix_mitos_lines(lines, args)
    elif in_fmt == "braker3":
        return fix_braker3(lines, args), []
    elif in_fmt in ("infernal", "trnascan-se"):
        return fix_gff_lines(lines, csv_fp, in_fmt, args), []
    elif in_fmt == ("gb"):
        return fix_gb_origin_lines(lines, args), []
    else:
        return lines, []

def cleanup_for_fix(lines):
    header_count_before = sum(1 for line in lines if line.startswith("#"))
    out = []
    version_found = False
    for line in lines:
        if line.startswith("#!gff-spec-version"):
            continue
        # NEW: drop any MitoHiFi (or otherwise) source‑version header
        if line.startswith("##source-version"):
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
    header_count_after = sum(1 for line in final if line.startswith("#"))
    dropped_header_count = header_count_before - header_count_after
    return final, dropped_header_count, header_count_before

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
        # -- fix malformed coordinates (start > end) --
        start_val, end_val = parts[3], parts[4]
        if start_val.isdigit() and end_val.isdigit():
            if int(start_val) > int(end_val):
                parts[3], parts[4] = end_val, start_val  # swap
        rows.append(parts[:9])

    if not rows:
        return ["##gff-version 3"] + headers

    df = pd.DataFrame(rows, columns=['scaffold','source','feature','start','end','score','strand','phase','attributes'])
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df = df.dropna(subset=['start'])
    df['hierarchy'] = df.apply(lambda row: get_hierarchy_level(row['feature'], row['attributes']), axis=1)
    df['scaffold_key'] = df['scaffold'].apply(
        lambda x: tuple(int(text) if text.isdigit() else text.lower() 
                        for text in re.split(r'(\d+)', x))
    )
    df = df.sort_values(by=['scaffold_key','start','hierarchy'])
    df = df.drop(columns=['scaffold_key','hierarchy'])
    merged = ["\t".join(map(str, row)) for row in df.values]
    unique_hdr = []
    for h in headers:
        if h not in unique_hdr:
            unique_hdr.append(h)
    final = ["##gff-version 3"] + unique_hdr + merged
    return final

def process_file(fmt, fname, args):
    with open(fname, 'r') as f:
        orig_lines = f.readlines()
    input_count = len(orig_lines)
    orig_header_count = sum(1 for line in orig_lines if line.startswith("#"))
    dropped = []
    fmt_lower = fmt.lower()
    conv_count = None
    if fmt_lower == "infernal":
        converted, drops2 = convert_tblout_to_gff(orig_lines, args)
        conv_count = len(converted)
        lines = converted
        dropped = drops2
    elif fmt_lower in ("braker3", "mitos", "trnascan-se", "gb"):
        lines = orig_lines[:]
    fixed, local_dropped = fix_gff_lines_main(lines, args.csv, fmt_lower, args)
    dropped.extend(local_dropped)
    cleaned, _, _ = cleanup_for_fix(fixed)
    header_count_before = orig_header_count
    header_dropped = max(header_count_before - 1, 0)
    base = os.path.basename(fname)
    root, _ = os.path.splitext(base)
    out_fname = os.path.join(args.outdir, f"{root}_{fmt}_fix.gff")
    with open(out_fname, 'w') as outf:
        outf.write("\n".join(cleaned) + "\n")
    print(f"----\nWrote intermediate all-fix file: {out_fname} with {len(cleaned)} lines.")
    return cleaned, dropped, input_count, conv_count, header_count_before

def process_all(args):
    non_mitos_lines = []
    mitos_lines = []
    all_drops = {}
    metrics = {}
    global_log = args.save_filtered_hits if args.save_filtered_hits else os.path.join(args.outdir, "combined_dropout.log")
    with open(global_log, 'w') as gl:
        gl.write("Combined dropout log\n\n")
    for (fmt, fname) in args.I:
        fixed, dropped, input_count, conv_count, orig_header_count = process_file(fmt, fname, args)
        cleaned, _, _ = cleanup_for_fix(fixed)
        header_count_before = orig_header_count
        header_dropped = max(header_count_before - 1, 0)
        base = os.path.basename(fname)
        root, _ = os.path.splitext(base)
        intermediate_fname = os.path.join(args.outdir, f"{root}_{fmt}_fix.gff")
        with open(intermediate_fname, 'r') as inf:
            intermediate_lines = [line.strip() for line in inf.readlines()]
        output_count = len(intermediate_lines)
        if fmt.lower() == "mitos":
            mitos_lines.extend(cleaned)
        else:
            non_mitos_lines.extend(cleaned)
        all_drops[fname] = dropped
        metrics[(fmt, fname)] = (input_count, conv_count, header_dropped, header_count_before)
        if dropped:
            with open(global_log, 'a') as logF:
                logF.write(f"==== Dropouts for {fname} (FORMAT: {fmt}) ====\n")
                for entry in dropped:
                    logF.write(f"{entry[0]}\t{entry[1]}\n")
                logF.write("==== End Dropouts for this file ====\n\n")
    merged_non_mitos = sort_merge_gff_lines(non_mitos_lines)
    mitos_clean = [line for line in mitos_lines if not line.startswith("##gff-version")]
    merged = merged_non_mitos + mitos_clean
    print()
    summary_by_type = defaultdict(list)
    for (fmt, fn) in args.I:
        input_count, conv_count, header_dropped, header_count_before = metrics[(fmt, fn)]
        base = os.path.basename(fn)
        root, _ = os.path.splitext(base)
        intermediate_fname = os.path.join(args.outdir, f"{root}_{fmt}_fix.gff")
        with open(intermediate_fname, 'r') as inf:
            intermediate_lines = [line.strip() for line in inf.readlines()]
        output_count = len(intermediate_lines)
        if fmt.lower() == "mitos":
            mrna_count = sum(1 for line in intermediate_lines if not line.startswith("#") and line.split("\t")[2].lower() == "mrna")
            cds_count = sum(1 for line in intermediate_lines if not line.startswith("#") and line.split("\t")[2].lower() == "cds")
            dropout_count = len(all_drops[fn])
            msg = (f"Mitos Processing\n\n"
                   f"Input File: {fn}\n"
                   f"Total Lines: {input_count}\n\n"
                   f"Output File: {os.path.basename(intermediate_fname)}\n"
                   f"Total Lines: {output_count}\n"
                   f"Dropped header lines: {header_dropped} (from {header_count_before} originally)\n"
                   f"Inserted {mrna_count} mRNA lines and {cds_count} CDS lines\n"
                   f"Dropped {dropout_count} entries (see log: {global_log})\n"
                   f"------------------------------------------------------------")
            summary_by_type["mitos"].append(msg)
        elif fmt.lower() in ("infernal", "trnascan-se"):
            dropped_count = len(all_drops[fn])
            if fmt.lower() == "infernal":
                parent_count = sum(1 for line in intermediate_lines if not line.startswith("#") and line.split("\t")[2].lower() == "gene")
                exon_count = sum(1 for line in intermediate_lines if not line.startswith("#") and line.split("\t")[2].lower() == "exon")
                msg = (f"Infernal Processing\n\n"
                       f"Input File: {fn}\n"
                       f"Total Lines: {input_count}\n\n"
                       f"Processing Steps:\n"
                       f"  - Converted to GFF: {conv_count} lines\n"
                       f"  - Dropped Entries: {dropped_count}\n"
                       f"  - Loaded Rfam CSV: {len(load_csv(args.csv))} rows\n\n"
                       f"Output File: {os.path.basename(intermediate_fname)}\n"
                       f"Total Lines: {output_count}\n"
                       f"Dropped header lines: {header_dropped} (from {header_count_before} originally)\n"
                       f"Drop Log: {global_log}\n\n"
                       f"Added {parent_count} parent lines, {exon_count} exon lines\n"
                       f"------------------------------------------------------------")
                summary_by_type["infernal"].append(msg)
            elif fmt.lower() == "trnascan-se":
                gene_count = sum(1 for line in intermediate_lines if not line.startswith("#") and line.split("\t")[2].lower() == "gene")
                msg = (f"tRNAscan-SE Processing\n\n"
                       f"Input File: {fn}\n"
                       f"Total Lines: {input_count}\n\n"
                       f"Output File: {os.path.basename(intermediate_fname)}\n"
                       f"Total Lines: {output_count}\n"
                       f"Dropped header lines: {header_dropped} (from {header_count_before} originally)\n"
                       f"Inserted {gene_count} gene entries\n"
                       f"------------------------------------------------------------")
                summary_by_type["trnascan-se"].append(msg)
        else:
            msg = (f"{fmt.upper()} Processing\n\n"
                   f"Input File: {fn}\n"
                   f"Total Lines: {input_count}\n\n"
                   f"Output File: {os.path.basename(intermediate_fname)}\n"
                   f"Total Lines: {output_count}\n"
                   f"------------------------------------------------------------")
            summary_by_type[fmt.lower()].append(msg)
    for file_type, summaries in summary_by_type.items():
        print(f"\n==== Summary for {file_type.upper()} Files ====\n")
        for msg in summaries:
            print(msg)
    return merged, all_drops, metrics

def validate_args(args):
    """Validate command-line arguments and ensure required combinations are present."""
    gb_formats = False
    for fmt, _ in args.I:
        if fmt.lower() == "gb":
            gb_formats = True
            break
    
    if gb_formats and not args.gbf:
        raise ValueError("The --gbf option is required when using the 'gb' format with -I")
    
    return args

def main():
    parser = argparse.ArgumentParser(
        description="GFF combine & fix tool. Processes all input files, writes intermediate _fix.gff files, then merges them into merged_fix.gff."
    )
    parser.add_argument("-I", nargs=2, metavar=("FORMAT","FILE"), action="append", required=True,
                        help="(FORMAT in {infernal, trnascan-se, braker3, mitos, gb})")
    parser.add_argument("--csv", type=str, default="Rfam_15_0.csv",
                        help="CSV for lookup (default=Rfam_15_0.csv).")
    parser.add_argument("--outdir", type=str, required=True,
                        help="Output directory. Also writes merged_fix.gff.")
    parser.add_argument("--save-filtered-hits", type=str, default=None,
                        help="Prefix (or filename) for dropped-hits log (global for all files).")
    parser.add_argument("--min-bit", "-T", type=float, default=None,
                        help="Min bit score for infernal.")
    parser.add_argument("--max-evalue", "-E", type=float, default=None,
                        help="Max evalue for infernal.")
    parser.add_argument("--fmt2", action="store_true", help="Infernal tblout with --fmt2 (>=27 fields).")
    parser.add_argument("--ignore", type=str, help="Skip specified hits in infernal. Separate multiple entries with ',' (e.g., tRNA,rRNA)")
    parser.add_argument("--basename", type=str, default=None,
                        help="Prepend to generated gene IDs, e.g. 'M_tai' => 'M_tai_rrna.g1'.")
    parser.add_argument("--source", type=str, default=None,
                        help="Specify source field for infernal-to-GFF conversion. If provided, this value is used in the second column; otherwise, default is 'cmscan'.")
    parser.add_argument("--gb-seq-name", type=str, default=None,
                    help="Override sequence name for gb_origin input GFF.")
    parser.add_argument("--gb-source", type=str, default=None,
                    help="Specify source field for GB input. If provided, this value is used in the second column; otherwise, default is '.'.")
    parser.add_argument("--gbf", type=str, default=None,
                    help="Path to GenBank file (.gb/.gbf) to extract strand information. Required when using '-I gb'.")
    
    args = parser.parse_args()

    # Validate the arguments
    try:
        args = validate_args(args)
    except ValueError as e:
        parser.error(str(e))  # This will print the error message and exit
    
    csv_data = load_csv(args.csv)
    merged, all_drops, metrics = process_all(args)
    merged_fname = os.path.join(args.outdir, "merged_fix.gff")
    with open(merged_fname, 'w') as outF:
        outF.write("\n".join(merged) + "\n")
    merged_count = len(merged)
    print(f"\nFinal Merging\n\nMerged File Created: {os.path.basename(merged_fname)}\nTotal Lines in Merged File: {merged_count}\n------------------------------------------------------------")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
