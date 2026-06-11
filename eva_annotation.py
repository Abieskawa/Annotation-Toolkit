#!/usr/bin/env python3
import os
import re
import sys
import subprocess
import argparse
import statistics
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# -----------------------------------------------------------------------
# Utility: Parse GFF attributes
# -----------------------------------------------------------------------
def parse_attributes(attr_str):
    """Parse GFF attributes in the form key=value;... and return a dict."""
    attrs = {}
    attr_str = attr_str.strip().rstrip(';')
    if not attr_str:
        return attrs
    for kv in attr_str.split(';'):
        kv = kv.strip()
        if '=' in kv:
            key, val = kv.split('=', 1)
            attrs[key] = val
    return attrs



# Backward-compatible alias (some helper code expects this name)
def _parse_gff_attributes(attr_str):
    return parse_attributes(attr_str)

# -----------------------------------------------------------------------
# Utility: Build protein_id -> gene_id mapping from an NCBI/RefSeq GFF3
# -----------------------------------------------------------------------
def _normalize_gff_id(raw_id):
    """Strip common NCBI prefixes (gene-/rna-) from a GFF3 ID string."""
    if raw_id is None:
        return None
    raw_id = raw_id.strip()
    for prefix in ("gene-", "rna-", "mrna-", "mRNA-", "transcript-"):
        if raw_id.startswith(prefix):
            return raw_id[len(prefix):]
    return raw_id

def _split_attr_list(val):
    """Split a comma-separated GFF3 attribute value into a clean list."""
    if not val:
        return []
    return [x.strip() for x in val.split(",") if x.strip()]


def build_protein_to_gene_map_from_gff(gff_file, unresolved_out=None):
    """
    Build a mapping of protein IDs to gene labels using GFF3 parentage.

    Robust strategy:
      1) Index all feature IDs and their Parent relationships (ID -> Parent(s)).
      2) Cache gene feature labels (for type=gene).
      3) For each CDS that has a protein_id, walk up the Parent chain (BFS) until a gene is found.
         This works even when intermediate feature types are non-standard (e.g., V_gene_segment).

    Notes about NCBI/RefSeq GFF3:
      - Many CDS lines have protein_id=<XP_...> (or similar).
      - Some biological features (e.g., immunoglobulin segments) may have CDS entries WITHOUT protein_id.
        Those cannot be mapped to protein FASTA/diamond qseqid; we list them separately when
        `unresolved_out` is provided.

    Gene label priority (from the gene feature attributes):
        gene  >  locus_tag  >  Name  >  ID (normalized by stripping 'gene-')

    Args:
        gff_file (str): Path to a GFF3 file.
        unresolved_out (str|None): If provided, write a TSV listing CDS entries that:
            - have protein_id but cannot be resolved to a gene, and/or
            - lack protein_id entirely (common for special cases like rearrangement-required segments).

    Returns:
        dict: {protein_id: gene_label}
    """
    id_to_parents = {}
    id_to_type = {}
    gene_id_to_label = {}
    gene_label_to_id = {}
    # Some non-gene features carry useful labels; keep as fallback.
    id_to_geneish_label = {}

    # protein_id -> list of CDS entry dicts
    cds_by_protein = {}
    # CDS entries that do not have protein_id
    cds_missing_protein = []

    def _norm_gene_id(feat_id):
        if not feat_id:
            return feat_id
        return feat_id[5:] if feat_id.startswith("gene-") else feat_id

    def _gene_label_from_attrs(attrs, feat_id=None):
        return attrs.get("gene") or attrs.get("locus_tag") or attrs.get("Name") or _norm_gene_id(feat_id)

    with open(gff_file, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attr_str = parts
            attrs = _parse_gff_attributes(attr_str)

            feat_id = attrs.get("ID")
            if feat_id:
                id_to_type[feat_id] = feature_type

                parents = _split_attr_list(attrs.get("Parent"))
                if parents:
                    id_to_parents[feat_id] = parents

                # Cache any gene-ish label as fallback for intermediate nodes
                lbl = attrs.get("gene") or attrs.get("locus_tag") or attrs.get("Name")
                if lbl:
                    id_to_geneish_label[feat_id] = lbl

            # Cache gene labels (authoritative) and reverse lookup from common labels
            if feature_type == "gene" and feat_id:
                glabel = _gene_label_from_attrs(attrs, feat_id)
                gene_id_to_label[feat_id] = glabel
                for key in (attrs.get("gene"), attrs.get("locus_tag"), attrs.get("Name"), glabel, _norm_gene_id(feat_id), feat_id):
                    if key and key not in gene_label_to_id:
                        gene_label_to_id[key] = feat_id
                continue

            # CDS entries: store protein_id and parent(s)
            if feature_type == "CDS":
                protein_id = attrs.get("protein_id")
                parents = _split_attr_list(attrs.get("Parent"))
                cds_id = attrs.get("ID")

                entry = {
                    "protein_id": protein_id or "",
                    "cds_id": cds_id or "",
                    "parent_id": parents[0] if parents else "",
                    "parent_ids": parents or [],
                    "seqid": seqid,
                    "start": int(start) if start.isdigit() else start,
                    "end": int(end) if end.isdigit() else end,
                    "strand": strand,
                    "attrs_raw": attr_str,
                    "gene_attr": attrs.get("gene") or attrs.get("locus_tag") or attrs.get("Name"),
                }

                if not protein_id:
                    cds_missing_protein.append(entry)
                else:
                    cds_by_protein.setdefault(protein_id, []).append(entry)

    def _find_gene_label_from_feature(start_id):
        """
        Walk parent relationships upwards until a gene feature is found.
        Returns (gene_label, found_gene_feat_id) or (None, None).
        """
        if not start_id:
            return None, None

        q = [start_id]
        visited = set()

        while q:
            cur = q.pop(0)
            if not cur or cur in visited:
                continue
            visited.add(cur)

            # Direct gene feature
            if cur in gene_id_to_label:
                return gene_id_to_label[cur], cur

            # Sometimes gene IDs exist but weren't cached as type=gene for odd files
            if cur.startswith("gene-"):
                # Try best-effort
                return id_to_geneish_label.get(cur) or _norm_gene_id(cur), cur

            # If this node carries a gene-ish label, keep it as potential fallback but
            # continue walking to prefer an actual gene feature.
            # (We don't return immediately because intermediate nodes may have gene=, but
            # the gene feature label might be better.)
            q.extend(id_to_parents.get(cur, []))

        return None, None

    protein_to_gene = {}
    unresolved_rows = []
    unresolved_count = 0

    for protein_id, entries in cds_by_protein.items():
        gene_label = None
        gene_feat = None
        reason = None

        # Try parents from all CDS entries (some proteins have multiple CDS parts)
        parent_ids = []
        for e in entries:
            parent_ids.extend(e.get("parent_ids", []))
        # Preserve order but unique
        seen = set()
        parent_ids = [p for p in parent_ids if p and not (p in seen or seen.add(p))]

        for pid in parent_ids:
            gene_label, gene_feat = _find_gene_label_from_feature(pid)
            if gene_label:
                break

        if not gene_label:
            # Fallback: use gene-like attribute carried on CDS itself (common in NCBI GFF)
            for e in entries:
                if e.get("gene_attr"):
                    gene_label = e["gene_attr"]
                    gene_feat = gene_label_to_id.get(gene_label)
                    if not gene_feat and gene_label and not str(gene_label).startswith("gene-"):
                        gene_feat = f"gene-{gene_label}"
                    reason = "Used CDS gene/locus_tag/Name attribute as fallback (no resolvable gene parentage)"
                    break

        if gene_label:
            protein_to_gene[protein_id] = gene_feat if gene_feat else gene_label
        else:
            unresolved_count += 1
            # Pick one representative CDS entry for reporting
            e = entries[0]
            parent0 = e.get("parent_id", "")
            if not parent0:
                reason = "CDS has no Parent attribute"
            elif parent0 in id_to_type:
                reason = f"Parent feature '{parent0}' exists (type='{id_to_type[parent0]}') but no path to a gene feature was found"
            else:
                reason = f"Parent feature '{parent0}' not found among GFF IDs"

            unresolved_rows.append(
                {
                    "protein_id": protein_id,
                    "cds_id": e.get("cds_id", ""),
                    "parent_id": parent0,
                    "seqid": e.get("seqid", ""),
                    "start": e.get("start", ""),
                    "end": e.get("end", ""),
                    "strand": e.get("strand", ""),
                    "reason": reason,
                    "cds_attributes": e.get("attrs_raw", ""),
                }
            )

    # Also report CDS lacking protein_id so users can inspect special cases
    missing_count = len(cds_missing_protein)

    # Create proxy mappings for CDS lacking protein_id (e.g., immune segments / rearrangement-required).
    # Use gene/locus_tag/Name carried on CDS (if present) as a "protein ID" key and map it to the gene feature ID.
    proxy_added = 0
    if missing_count > 0:
        for e in cds_missing_protein:
            proxy_key = e.get("gene_attr") or ""
            if not proxy_key:
                continue

            # Try to resolve gene via Parent chain starting from the CDS's parent feature.
            gene_label, gene_feat = _find_gene_label_from_feature(e.get("parent_id"))
            if not gene_feat:
                # Reverse-lookup from gene label -> gene feature ID
                gene_feat = gene_label_to_id.get(proxy_key) or gene_label_to_id.get(gene_label)
            if not gene_feat and proxy_key and not str(proxy_key).startswith("gene-"):
                gene_feat = f"gene-{proxy_key}"

            if gene_feat and proxy_key not in protein_to_gene:
                protein_to_gene[proxy_key] = gene_feat
                proxy_added += 1
    if unresolved_out and (unresolved_count > 0 or missing_count > 0):
        try:
            with open(unresolved_out, "w", encoding="utf-8") as out:
                out.write("\t".join(
                    ["protein_id", "cds_id", "parent_id", "seqid", "start", "end", "strand", "reason", "cds_attributes"]
                ) + "\n")

                for r in unresolved_rows:
                    out.write("\t".join([
                        r["protein_id"],
                        r["cds_id"],
                        r["parent_id"],
                        r["seqid"],
                        str(r["start"]),
                        str(r["end"]),
                        r["strand"],
                        r["reason"],
                        r["cds_attributes"],
                    ]) + "\n")

                for e in cds_missing_protein:
                    out.write("\t".join([
                        "",  # protein_id missing
                        e.get("cds_id", ""),
                        e.get("parent_id", ""),
                        e.get("seqid", ""),
                        str(e.get("start", "")),
                        str(e.get("end", "")),
                        e.get("strand", ""),
                        "CDS has no protein_id attribute (common for special/rearrangement-required annotations)",
                        e.get("attrs_raw", ""),
                    ]) + "\n")

            print(f"[INFO] Wrote GFF mapping issues to: {unresolved_out}", file=sys.stderr)
        except Exception as e:
            print(f"[WARN] Failed to write unresolved list to '{unresolved_out}': {e}", file=sys.stderr)

    if missing_count > 0:
        print(
            f"[WARN] {missing_count} CDS feature(s) lack protein_id in: {gff_file} "
            f"(true protein_id is missing; created {proxy_added} proxy mapping(s) using gene/locus_tag/Name as protein IDs).",
            file=sys.stderr,
        )

    if unresolved_count > 0:
        print(
            f"[WARN] {unresolved_count} protein_id(s) from CDS could not be resolved back to a gene via Parent chains in: {gff_file}",
            file=sys.stderr,
        )

    return protein_to_gene

def map_protein_to_gene(protein_id, gff_map=None, fallback=None):
    """Resolve a gene ID from a protein ID using an optional GFF-derived mapping."""
    if gff_map and protein_id in gff_map:
        return gff_map[protein_id]
    return fallback(protein_id) if fallback else protein_id

# -----------------------------------------------------------------------
# Helper: Ensure Unique Entries Based on 'ID'
# -----------------------------------------------------------------------
def ensure_unique_entries(df, entry_type):
    """
    Checks the DataFrame for duplicate IDs, reports them,
    and returns the original DataFrame without aggregation.
    """
    dup = df[df['ID'].duplicated(keep=False)]
    warning_msg = ""
    if not dup.empty:
        warning_msg = f"Warning: Duplicate {entry_type} IDs found: {dup['ID'].unique().tolist()}\n"
    return df, warning_msg

# -----------------------------------------------------------------------
# Helper: Custom binning with non-overlapping ranges and custom labels
# -----------------------------------------------------------------------
def custom_binning(lengths, bin_size, thresh):
    """
    Bins the data with the following logic:
      - First bin: values < bin_size; label: "<{bin_size}"
      - Subsequent bins: each covers values from (previous bin's maximum + 1) to the current boundary.
      - An extra bin: values > thresh; label: ">{thresh}"
    Returns three lists: centers (for bar placement), counts, and labels.
    """
    counts = []
    labels = []
    centers = []
    first_count = sum(1 for x in lengths if x < bin_size)
    counts.append(first_count)
    labels.append(f"<{bin_size}")
    first_center = (0 + (bin_size - 1)) / 2
    centers.append(first_center)
    
    lower = bin_size + 1
    upper = bin_size * 2
    while upper <= thresh:
        cnt = sum(1 for x in lengths if lower <= x <= upper)
        counts.append(cnt)
        labels.append(f"{lower}-{upper}")
        center = (lower + upper) / 2
        centers.append(center)
        lower = upper + 1
        upper += bin_size
        
    extra_count = sum(1 for x in lengths if x > thresh)
    counts.append(extra_count)
    labels.append(f">{thresh}")
    center_extra = thresh + bin_size / 2
    centers.append(center_extra)
    
    return centers, counts, labels

# -----------------------------------------------------------------------
# Utility: Plot a length histogram with custom bin labels
# -----------------------------------------------------------------------
def plot_bar(lengths, med, out_file, bin_size, thresh, xlabel):
    """
    Plots a histogram of 'lengths' using custom binning and saves the plot.
    """
    if not lengths:
        print("No lengths to plot.")
        return

    centers, counts, labels = custom_binning(lengths, bin_size, thresh)
    plt.figure(figsize=(8,6))
    bar_width = bin_size * 0.8
    plt.bar(centers, counts, width=bar_width, align='center', alpha=0.7, edgecolor='black')
    plt.axvline(med, color='red', linestyle='dashed', linewidth=2, label=f"Median: {int(med)}")
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.title("Length Distribution")
    plt.legend()
    plt.xticks(centers, labels)
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()

# -----------------------------------------------------------------------
# Utility: Check start/stop codons for an mRNA
# -----------------------------------------------------------------------
def check_codon(mrna, grp, codon_type):
    """
    Checks that an mRNA has exactly one codon of 'codon_type'.
    """
    mrna_id     = mrna["ID"]
    mrna_seqid  = mrna["seqid"]
    mrna_start  = mrna["start"]
    mrna_end    = mrna["end"]
    mrna_strand = mrna["strand"]

    if grp is None or mrna_id not in grp.groups:
        return f"mRNA {mrna_id} has no {codon_type} codon defined."

    subset = grp.get_group(mrna_id)
    count = subset[
        (subset["start"] >= mrna_start) &
        (subset["end"]   <= mrna_end)   &
        (subset["strand"] == mrna_strand) &
        (subset["seqid"]  == mrna_seqid)
    ].shape[0]

    if count == 0:
        return f"mRNA {mrna_id} has no {codon_type} codon defined."
    if count > 1:
        return f"mRNA {mrna_id} has multiple {codon_type} codons ({count})."
    return ""

# -----------------------------------------------------------------------
# New Helper: Merge overlapping intervals
# -----------------------------------------------------------------------
def merge_intervals(intervals):
    """Merge overlapping intervals from a list of (start, end) tuples."""
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev = merged[-1]
        if current[0] <= prev[1] + 1:
            merged[-1] = (prev[0], max(prev[1], current[1]))
        else:
            merged.append(current)
    return merged

# -----------------------------------------------------------------------
# New Helper: Calculate total base-pairs from a DataFrame of intervals
# -----------------------------------------------------------------------
def calculate_total_interval_bp(df_intervals):
    """
    Calculate total bp covered by intervals (merging overlapping ones) in a dataframe.
    The dataframe should have columns: 'seqid', 'start', and 'end'.
    """
    total_bp = 0
    for seqid, group in df_intervals.groupby("seqid"):
        intervals = list(zip(group['start'], group['end']))
        merged = merge_intervals(intervals)
        for start, end in merged:
            total_bp += (end - start + 1)
    return total_bp

# -----------------------------------------------------------------------
# New Helper: Calculate genome size from a FASTA file, with optional filtering.
# -----------------------------------------------------------------------
def calculate_genome_size(fasta_file, include_set=None):
    """
    Calculate total genome bp by summing the lengths of sequence lines in a FASTA file.
    
    If include_set is provided (a set of scaffold IDs), only sequences whose header 
    (the first token after '>') is in include_set will be included.
    """
    total_size = 0
    current_seq_id = None
    seq_lines = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_seq_id is not None:
                    if include_set is None or current_seq_id in include_set:
                        total_size += len("".join(seq_lines))
                current_seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if current_seq_id is not None:
            if include_set is None or current_seq_id in include_set:
                total_size += len("".join(seq_lines))
    return total_size

# -----------------------------------------------------------------------
# Structural Annotation Mode 
# -----------------------------------------------------------------------
def run_structural(args):
    cols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    df = pd.read_csv(args.gff, sep='\t', comment='#', names=cols, header=None,
                     dtype={'start': int, 'end': int}, na_filter=False)

    df['parsed_attrs'] = df['attributes'].apply(parse_attributes)
    df['ID'] = df['parsed_attrs'].apply(lambda d: d.get('ID'))
    df['Parent'] = df['parsed_attrs'].apply(lambda d: d.get('Parent'))

    # Subset feature DataFrames and ensure unique entries.
    genes_df = df[df['type'] == "gene"].copy()
    genes_df, gene_dup_msg = ensure_unique_entries(genes_df, "gene")
    pseudos_df = df[df['type'] == "pseudogene"].copy()
    pseudos_df, pseudo_dup_msg = ensure_unique_entries(pseudos_df, "pseudogene")
    mRNA_df = df[df['type'] == 'mRNA'].copy()
    mRNA_df, mrna_dup_msg = ensure_unique_entries(mRNA_df, "mRNA")
    exon_df = df[df['type'] == "exon"].copy()
    exon_df, exon_dup_msg = ensure_unique_entries(exon_df, "exon")
    
    n_genes = genes_df.shape[0]
    n_pseudos = pseudos_df.shape[0]

    gene_biotype_available = False
    if not genes_df.empty:
        genes_df["gene_biotype"] = genes_df["parsed_attrs"].apply(lambda d: d.get("gene_biotype", "NA"))
        if set(genes_df["gene_biotype"].unique()) != {"NA"}:
            gene_biotype_available = True

    skip_codon_check = False
    if "Gnomon" in df["source"].unique():
        skip_codon_check = True

    start_errs = []
    stop_errs = []
    if not skip_codon_check:
        if gene_biotype_available:
            gene_info = genes_df[['ID', 'gene_biotype']].rename(columns={'ID': 'gene'})
            mRNA_df = mRNA_df.merge(gene_info, left_on='Parent', right_on='gene', how='left')
        start_codons = df[df['type'] == 'start_codon']
        stop_codons = df[df['type'] == 'stop_codon']
        start_grp = start_codons.groupby("Parent") if not start_codons.empty else None
        stop_grp = stop_codons.groupby("Parent") if not stop_codons.empty else None
        if gene_biotype_available:
            pc_mrnas = mRNA_df[mRNA_df["gene_biotype"] == "protein_coding"]
            for _, r in pc_mrnas.iterrows():
                err = check_codon(r, start_grp, "start")
                if err:
                    start_errs.append(err)
                err = check_codon(r, stop_grp, "stop")
                if err:
                    stop_errs.append(err)
        else:
            for _, r in mRNA_df.iterrows():
                err = check_codon(r, start_grp, "start")
                if err:
                    start_errs.append(err)
                err = check_codon(r, stop_grp, "stop")
                if err:
                    stop_errs.append(err)

    if gene_biotype_available:
        gene_info = genes_df[['ID', 'gene_biotype']].rename(columns={'ID': 'gene'})
        mRNA_biotype = df[df['type'] == 'mRNA'].merge(gene_info, left_on='Parent', right_on='gene', how='left')
        pc_mrnas = mRNA_biotype[mRNA_biotype['gene_biotype'] == 'protein_coding'].copy()
        pc_exons = exon_df.merge(pc_mrnas[['ID']], left_on='Parent', right_on='ID', how='inner')
        n_exon = pc_exons.shape[0]
        pc_mrnas_map = pc_mrnas[['ID', 'Parent']].rename(columns={'Parent': 'Gene_ID'})
        pc_exons_merged = pc_exons.merge(pc_mrnas_map, left_on='Parent', right_on='ID', how='left', suffixes=('_exon','_mRNA'))
        pc_exon_counts = pc_exons_merged.groupby('Gene_ID').size()
        single_exon_genes = 0
        pc_genes = genes_df[genes_df["gene_biotype"] == "protein_coding"]
        for gid in pc_genes['ID']:
            if pc_exon_counts.get(gid, 0) == 1:
                single_exon_genes += 1
    else:
        n_exon = exon_df.shape[0]
        mRNA_gene_map = mRNA_df[['ID','Parent']].rename(columns={'Parent':'Gene_ID'})
        exons_merged = exon_df.merge(mRNA_gene_map, left_on='Parent', right_on='ID', how='left', suffixes=('_exon','_mRNA'))
        exon_counts = exons_merged.groupby('Gene_ID').size()
        unique_genes = genes_df.copy()
        single_exon_genes = 0
        for gid in unique_genes['ID']:
            if exon_counts.get(gid, 0) == 1:
                single_exon_genes += 1

    n_mrna = mRNA_df.shape[0]

    # Compute gene lengths (avoid SettingWithCopyWarning).
    all_genes = genes_df.copy()
    if not all_genes.empty:
        all_genes = all_genes.copy()
        all_genes.loc[:, "gene_length"] = (all_genes["end"] - all_genes["start"]).abs()
        gene_lengths = all_genes["gene_length"].tolist()
        median_len = statistics.median(gene_lengths) if gene_lengths else 0
    else:
        gene_lengths = []
        median_len = 0

    # Compute full genome size from FASTA (using all scaffolds in the GFF).
    all_scaffolds = set(df["seqid"].unique())
    full_genome_size = calculate_genome_size(args.fasta, include_set=all_scaffolds)

    # If split scaffolds are provided, compute group-specific genome sizes.
    split_set = set(args.split_scaffolds) if args.split_scaffolds else None

    if split_set:
        included_scaffolds = all_scaffolds - split_set
        genome_size_included = calculate_genome_size(args.fasta, include_set=included_scaffolds)
        genome_size_split = calculate_genome_size(args.fasta, include_set=split_set)

        genes_df_included = genes_df[~genes_df["seqid"].isin(split_set)].copy()
        pseudos_df_included = pseudos_df[~pseudos_df["seqid"].isin(split_set)].copy()
        mRNA_df_included = mRNA_df[~mRNA_df["seqid"].isin(split_set)].copy()
        exon_df_included = exon_df[~exon_df["seqid"].isin(split_set)].copy()
        n_genes_included = genes_df_included.shape[0]
        n_pseudos_included = pseudos_df_included.shape[0]
        n_mrna_included = mRNA_df_included.shape[0]
        n_exon_included = exon_df_included.shape[0]
        if not genes_df_included.empty:
            genes_df_included = genes_df_included.copy()
            genes_df_included.loc[:, "gene_length"] = (genes_df_included["end"] - genes_df_included["start"]).abs()
            lengths_included = genes_df_included["gene_length"].tolist()
            median_included = statistics.median(lengths_included) if lengths_included else 0
        else:
            median_included = 0

        # Calculate protein-coding average length for Non-Split scaffolds.
        if gene_biotype_available:
            pc_included = genes_df_included[genes_df_included["gene_biotype"]=="protein_coding"]
            if not pc_included.empty:
                avg_length_pc_included = pc_included["gene_length"].mean()
            else:
                avg_length_pc_included = 0

        # (Pseudogene average length for Non-Split scaffolds)
        if not pseudos_df_included.empty:
            avg_length_pseudogene_included = pseudos_df_included["end"].sub(pseudos_df_included["start"]).abs().mean()
        else:
            avg_length_pseudogene_included = 0

        genes_df_split = genes_df[genes_df["seqid"].isin(split_set)].copy()
        pseudos_df_split = pseudos_df[pseudos_df["seqid"].isin(split_set)].copy()
        mRNA_df_split = mRNA_df[mRNA_df["seqid"].isin(split_set)].copy()
        exon_df_split = exon_df[exon_df["seqid"].isin(split_set)].copy()
        n_genes_split = genes_df_split.shape[0]
        n_pseudos_split = pseudos_df_split.shape[0]
        n_mrna_split = mRNA_df_split.shape[0]
        n_exon_split = exon_df_split.shape[0]
        if not genes_df_split.empty:
            genes_df_split = genes_df_split.copy()
            genes_df_split.loc[:, "gene_length"] = (genes_df_split["end"] - genes_df_split["start"]).abs()
            lengths_split = genes_df_split["gene_length"].tolist()
            median_split = statistics.median(lengths_split) if lengths_split else 0
        else:
            median_split = 0

        # Calculate protein-coding average length for Split scaffolds.
        if gene_biotype_available:
            pc_split = genes_df_split[genes_df_split["gene_biotype"]=="protein_coding"]
            if not pc_split.empty:
                avg_length_pc_split = pc_split["gene_length"].mean()
            else:
                avg_length_pc_split = 0
        # (Pseudogene average length for Split scaffolds.)
        if not pseudos_df_split.empty:
            avg_length_pseudogene_split = pseudos_df_split["end"].sub(pseudos_df_split["start"]).abs().mean()
        else:
            avg_length_pseudogene_split = 0
    else:
        n_genes_included = n_genes
        n_pseudos_included = n_pseudos
        n_mrna_included = n_mrna
        n_exon_included = n_exon
        median_included = median_len
        genome_size_included = full_genome_size
        genes_df_included = genes_df.copy()
        if gene_biotype_available:
            pc_included = genes_df.copy()  # all genes_df
            pc_included = pc_included[pc_included["gene_biotype"]=="protein_coding"]
            if not pc_included.empty:
                avg_length_pc_included = pc_included["end"].sub(pc_included["start"]).abs().mean()
            else:
                avg_length_pc_included = 0
        if not pseudos_df.empty:
            avg_length_pseudogene_included = pseudos_df["end"].sub(pseudos_df["start"]).abs().mean()
        else:
            avg_length_pseudogene_included = 0

    base = args.basename
    outdir = args.outdir
    qc_output_file = os.path.join(outdir, f"{base}_qc.txt")
    qc_stat_file   = os.path.join(outdir, f"{base}_statistic.txt")
    qc_extra_file  = os.path.join(outdir, f"{base}_additional_validation.txt")

    subprocess.run([
        "gff3_QC",
        "--gff", args.gff,
        "--fasta", args.fasta,
        "--output", qc_output_file,
        "--statistic", qc_stat_file
    ], check=True)

    with open(qc_extra_file, 'w') as f:
        # Overall header
        f.write("==== SUMMARY OVERALL ====\n\n")
        if split_set:
            f.write(">> SUMMARY (Non-Split scaffolds):\n")
            f.write(f"Genes (unique): {n_genes_included}\n")
            f.write(f"Pseudogenes (unique): {n_pseudos_included}\n")
            f.write(f"Single-exon genes (Protein-coding): {single_exon_genes}\n")
            f.write(f"mRNAs (unique): {n_mrna_included}\n")
            f.write(f"Exons (unique, Protein-coding if only available): {n_exon_included}\n")
            f.write(f"Median gene length (all genes): {median_included}\n")
            f.write(f"Genome size (Non-Split): {genome_size_included}\n")
            if gene_biotype_available:
                f.write(f"Protein-coding genes average length: {avg_length_pc_included:.4f} bp\n")
            f.write("\n====\n\n")
            f.write(">> SUMMARY (Split scaffolds):\n")
            f.write(f"Genes (unique): {n_genes_split}\n")
            f.write(f"Pseudogenes (unique): {n_pseudos_split}\n")
            f.write(f"mRNAs (unique): {n_mrna_split}\n")
            f.write(f"Exons (unique): {n_exon_split}\n")
            f.write(f"Median gene length (all genes): {median_split}\n")
            f.write(f"Genome size (Split): {genome_size_split}\n")
            if gene_biotype_available:
                f.write(f"Protein-coding genes average length: {avg_length_pc_split:.4f} bp\n")
            f.write(f"Pseudogenes average length: {avg_length_pseudogene_split:.4f} bp\n")
        else:
            f.write(f"Genes (unique): {n_genes_included}\n")
            f.write(f"Pseudogenes (unique): {n_pseudos_included}\n")
            f.write(f"Single-exon genes (Protein-coding if available): {single_exon_genes}\n")
            f.write(f"mRNAs (unique): {n_mrna_included}\n")
            f.write(f"Exons (unique, Protein-coding if available): {n_exon_included}\n")
            f.write(f"Median gene length (all genes): {median_included}\n")
            f.write(f"Genome size: {genome_size_included}\n")
            if gene_biotype_available:
                f.write(f"Protein-coding genes average length: {avg_length_pc_included:.4f} bp\n")
            f.write(f"Pseudogenes average length: {avg_length_pseudogene_included:.4f} bp\n")
        
        f.write("\n-------------------------\n")
        if gene_dup_msg:
            f.write(gene_dup_msg)
        if pseudo_dup_msg:
            f.write(pseudo_dup_msg)
        if mrna_dup_msg:
            f.write(mrna_dup_msg)
        if exon_dup_msg:
            f.write(exon_dup_msg)
        f.write("\n-------------------------\n")
        if not skip_codon_check:
            f.write(f"\nStart codon issues: {len(start_errs)}\n")
            f.write(f"Stop codon issues: {len(stop_errs)}\n")
            f.write("****** Details ******\n")
            f.write("Start codon issues\n")
            for err in start_errs:
                f.write(err + "\n")
            f.write("-------------------------\n")
            f.write("Stop codon issues\n")
            for err in stop_errs:
                f.write(err + "\n")
            f.write("-------------------------\n")
        else:
            f.write("Skipping codon checks (gene_biotype present or source is 'Gnomon').\n")
        
        f.write("\n==== Gene Biotype & Coverage Summary ====\n")
        if gene_biotype_available:
            # Non-Split gene biotype info.
            genes_df_excl = genes_df_included.copy()
            genes_df_excl.loc[:, "gene_length"] = (genes_df_excl["end"] - genes_df_excl["start"]).abs()
            biotype_counts_excl = genes_df_excl["gene_biotype"].value_counts()
            bp_by_biotype = genes_df_excl.groupby("gene_biotype")["gene_length"].sum()
            f.write("\nGene biotype counts and total bp (Non-Split):\n")
            for bt in biotype_counts_excl.index:
                count = biotype_counts_excl[bt]
                total_bp = bp_by_biotype.get(bt, 0)
                ratio_biotype = (total_bp / genome_size_included * 100) if genome_size_included else 0
                f.write(f"  {bt}: {count} (Total bp: {total_bp}, Ratio: {ratio_biotype:.4f}%)\n")
            
            f.write("\nAverage gene length for non-protein-coding genes (Non-Split):\n")
            non_pc = genes_df_excl[genes_df_excl["gene_biotype"] != "protein_coding"]
            if not non_pc.empty:
                for bt in biotype_counts_excl.index:
                    if bt == "protein_coding":
                        continue
                    subset = non_pc[non_pc["gene_biotype"] == bt]
                    if not subset.empty:
                        avg_len = subset["gene_length"].mean()
                        f.write(f"  {bt}: {avg_len:.4f} bp\n")
            else:
                f.write("  No non-protein-coding genes available.\n")
            
            if split_set:
                genes_df_split_copy = genes_df_split.copy()
                genes_df_split_copy.loc[:, "gene_length"] = (genes_df_split_copy["end"] - genes_df_split_copy["start"]).abs()
                biotype_counts_split = genes_df_split_copy["gene_biotype"].value_counts()
                bp_by_biotype_split = genes_df_split_copy.groupby("gene_biotype")["gene_length"].sum()
                f.write("\nGene biotype counts and total bp (Split):\n")
                for bt in biotype_counts_split.index:
                    count = biotype_counts_split[bt]
                    total_bp = bp_by_biotype_split.get(bt, 0)
                    ratio_biotype = (total_bp / genome_size_split * 100) if genome_size_split else 0
                    f.write(f"  {bt}: {count} (Total bp: {total_bp}, Ratio: {ratio_biotype:.4f}%)\n")
                f.write("\nAverage gene length for non-protein-coding genes (Split):\n")
                non_pc_split = genes_df_split_copy[genes_df_split_copy["gene_biotype"] != "protein_coding"]
                if not non_pc_split.empty:
                    for bt in biotype_counts_split.index:
                        if bt == "protein_coding":
                            continue
                        subset = non_pc_split[non_pc_split["gene_biotype"] == bt]
                        if not subset.empty:
                            avg_len = subset["gene_length"].mean()
                            f.write(f"  {bt}: {avg_len:.4f} bp\n")
                else:
                    f.write("  No non-protein-coding genes available.\n")
        else:
            pass

        # ---- Non-gene Feature Coverage Section (Pseudogene and Biological Region) ----
        f.write("\n-------------------------\n")
        f.write("---- Non-gene Feature Coverage ----\n")
        # Pseudogene coverage (only for non-split scaffolds)
        pseudos_non_split = pseudos_df_included if split_set else pseudos_df
        if not pseudos_non_split.empty:
            total_bp_pseudogene = calculate_total_interval_bp(pseudos_non_split)
            ratio_pseudogene = (total_bp_pseudogene / genome_size_included * 100) if genome_size_included else 0
            avg_length_pseudogene = pseudos_non_split["end"].sub(pseudos_non_split["start"]).abs().mean() if not pseudos_non_split.empty else 0
            f.write("Pseudogene Coverage:\n")
            f.write(f"  Count: {pseudos_non_split.shape[0]}\n")
            f.write(f"  Total bp: {total_bp_pseudogene}\n")
            f.write(f"  Ratio: {ratio_pseudogene:.4f}%\n")
            f.write(f"  Average length: {avg_length_pseudogene:.4f} bp\n")
        else:
            f.write("Pseudogene Coverage: No pseudogenes found\n")
        
        # Biological_region coverage (separated because these features are not genes)
        bio_df = df[df['type'] == "biological_region"].copy()
        if split_set:
            bio_df_included = bio_df[~bio_df["seqid"].isin(split_set)].copy()
        else:
            bio_df_included = bio_df
        if not bio_df_included.empty:
            total_bp_bio = calculate_total_interval_bp(bio_df_included)
            ratio_bio = (total_bp_bio / genome_size_included * 100) if genome_size_included else 0
            avg_length_bio = bio_df_included["end"].sub(bio_df_included["start"]).abs().mean() if not bio_df_included.empty else 0
            f.write("\nBiological Region Coverage:\n")
            f.write(f"  Count: {bio_df_included.shape[0]}\n")
            f.write(f"  Total bp: {total_bp_bio}\n")
            f.write(f"  Ratio: {ratio_bio:.4f}%\n")
            f.write(f"  Average length: {avg_length_bio:.4f} bp\n")
        else:
            f.write("\nBiological Region Coverage: No biological_region features found\n")
        f.write("-------------------------\n")
        
        if args.protein:
            try:
                with open(args.protein, 'r') as pf:
                    protein_count = sum(1 for line in pf if line.startswith('>'))
                f.write(f"\nProtein FASTA total sequences: {protein_count}\n")
            except Exception as e:
                f.write(f"\nError reading protein FASTA file: {e}\n")
        f.write("\n==== END SUMMARY ====\n")
    # End of summary block

    # BUSCO and OMArk sections remain unchanged.
    busco_basename = args.basename
    if (args.protein and args.busco_lineage and args.threads is not None and args.busco_out_path):
        with open(args.protein, 'r') as pf:
            protein_count_dummy = sum(1 for line in pf if line.startswith('>'))
        print("Running BUSCO analysis using Docker...")
        protein_abs = os.path.abspath(args.protein)
        protein_dir = os.path.dirname(protein_abs)
        host_busco_out = os.path.abspath(args.busco_out_path)
        if args.busco_workdir:
            mount_dir = os.path.abspath(args.busco_workdir)
            container_busco_out = os.path.basename(host_busco_out)
        else:
            mount_dir = os.path.dirname(host_busco_out)
            container_busco_out = os.path.basename(host_busco_out)
        full_host_busco_out = os.path.join(mount_dir, container_busco_out)
        os.makedirs(full_host_busco_out, exist_ok=True)
        docker_cmd = [
            "docker", "run", "--rm",
            "-v", f"{os.path.abspath(protein_dir)}:/data",
            "-v", f"{mount_dir}:/output",
            "-u", "root",
            "-w", "/output",
            args.busco_docker_image,
            "busco",
            "-i", f"/data/{os.path.basename(protein_abs)}",
            "-m", "proteins",
            "-l", args.busco_lineage,
            "-c", str(args.threads),
            "-o", busco_basename,
            "--out_path", f"/output/{container_busco_out}"
        ]
        result_busco = subprocess.run(docker_cmd, check=True, capture_output=True, text=True)
        with open(qc_extra_file, 'a') as f:
            f.write("\n=== BUSCO Analysis ===\n")
            f.write(result_busco.stdout)
        print("BUSCO analysis completed.")
    if (args.omark_dir and args.protein and args.threads is not None):
        print("Running OMArk analysis using conda run...")
        script_dir = os.path.abspath(os.path.dirname(__file__))
        run_omark_path = os.path.join(script_dir, "run_omark.sh")
        omark_cmd = [
            "conda", "run", "-n", "omark_omamer", "bash", "-c",
            f"{run_omark_path} -i {args.protein} -base {busco_basename} -d {args.omark_dir}/{busco_basename} -t {args.threads}"
        ]
        result_omark = subprocess.run(omark_cmd, check=True, capture_output=True, text=True)
        with open(qc_extra_file, 'a') as f:
            f.write("\n=== OMArk Analysis ===\n")
            f.write(result_omark.stdout)
        print("OMArk analysis completed.")

    print(f"Done. Additional checks in '{qc_extra_file}'")

def run_interpro(args):
    import re
    gff_map = None
    if getattr(args, 'gff', None):
        unresolved_out = f"{args.outprefix}_gff_unresolved_protein_ids.tsv"
        gff_map = build_protein_to_gene_map_from_gff(args.gff, unresolved_out=unresolved_out)
    fallback_gene = lambda pid: re.sub(r"\.p\d+$", "", pid)
    proteins = set()
    genes = set()
    with open(args.protein, "r") as f:
        for line in f:
            if line.startswith(">"):
                protein_id = line.strip().split()[0][1:]
                proteins.add(protein_id)
                gene_id = map_protein_to_gene(protein_id, gff_map, fallback_gene)
                genes.add(gene_id)
    annotated_proteins_ipr = set()
    annotated_genes_ipr    = set()
    annotated_proteins_go  = set()
    annotated_genes_go     = set()
    with open(args.interpro_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            protein_name = cols[0].strip()
            gene_name = map_protein_to_gene(protein_name, gff_map, fallback_gene)
            found_ipr = any(col.strip().startswith("IPR") for col in cols[1:] if col.strip())
            found_go = any("GO:" in col for col in cols[1:] if col.strip())
            if found_ipr:
                annotated_proteins_ipr.add(protein_name)
                annotated_genes_ipr.add(gene_name)
            if found_go:
                annotated_proteins_go.add(protein_name)
                annotated_genes_go.add(gene_name)
    num_proteins = len(proteins)
    num_genes    = len(genes)
    num_ipr_prot = len(annotated_proteins_ipr)
    num_ipr_gene = len(annotated_genes_ipr)
    num_go_prot  = len(annotated_proteins_go)
    num_go_gene  = len(annotated_genes_go)
    ipr_prot_percent = (num_ipr_prot / num_proteins * 100) if num_proteins else 0
    ipr_gene_percent = (num_ipr_gene / num_genes * 100) if num_genes else 0
    go_prot_percent  = (num_go_prot  / num_proteins * 100) if num_proteins else 0
    go_gene_percent  = (num_go_gene  / num_genes * 100) if num_genes else 0
    out_prefix   = args.outprefix
    prot_ipr_out = f"{out_prefix}_annotated_protein_ipr.tsv"
    gene_ipr_out = f"{out_prefix}_annotated_gene_ipr.tsv"
    prot_go_out  = f"{out_prefix}_annotated_protein_go.tsv"
    gene_go_out  = f"{out_prefix}_annotated_gene_go.tsv"
    with open(prot_ipr_out, "w") as f_out:
        for p in sorted(annotated_proteins_ipr):
            f_out.write(f"{p}\tIPR annotated\n")
    with open(gene_ipr_out, "w") as f_out:
        for g in sorted(annotated_genes_ipr):
            f_out.write(f"{g}\tIPR annotated\n")
    with open(prot_go_out, "w") as f_out:
        for p in sorted(annotated_proteins_go):
            f_out.write(f"{p}\tGO annotated\n")
    with open(gene_go_out, "w") as f_out:
        for g in sorted(annotated_genes_go):
            f_out.write(f"{g}\tGO annotated\n")
    summary_text = (
        "[Summary]\n"
        f"Total proteins: {num_proteins}\n"
        f"Total genes:    {num_genes}\n\n"
        "[IPR Annotation]\n"
        f"Proteins with IPR: {num_ipr_prot} ({ipr_prot_percent:.4f}%)\n"
        f"Genes with IPR:    {num_ipr_gene} ({ipr_gene_percent:.4f}%)\n\n"
        "[GO Annotation]\n"
        f"Proteins with GO:  {num_go_prot} ({go_prot_percent:.4f}%)\n"
        f"Genes with GO:     {num_go_gene} ({go_gene_percent:.4f}%)\n\n"
        "[Output Files]\n"
        f"IPR-annotated proteins: {prot_ipr_out}\n"
        f"IPR-annotated genes:    {gene_ipr_out}\n"
        f"GO-annotated proteins:  {prot_go_out}\n"
        f"GO-annotated genes:     {gene_go_out}\n"
    )
    print(summary_text)
    summary_file = f"{out_prefix}_summary.txt"
    with open(summary_file, "w") as f_out:
        f_out.write(summary_text)

def run_kaas(args):
    file_path = args.query
    out_dir = args.outdir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    gff_map = None
    if getattr(args, 'gff', None):
        unresolved_out = os.path.join(out_dir, "kaas_gff_unresolved_protein_ids.tsv")
        gff_map = build_protein_to_gene_map_from_gff(args.gff, unresolved_out=unresolved_out)

    fallback_gene = lambda pid: pid.split(".", 1)[0]
    f1_count = 0
    f2_count = 0
    all_genes = set()
    genes_with_f2 = set()
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")
            if len(fields) >= 1 and fields[0].strip():
                f1_count += 1
                protein_name = fields[0].strip()
                gene_name = map_protein_to_gene(protein_name, gff_map, fallback_gene)
                all_genes.add(gene_name)
                if len(fields) >= 2 and fields[1].strip():
                    f2_count += 1
                    genes_with_f2.add(gene_name)
    ratio_lines = (f2_count / f1_count * 100) if f1_count else 0
    total_genes = len(all_genes)
    annotated_genes = len(genes_with_f2)
    ratio_genes = (annotated_genes / total_genes * 100) if total_genes else 0
    summary_lines = []
    summary_lines.append(f"File: {file_path}")
    summary_lines.append(f"Protein sequences:               {f1_count}")
    summary_lines.append(f"KAAS annotated Protein seq:      {f2_count} ({ratio_lines:.4f}%)")
    summary_lines.append(f"Total unique genes (first col):  {total_genes}")
    summary_lines.append(f"Annotated unique genes:          {annotated_genes} ({ratio_genes:.4f}%)")
    for line in summary_lines:
        print(line)
    gene_out_file = os.path.join(out_dir, "kaas_annotated_genes.tsv")
    with open(gene_out_file, "w", encoding="utf-8") as out_gene:
        for g in sorted(genes_with_f2):
            out_gene.write(f"{g}\tKO annotated\n")
    print(f"Annotated gene list saved to: {gene_out_file}")
    summary_file = os.path.join(out_dir, "kaas_summary.txt")
    with open(summary_file, "w", encoding="utf-8") as sf:
        sf.write("\n".join(summary_lines))
    print(f"Summary file saved to: {summary_file}")

def run_nr(args):
    gff_map = None
    if getattr(args, 'gff', None):
        unresolved_out = f"{args.outprefix}_gff_unresolved_protein_ids.tsv"
        gff_map = build_protein_to_gene_map_from_gff(args.gff, unresolved_out=unresolved_out)
    fallback_gene = lambda pid: re.sub(r"\.p\d+$", "", pid)
    proteins = set()
    genes = set()
    with open(args.protein, "r") as f:
        for line in f:
            if line.startswith(">"):
                protein_id = line.strip().split()[0][1:]
                proteins.add(protein_id)
                gene_id = map_protein_to_gene(protein_id, gff_map, fallback_gene)
                genes.add(gene_id)
    hits_protein = set()
    hits_gene = set()
    unchar_proteins = set()
    hypot_proteins = set()
    low_quality_proteins = set()
    unnamed_proteins = set()
    unchar_genes = set()
    hypot_genes = set()
    low_quality_genes = set()
    unnamed_genes = set()
    species_counter = Counter()
    with open(args.diamond_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("qseqid"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            qseqid = cols[0].strip()
            stitle = cols[1].strip()
            hits_protein.add(qseqid)
            gene_id = map_protein_to_gene(qseqid, gff_map, fallback_gene)
            hits_gene.add(gene_id)
            stitle_lower = stitle.lower()
            if "uncharacterized" in stitle_lower:
                unchar_proteins.add(qseqid)
                unchar_genes.add(gene_id)
            if "hypothetical protein" in stitle_lower:
                hypot_proteins.add(qseqid)
                hypot_genes.add(gene_id)
            if "low quality protein" in stitle_lower:
                low_quality_proteins.add(qseqid)
                low_quality_genes.add(gene_id)
            if "unnamed protein product" in stitle_lower:
                unnamed_proteins.add(qseqid)
                unnamed_genes.add(gene_id)
            match = re.search(r"\[(.*?)\]", stitle)
            if match:
                species = match.group(1).strip()
                species_counter[species] += 1
    num_proteins = len(proteins)
    num_hits_protein = len(hits_protein)
    pct_protein_hits = (num_hits_protein / num_proteins * 100) if num_proteins else 0
    num_genes = len(genes)
    num_hits_gene = len(hits_gene)
    pct_gene_hits = (num_hits_gene / num_genes * 100) if num_genes else 0
    def pct(count, denom):
        return (count / denom * 100) if denom else 0
    unchar_count_p = len(unchar_proteins)
    hypot_count_p = len(hypot_proteins)
    lowq_count_p = len(low_quality_proteins)
    unnamed_count_p = len(unnamed_proteins)
    unchar_pct_all_p = pct(unchar_count_p, num_proteins)
    hypot_pct_all_p = pct(hypot_count_p, num_proteins)
    lowq_pct_all_p = pct(lowq_count_p, num_proteins)
    unnamed_pct_all_p = pct(unnamed_count_p, num_proteins)
    unchar_pct_hits_p = pct(unchar_count_p, num_hits_protein)
    hypot_pct_hits_p = pct(hypot_count_p, num_hits_protein)
    lowq_pct_hits_p = pct(lowq_count_p, num_hits_protein)
    unnamed_pct_hits_p = pct(unnamed_count_p, num_hits_protein)
    unchar_count_g = len(unchar_genes)
    hypot_count_g = len(hypot_genes)
    lowq_count_g = len(low_quality_genes)
    unnamed_count_g = len(unnamed_genes)
    unchar_pct_all_g = pct(unchar_count_g, num_genes)
    hypot_pct_all_g = pct(hypot_count_g, num_genes)
    lowq_pct_all_g = pct(lowq_count_g, num_genes)
    unnamed_pct_all_g = pct(unnamed_count_g, num_genes)
    unchar_pct_hits_g = pct(unchar_count_g, num_hits_gene)
    hypot_pct_hits_g = pct(hypot_count_g, num_hits_gene)
    lowq_pct_hits_g = pct(lowq_count_g, num_hits_gene)
    unnamed_pct_hits_g = pct(unnamed_count_g, num_hits_gene)
    output_lines = []
    output_lines.append(f"Total proteins in FASTA: {num_proteins}")
    output_lines.append(f"Proteins with Diamond BLASTP hits: {num_hits_protein}")
    output_lines.append(f"Protein-level annotation percentage: {pct_protein_hits:.4f}%")
    output_lines.append("")
    output_lines.append(f"Total genes in FASTA: {num_genes}")
    output_lines.append(f"Genes with Diamond BLASTP hits: {num_hits_gene}")
    output_lines.append(f"Gene-level annotation percentage: {pct_gene_hits:.4f}%")
    output_lines.append("")
    output_lines.append("[Uncharacterized / Hypothetical / Low-quality / Unnamed] (Protein level)")
    output_lines.append("Counts among all proteins, plus two percentages:")
    output_lines.append("  1) % among ALL proteins")
    output_lines.append("  2) % among proteins WITH hits")
    output_lines.append("")
    output_lines.append("Uncharacterized:")
    output_lines.append(f"  Count: {unchar_count_p}")
    output_lines.append(f"  Among ALL proteins: {unchar_pct_all_p:.4f}%")
    output_lines.append(f"  Among proteins WITH hits: {unchar_pct_hits_p:.4f}%")
    output_lines.append("")
    output_lines.append("Hypothetical protein:")
    output_lines.append(f"  Count: {hypot_count_p}")
    output_lines.append(f"  Among ALL proteins: {hypot_pct_all_p:.4f}%")
    output_lines.append(f"  Among proteins WITH hits: {hypot_pct_hits_p:.4f}%")
    output_lines.append("")
    output_lines.append("Low quality protein:")
    output_lines.append(f"  Count: {lowq_count_p}")
    output_lines.append(f"  Among ALL proteins: {lowq_pct_all_p:.4f}%")
    output_lines.append(f"  Among proteins WITH hits: {lowq_pct_hits_p:.4f}%")
    output_lines.append("")
    output_lines.append("Unnamed protein product:")
    output_lines.append(f"  Count: {unnamed_count_p}")
    output_lines.append(f"  Among ALL proteins: {unnamed_pct_all_p:.4f}%")
    output_lines.append(f"  Among proteins WITH hits: {unnamed_pct_hits_p:.4f}%")
    output_lines.append("")
    output_lines.append("[Uncharacterized / Hypothetical / Low-quality / Unnamed] (Gene level)")
    output_lines.append("Counts among all genes, plus two percentages:")
    output_lines.append("  1) % among ALL genes")
    output_lines.append("  2) % among genes WITH hits")
    output_lines.append("")
    output_lines.append("Uncharacterized genes:")
    output_lines.append(f"  Count: {unchar_count_g}")
    output_lines.append(f"  Among ALL genes: {unchar_pct_all_g:.4f}%")
    output_lines.append(f"  Among genes WITH hits: {unchar_pct_hits_g:.4f}%")
    output_lines.append("")
    output_lines.append("Hypothetical genes:")
    output_lines.append(f"  Count: {hypot_count_g}")
    output_lines.append(f"  Among ALL genes: {hypot_pct_all_g:.4f}%")
    output_lines.append(f"  Among genes WITH hits: {hypot_pct_hits_g:.4f}%")
    output_lines.append("")
    output_lines.append("Low-quality genes:")
    output_lines.append(f"  Count: {lowq_count_g}")
    output_lines.append(f"  Among ALL genes: {lowq_pct_all_g:.4f}%")
    output_lines.append(f"  Among genes WITH hits: {lowq_pct_hits_g:.4f}%")
    output_lines.append("")
    output_lines.append("Unnamed-product genes:")
    output_lines.append(f"  Count: {unnamed_count_g}")
    output_lines.append(f"  Among ALL genes: {unnamed_pct_all_g:.4f}%")
    output_lines.append(f"  Among genes WITH hits: {unnamed_pct_hits_g:.4f}%")
    output_lines.append("")
    output_lines.append("[Informative Ratio Among Proteins With Hits]")
    informative_hits_p = hits_protein - (unchar_proteins | hypot_proteins | low_quality_proteins | unnamed_proteins)
    count_informative_p = len(informative_hits_p)
    ratio_informative_all_p = pct(count_informative_p, num_proteins)
    ratio_informative_hits_p = pct(count_informative_p, num_hits_protein)
    output_lines.append(f"  Informative proteins (hits): {count_informative_p}")
    output_lines.append(f"  Among ALL proteins: {ratio_informative_all_p:.4f}%")
    output_lines.append(f"  Among proteins WITH hits: {ratio_informative_hits_p:.4f}%")
    output_lines.append("")
    output_lines.append("[Informative Ratio Among Genes With Hits]")
    informative_hits_g = hits_gene - (unchar_genes | hypot_genes | low_quality_genes | unnamed_genes)
    count_informative_g = len(informative_hits_g)
    ratio_informative_all_g = pct(count_informative_g, num_genes)
    ratio_informative_hits_g = pct(count_informative_g, num_hits_gene)
    output_lines.append(f"  Informative genes (hits): {count_informative_g}")
    output_lines.append(f"  Among ALL genes: {ratio_informative_all_g:.4f}%")
    output_lines.append(f"  Among genes WITH hits: {ratio_informative_hits_g:.4f}%")
    output_lines.append("")
    informative_genes_file = f"{args.outprefix}.informative_gene_list.tsv"
    output_lines.append(f"Informative gene list written to: {informative_genes_file}")
    output_lines.append("---------------------------")
    if species_counter:
        most_common_3 = species_counter.most_common(3)
        top_three_sum = sum(count for spec, count in most_common_3)
        top_three_ratio_all_p = pct(top_three_sum, num_proteins)
        top_three_ratio_hits_p = pct(top_three_sum, num_hits_protein)
        output_lines.append(f"Top three species total hits: {top_three_sum}")
        output_lines.append(f"  Among ALL proteins: {top_three_ratio_all_p:.4f}%")
        output_lines.append(f"  Among proteins WITH hits:  {top_three_ratio_hits_p:.4f}%")
        output_lines.append("")

    output_lines.append("Top 3 species among hits (and ratio among ALL proteins):")
    if species_counter:
        for i, (spec, count) in enumerate(species_counter.most_common(3), start=1):
            species_ratio_all = pct(count, num_proteins)
            output_lines.append(f"{i}. {spec} ({count} hits) => {species_ratio_all:.4f}% of all proteins")
    else:
        output_lines.append("No species found in annotation (check your input).")
    summary_text = "\n".join(output_lines)
    summary_filename = f"{args.outprefix}_diamond_nr_summary.txt"
    with open(summary_filename, "w") as out_file:
        out_file.write(summary_text)
    with open(informative_genes_file, "w") as out_g:
        for g in sorted(informative_hits_g):
            out_g.write(f"{g}\tInformative diamond blastp\n")
    print(summary_text)
    print(f"\nDetailed summary written to: {summary_filename}")

def main():
    parser = argparse.ArgumentParser(
        description="Combined script for Structural Annotation, InterPro, KAAS, and Diamond BLASTP (NR) processing."
    )
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    structural_parser = subparsers.add_parser("structural",
        help="Structural Annotation mode (process GFF QC, gene biotype bp and average length, and non-coding/pseudogene coverage).")
    structural_parser.add_argument("-g", "--gff", required=True, help="Input GFF file")
    structural_parser.add_argument("-f", "--fasta", required=True, help="Input genome FASTA file (.fna)")
    structural_parser.add_argument("-b", "--basename", required=True, help="Base name for output files")
    structural_parser.add_argument("-o", "--outdir", default=".", help="Output directory for Structural mode files (default: current directory)")
    structural_parser.add_argument("-p", "--protein", help="Protein sequence provided for BUSCO (ex. braker.aa path)")
    structural_parser.add_argument("--busco_lineage", help="BUSCO lineage (e.g. busco lineage)")
    structural_parser.add_argument("-t", "--threads", type=int, help="Number of threads for BUSCO and OMArk")
    structural_parser.add_argument("--busco_out_path", help="Full host path for BUSCO results output (this directory will be mounted under /output in the container)")
    structural_parser.add_argument("--busco_docker_image", default="busco/busco:latest", help="Docker image for BUSCO analysis (default: busco/busco:latest)")
    structural_parser.add_argument("--busco_workdir", help="Optional: Docker working directory to mount as /output. If provided, it is used as the container's working directory.")
    structural_parser.add_argument("--omark_dir", help="Directory for OMArk run output (for -d option)")
    structural_parser.add_argument("--split_scaffolds", nargs="+",
        help="List of scaffolds/chromosomes to split out and report separately in summary statistics (e.g. mitochondria).")
    structural_parser.set_defaults(func=run_structural)

    ip = subparsers.add_parser("interpro", help="Process an InterProScan TSV file.")
    ip.add_argument("-p", "--protein", required=True, help="Input protein FASTA file")
    ip.add_argument("-i", "--interpro_tsv", required=True, help="InterProScan TSV file")
    ip.add_argument("-g", "--gff", help="Optional GFF3 file to map protein IDs to gene IDs via CDS->mRNA->gene (NCBI/RefSeq style).")
    ip.add_argument("-o", "--outprefix", default="interpro_output", help="Prefix for output annotation files.")
    ip.set_defaults(func=run_interpro)

    ka = subparsers.add_parser("kaas",
        help="Process a KAAS result file and generate a summary file. Usage: python eva_annotation.py kaas query.ko [outdir]")
    ka.add_argument("query", help="KAAS result file (e.g., query.ko)")
    ka.add_argument("outdir", nargs="?", default=".", help="Output directory (default: current directory)")
    ka.add_argument("-g", "--gff", help="Optional GFF3 file to map first-column protein IDs to gene IDs via CDS->mRNA->gene (NCBI/RefSeq style).")
    ka.set_defaults(func=run_kaas)

    dn = subparsers.add_parser("diamond-nr", help="Evaluate Diamond BLASTP TSV output (NR database mode).")
    dn.add_argument("-p", "--protein", required=True, help="Path to the input protein FASTA file.")
    dn.add_argument("-d", "--diamond_tsv", required=True, help="Path to the Diamond BLASTP TSV output file.")
    dn.add_argument("-g", "--gff", help="Optional GFF3 file to map protein IDs to gene IDs via CDS->mRNA->gene (NCBI/RefSeq style).")
    dn.add_argument("-o", "--outprefix", default="output", help="Prefix for output files (default: 'output').")
    dn.set_defaults(func=run_nr)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
