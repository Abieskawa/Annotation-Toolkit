#!/usr/bin/env python3
import os
import re
import sys
import subprocess
import argparse
import statistics
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------
# 1) Utility: Parse GFF attributes
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

# -----------------------------------------------------------------------
# 2) Utility: Plot a length histogram (for genes or sequences)
# -----------------------------------------------------------------------
def plot_bar(lengths, med, out_file, bin_size, thresh, xlabel):
    """
    Plots a histogram of 'lengths' with bins of width 'bin_size' up to 'thresh',
    and one extra bin for values >= thresh.
    
    For two specific cases, the last bin is labeled specially:
      - If bin_size==2000 and thresh==30000, the last bin is labeled ">30000"
      - If bin_size==100  and thresh==2000,  the last bin is labeled ">2000"
      
    A vertical dashed red line is drawn at the median.
    The plot is saved to 'out_file'.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    if not lengths:
        print("No lengths to plot.")
        return

    # Create bins from 0 to thresh (non-inclusive) and then a final bin covering [thresh, +inf)
    bins = list(range(0, thresh, bin_size)) + [float('inf')]
    counts, edges = np.histogram(lengths, bins=bins)

    centers = []
    widths = []
    for i in range(len(edges) - 1):
        left_edge, right_edge = edges[i], edges[i+1]
        if not np.isinf(right_edge):
            centers.append((left_edge + right_edge) / 2)
            widths.append(right_edge - left_edge)
        else:
            # For the last bin ([thresh, +inf)), use a center offset by half the bin size.
            centers.append(thresh + bin_size / 2)
            widths.append(bin_size)

    plt.figure(figsize=(8,6))
    plt.bar(centers, counts, width=widths, align='center', alpha=0.7, edgecolor='black')

    # If the median exceeds 'thresh', clamp it to appear on the last bin
    med_x = min(med, thresh + bin_size/2)
    plt.axvline(med_x, color='red', linestyle='dashed', linewidth=2, label=f"Median: {int(med)}")

    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.title("Length Distribution")
    plt.legend()

    # Create tick labels from the centers; then override the last label if in one of the special cases.
    tick_labels = [f"{int(c)}" for c in centers]
    if bin_size == 2000 and thresh == 30000:
        tick_labels[-1] = ">30000"
    elif bin_size == 100 and thresh == 2000:
        tick_labels[-1] = ">2000"
        
    plt.xticks(centers, tick_labels)

    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()

# -----------------------------------------------------------------------
# 3) Utility: Check start/stop codons for an mRNA
# -----------------------------------------------------------------------
def check_codon(mrna, grp, codon_type):
    """
    Return an error message if the mRNA lacks exactly one codon of 'codon_type'.
    We check the start/stop codons that are within the mRNA's [start, end],
    same seqid, same strand, same Parent ID.
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
# 4) Sub-command: GFF Processing + QC
# -----------------------------------------------------------------------
def run_gff(args):
    # Read GFF file
    cols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    df = pd.read_csv(args.gff, sep='\t', comment='#', names=cols, header=None,
                     dtype={'start': int, 'end': int}, na_filter=False)
    df = df[(df['seqid'] != '') & (df['type'] != '')]

    # Parse attributes
    df['parsed_attrs'] = df['attributes'].apply(parse_attributes)
    df['ID']     = df['parsed_attrs'].apply(lambda d: d.get('ID'))
    df['Parent'] = df['parsed_attrs'].apply(lambda d: d.get('Parent'))

    # Identify gene/pseudogene features
    genes_df     = df[df['type'] == "gene"].copy()
    pseudos_df   = df[df['type'] == "pseudogene"].copy()
    n_genes      = genes_df.shape[0]
    n_pseudos    = pseudos_df.shape[0]

    # Possibly parse gene_biotype if present
    gene_biotype_available = False
    if not genes_df.empty:
        # Extract gene_biotype, defaulting to "NA"
        genes_df["gene_biotype"] = genes_df["parsed_attrs"].apply(lambda d: d.get("gene_biotype", "NA"))
        # If everything is "NA", we treat it as unavailable
        unique_biotypes = set(genes_df["gene_biotype"].unique())
        if unique_biotypes != {"NA"}:
            gene_biotype_available = True

    # Skip codon checks if:
    #  1) gene_biotype info is available, OR
    #  2) 'Gnomon' found in the source column
    skip_codon_check = False
    if gene_biotype_available or ("Gnomon" in df["source"].unique()):
        skip_codon_check = True

    # Start/Stop Codon Checks (only if not skipping)
    start_errs = []
    stop_errs  = []

    mRNA_df = df[df['type'] == 'mRNA'].copy()
    if not skip_codon_check:
        if gene_biotype_available:
            # Merge gene info to see which mRNAs are from protein-coding genes
            gene_info = genes_df[['ID', 'gene_biotype']].rename(columns={'ID': 'gene'})
            mRNA_df = mRNA_df.merge(gene_info, left_on='Parent', right_on='gene', how='left')
        start_codons = df[df['type'] == 'start_codon']
        stop_codons  = df[df['type'] == 'stop_codon']
        start_grp = start_codons.groupby("Parent") if not start_codons.empty else None
        stop_grp  = stop_codons.groupby("Parent")  if not stop_codons.empty  else None

        if gene_biotype_available:
            # Filter to only protein_coding mRNAs
            pc_mrnas = mRNA_df[mRNA_df["gene_biotype"] == "protein_coding"]
            for _, r in pc_mrnas.iterrows():
                err = check_codon(r, start_grp, "start")
                if err:
                    start_errs.append(err)
                err = check_codon(r, stop_grp, "stop")
                if err:
                    stop_errs.append(err)
        else:
            # Check for all mRNAs
            for _, r in mRNA_df.iterrows():
                err = check_codon(r, start_grp, "start")
                if err:
                    start_errs.append(err)
                err = check_codon(r, stop_grp, "stop")
                if err:
                    stop_errs.append(err)

    # When calculating total exon number and single-exon genes,
    # refer to protein-coding genes if gene_biotype is available.
    exon_df = df[df['type'] == "exon"].copy()
    if gene_biotype_available:
        # 1) Identify protein-coding mRNAs
        gene_info = genes_df[['ID', 'gene_biotype']].rename(columns={'ID': 'gene'})
        mRNA_biotype = df[df['type'] == 'mRNA'].merge(gene_info, left_on='Parent', right_on='gene', how='left')
        pc_mRNAs = mRNA_biotype[mRNA_biotype['gene_biotype'] == 'protein_coding'].copy()
        # 2) Filter exons to only those belonging to protein-coding mRNAs
        pc_exons = exon_df.merge(pc_mRNAs[['ID']], left_on='Parent', right_on='ID', how='inner')
        # Count total exons (protein-coding only)
        n_exon = pc_exons.shape[0]
        # 3) Count single-exon genes (protein-coding only)
        pc_mRNAs_map = pc_mRNAs[['ID', 'Parent']].rename(columns={'Parent': 'Gene_ID'})
        pc_exons_merged = pc_exons.merge(pc_mRNAs_map, left_on='Parent', right_on='ID', how='left', suffixes=('_exon','_mRNA'))
        pc_exon_counts = pc_exons_merged.groupby('Gene_ID').size()
        single_exon_genes = 0
        pc_genes = genes_df[genes_df["gene_biotype"] == "protein_coding"]
        for gid in pc_genes['ID']:
            if pc_exon_counts.get(gid, 0) == 1:
                single_exon_genes += 1
    else:
        # Old method if gene_biotype not available
        n_exon = exon_df.shape[0]
        mRNA_gene_map = df[df['type'] == 'mRNA'][['ID','Parent']].rename(columns={'Parent':'Gene_ID'})
        exons_merged = exon_df.merge(mRNA_gene_map, left_on='Parent', right_on='ID', how='left', suffixes=('_exon','_mRNA'))
        exon_counts = exons_merged.groupby('Gene_ID').size()

        # Consider only genes (excluding pseudogenes) for single-exon gene calculation
        all_genes = genes_df.copy()
        single_exon_genes = 0
        for gid in all_genes['ID']:
            if exon_counts.get(gid, 0) == 1:
                single_exon_genes += 1

    # Count mRNAs in the final mRNA_df (may be filtered if codon checks were applied)
    n_mrna = mRNA_df.shape[0]

    # Compute gene lengths and median
    all_genes = genes_df.copy()
    if not all_genes.empty:
        all_genes["gene_length"] = (all_genes["end"] - all_genes["start"]).abs()
        gene_lengths = all_genes["gene_length"].tolist()
        gene_lengths.sort()
        if gene_lengths:
            mid = len(gene_lengths) // 2
            if len(gene_lengths) % 2 == 1:
                median_len = gene_lengths[mid]
            else:
                median_len = (gene_lengths[mid-1] + gene_lengths[mid]) / 2
        else:
            median_len = 0
    else:
        gene_lengths = []
        median_len = 0

    # Additional stats if gene_biotype is available
    biotype_counts = None
    if gene_biotype_available:
        biotype_counts = genes_df["gene_biotype"].value_counts()
        pc_genes = genes_df[genes_df["gene_biotype"] == "protein_coding"]
        npc_genes = genes_df[genes_df["gene_biotype"] != "protein_coding"]
        lengths_pc = (pc_genes["end"] - pc_genes["start"]).abs().tolist() if not pc_genes.empty else []
        med_len_pc = int(statistics.median(lengths_pc)) if lengths_pc else 0
        lengths_npc = (npc_genes["end"] - npc_genes["start"]).abs().tolist() if not npc_genes.empty else []
        med_len_npc = int(statistics.median(lengths_npc)) if lengths_npc else 0
    else:
        med_len_pc = 0
        med_len_npc = 0

    # Run external gff3_QC command
    base = args.basename
    qc_output_file = f"{base}_qc.txt"
    qc_stat_file   = f"{base}_statistic.txt"
    qc_extra_file  = f"{base}_additional_validation.txt"
    qc_plot_file   = f"{base}_gene_length_histogram.png"

    subprocess.run([
        "gff3_QC",
        "--gff", args.gff,
        "--fasta", args.fasta,
        "--output", qc_output_file,
        "--statistic", qc_stat_file
    ], check=True)

    # Plot histogram if gene lengths exist (using 2000-bp bins with a final bin labeled ">30000")
    if gene_lengths:
        plot_bar(gene_lengths, median_len, qc_plot_file, bin_size=2000, thresh=30000, xlabel="Gene Length (bp)")

    # Write additional QC info to file
    with open(qc_extra_file, 'w') as f:
        f.write("=== Additional Validation Checks ===\n\n")
        if not skip_codon_check:
            f.write(f"Stop codon issues: {len(stop_errs)}\n")
            for err in stop_errs:
                f.write(err + "\n")
            f.write(f"\nStart codon issues: {len(start_errs)}\n")
            for err in start_errs:
                f.write(err + "\n")
        else:
            f.write("Skipping codon checks (gene_biotype present or source is 'Gnomon').\n")
        f.write("\n=== Summary Stats ===\n")
        f.write(f"Genes: {n_genes}")
        f.write(f"Pseudogenes: {n_pseudos})\n")
        f.write(f"Single-exon genes (Protein-coding if available): {single_exon_genes}\n")
        f.write(f"mRNAs: {n_mrna}\n")
        f.write(f"Exons (Protein-coding if available): {n_exon}\n")
        f.write(f"Median gene length (all genes): {median_len}\n")
        if gene_biotype_available and biotype_counts is not None:
            f.write("\nSeparate median lengths:\n")
            f.write(f"  Protein-coding:     {med_len_pc}\n")
            f.write(f"  Non-protein-coding: {med_len_npc}\n")
            f.write("\nGene biotype counts:\n")
            for bt, count in biotype_counts.items():
                f.write(f"  {bt}: {count}\n")

    print(f"Done. Additional checks in '{qc_extra_file}', histogram in '{qc_plot_file}'.")

# -----------------------------------------------------------------------
# 5) Sub-command: InterPro Processing
# -----------------------------------------------------------------------
def run_interpro(args):
    """
    Processes an InterProScan TSV file to summarize IPR and GO annotations,
    and writes out annotated protein and gene lists.
    """
    import re

    # 1. Collect protein and gene IDs from FASTA
    proteins = set()
    genes = set()
    with open(args.fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                protein_id = line.strip().split()[0][1:]
                proteins.add(protein_id)
                gene_id = re.sub(r"\.p\d+$", "", protein_id)
                genes.add(gene_id)

    # 2. Parse InterProScan TSV
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
            # Skip lines that don't have any columns beyond the protein name.
            if len(cols) < 2:
                continue

            # The first column is the protein name.
            protein_name = cols[0].strip()
            gene_name = re.sub(r"\.p\d+$", "", protein_name)

            # Look for an IPR field in the remaining columns.
            found_ipr = any(col.strip().startswith("IPR") for col in cols[1:] if col.strip())
            # Look for any GO term in any column (if the field contains "GO:")
            found_go = any("GO:" in col for col in cols[1:] if col.strip())

            if found_ipr:
                annotated_proteins_ipr.add(protein_name)
                annotated_genes_ipr.add(gene_name)
            if found_go:
                annotated_proteins_go.add(protein_name)
                annotated_genes_go.add(gene_name)
    
    # 3. Summary stats
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

    # 4. Write annotation outputs
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

    # 5. Print summary
    print("[Summary]")
    print(f"Total proteins: {num_proteins}")
    print(f"Total genes:    {num_genes}")
    print("\n[IPR Annotation]")
    print(f"Proteins with IPR: {num_ipr_prot} ({ipr_prot_percent:.2f}%)")
    print(f"Genes with IPR:    {num_ipr_gene} ({ipr_gene_percent:.2f}%)")
    print("\n[GO Annotation]")
    print(f"Proteins with GO:  {num_go_prot} ({go_prot_percent:.2f}%)")
    print(f"Genes with GO:     {num_go_gene} ({go_gene_percent:.2f}%)")
    print("\n[Output Files]")
    print(f"IPR-annotated proteins: {prot_ipr_out}")
    print(f"IPR-annotated genes:    {gene_ipr_out}")
    print(f"GO-annotated proteins:  {prot_go_out}")
    print(f"GO-annotated genes:     {gene_go_out}")

# -----------------------------------------------------------------------
# 6) Sub-command: KAAS Processing
# -----------------------------------------------------------------------
def run_kaas(args):
    """
    Processes a KAAS result file to calculate annotation stats
    and outputs a list of annotated genes.
    """
    file_path = args.file
    out_dir = args.outdir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    f1_count = 0  # Lines with non-empty first field
    f2_count = 0  # Lines with non-empty second field
    all_genes = set()
    genes_with_f2 = set()

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")
            if len(fields) >= 1 and fields[0].strip():
                f1_count += 1
                protein_name = fields[0].strip()
                gene_name = protein_name.split(".", 1)[0]
                all_genes.add(gene_name)
                if len(fields) >= 2 and fields[1].strip():
                    f2_count += 1
                    genes_with_f2.add(gene_name)

    ratio_lines = (f2_count / f1_count * 100) if f1_count else 0
    total_genes = len(all_genes)
    annotated_genes = len(genes_with_f2)
    ratio_genes = (annotated_genes / total_genes * 100) if total_genes else 0

    print(f"File: {file_path}")
    print(f"Protein sequences:               {f1_count}")
    print(f"KAAS annotated Protein seq:      {f2_count} ({ratio_lines:.2f}%)")
    print(f"Total unique genes (first col):  {total_genes}")
    print(f"Annotated unique genes:          {annotated_genes} ({ratio_genes:.2f}%)")

    gene_out_file = os.path.join(out_dir, "kaas_annotated_genes.tsv")
    with open(gene_out_file, "w", encoding="utf-8") as out_gene:
        for g in sorted(genes_with_f2):
            out_gene.write(f"{g}\tKO annotated\n")

    print(f"Annotated gene list saved to: {gene_out_file}")

# -----------------------------------------------------------------------
# 7) Sub-command: FASTA Histogram
# -----------------------------------------------------------------------
def run_hist(args):
    """
    Plots a length histogram from a FASTA file of DNA or protein sequences.
    """
    lengths = []
    seq = ""
    header = None

    with open(args.fasta, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    lengths.append(len(seq))
                header = line
                seq = ""
            else:
                seq += line
    if header is not None:
        lengths.append(len(seq))

    if not lengths:
        print("No sequences found in FASTA file. Exiting.")
        return

    median_len = statistics.median(lengths)
    ext = os.path.splitext(args.fasta)[1].lower()
    if ext == ".aa":
        default_bin, default_thresh, default_xlabel = (100, 2000, "Protein Length (aa)")
    else:
        default_bin, default_thresh, default_xlabel = (1000, 20000, "Sequence Length (bp)")

    bin_size = args.bin_size if args.bin_size is not None else default_bin
    thresh = args.threshold if args.threshold is not None else default_thresh
    xlabel = args.xlabel if args.xlabel is not None else default_xlabel

    # Plot the bar chart using the plot_bar function.
    plot_bar(lengths, median_len, args.output, bin_size, thresh, xlabel)
    print(f"Histogram saved to {args.output}")

# -----------------------------------------------------------------------
# Main Parser
# -----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Combined script for GFF QC, InterPro, KAAS, and histogram plotting.")
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    # 1. GFF sub-command
    pg = subparsers.add_parser("gff", help="Process a GFF file and run QC checks (including gene_biotype if available).")
    pg.add_argument("-g", "--gff", required=True, help="Input GFF file")
    pg.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    pg.add_argument("-b", "--basename", required=True, help="Base name for output files")
    pg.set_defaults(func=run_gff)

    # 2. InterPro sub-command
    pi = subparsers.add_parser("interpro", help="Process an InterProScan TSV (molas) file.")
    pi.add_argument("-f", "--fasta", required=True, help="Input protein FASTA file")
    pi.add_argument("-i", "--interpro_tsv", required=True, help="InterProScan TSV file")
    pi.add_argument("-o", "--outprefix", default="interpro_output", help="Prefix for output annotation files.")
    pi.set_defaults(func=run_interpro)

    # 3. KAAS sub-command
    pk = subparsers.add_parser("kaas", help="Process a KAAS result file.")
    pk.add_argument("file", help="Input KAAS result file")
    pk.add_argument("outdir", nargs="?", default=".", help="Output directory (default: current dir)")
    pk.set_defaults(func=run_kaas)

    # 4. Hist sub-command
    ph = subparsers.add_parser("hist", help="Plot a length histogram from a FASTA file.")
    ph.add_argument("-f", "--fasta", required=True, help="Input FASTA file (.fa, .fasta, or .aa)")
    ph.add_argument("-o", "--output", required=True, help="Output image file (PNG, etc.)")
    ph.add_argument("--bin_size", type=int, default=None, help="Bin size (default depends on file extension).")
    ph.add_argument("--threshold", type=int, default=None, help="Upper threshold for histogram.")
    ph.add_argument("--xlabel", type=str, default=None, help="Label for x-axis.")
    ph.set_defaults(func=run_hist)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
