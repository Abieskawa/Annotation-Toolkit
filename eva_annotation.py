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
# Protein Structural Annotation (PSA) Mode
# -----------------------------------------------------------------------
def run_psa(args):
    cols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    df = pd.read_csv(args.gff, sep='\t', comment='#', names=cols, header=None,
                     dtype={'start': int, 'end': int}, na_filter=False)
    df = df[(df['seqid'] != '') & (df['type'] != '')]

    df['parsed_attrs'] = df['attributes'].apply(parse_attributes)
    df['ID'] = df['parsed_attrs'].apply(lambda d: d.get('ID'))
    df['Parent'] = df['parsed_attrs'].apply(lambda d: d.get('Parent'))

    # Subset feature-specific dataframes and ensure unique entries
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
        unique_biotypes = set(genes_df["gene_biotype"].unique())
        if unique_biotypes != {"NA"}:
            gene_biotype_available = True

    skip_codon_check = False
    if gene_biotype_available or ("Gnomon" in df["source"].unique()):
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

    # Calculate single exon genes using unique mRNA and exon entries.
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

    # Process gene lengths using unique gene entries.
    all_genes = genes_df.copy()
    if not all_genes.empty:
        all_genes["gene_length"] = (all_genes["end"] - all_genes["start"]).abs()
        gene_lengths = all_genes["gene_length"].tolist()
        median_len = statistics.median(gene_lengths) if gene_lengths else 0
    else:
        gene_lengths = []
        median_len = 0

    biotype_counts = None
    if gene_biotype_available:
        biotype_counts = genes_df["gene_biotype"].value_counts()
        pc_genes = genes_df[genes_df["gene_biotype"] == "protein_coding"]
        npc_genes = genes_df[genes_df["gene_biotype"] != "protein_coding"]
        lengths_pc = (pc_genes["end"] - pc_genes["start"]).abs().tolist() if not pc_genes.empty else []
        med_len_pc = statistics.median(lengths_pc) if lengths_pc else 0
        lengths_npc = (npc_genes["end"] - npc_genes["start"]).abs().tolist() if not npc_genes.empty else []
        med_len_npc = statistics.median(lengths_npc) if lengths_npc else 0
    else:
        med_len_pc = 0
        med_len_npc = 0

    base = args.basename
    outdir = args.outdir
    qc_output_file = os.path.join(outdir, f"{base}_qc.txt")
    qc_stat_file   = os.path.join(outdir, f"{base}_statistic.txt")
    qc_extra_file  = os.path.join(outdir, f"{base}_additional_validation.txt")
    qc_plot_file   = os.path.join(outdir, f"{base}_gene_length_histogram.png")

    subprocess.run([
        "gff3_QC",
        "--gff", args.gff,
        "--fasta", args.fasta,
        "--output", qc_output_file,
        "--statistic", qc_stat_file
    ], check=True)

    if gene_lengths:
        plot_bar(gene_lengths, median_len, qc_plot_file, bin_size=2000, thresh=30000, xlabel="Gene Length (bp)")

    with open(qc_extra_file, 'w') as f:
        f.write("=== Additional Validation Checks ===\n\n")
        # Report any duplicate warnings for each feature type.
        if gene_dup_msg:
            f.write(gene_dup_msg)
        if pseudo_dup_msg:
            f.write(pseudo_dup_msg)
        if mrna_dup_msg:
            f.write(mrna_dup_msg)
        if exon_dup_msg:
            f.write(exon_dup_msg)
        f.write("\n")
        if not skip_codon_check:
            f.write(f"\nStart codon issues: {len(start_errs)}\n")
            f.write(f"Stop codon issues: {len(stop_errs)}\n")
            f.write(f"******Details******\n")
            f.write(f"Start codon issues\n")
            for err in start_errs:
                f.write(err + "\n")
            f.write(f"==================================\n")
            f.write(f"Stop codon issues\n")
            for err in stop_errs:
                f.write(err + "\n")
            f.write(f"==================================\n")
        else:
            f.write("Skipping codon checks (gene_biotype present or source is 'Gnomon').\n")
        f.write("\n=== Summary Stats ===\n")
        f.write(f"Genes (unique): {n_genes}\n")
        f.write(f"Pseudogenes (unique): {n_pseudos}\n")
        f.write(f"Single-exon genes (Protein-coding if available): {single_exon_genes}\n")
        f.write(f"mRNAs (unique): {n_mrna}\n")
        f.write(f"Exons (unique, Protein-coding if available): {n_exon}\n")
        f.write(f"Median gene length (all genes): {median_len}\n")
        if gene_biotype_available and biotype_counts is not None:
            f.write("\nSeparate median lengths:\n")
            f.write(f"  Protein-coding:     {med_len_pc}\n")
            f.write(f"  Non-protein-coding: {med_len_npc}\n")
            f.write("\nGene biotype counts:\n")
            for bt, count in biotype_counts.items():
                f.write(f"  {bt}: {count}\n")
        # If a protein FASTA file is provided, count the total sequences.
        if args.protein:
            try:
                with open(args.protein, 'r') as pf:
                    protein_count = sum(1 for line in pf if line.startswith('>'))
                f.write(f"\nProtein FASTA total sequences: {protein_count}\n")
            except Exception as e:
                f.write(f"\nError reading protein FASTA file: {e}\n")

    # Use the basename as the BUSCO basename
    busco_basename = args.basename

    # --- Begin BUSCO analysis via Docker (if BUSCO parameters are provided) ---
    if (args.protein and args.busco_lineage and args.threads is not None and args.busco_out_path):
        with open(args.protein, 'r') as pf:
            protein_count_dummy = sum(1 for line in pf if line.startswith('>'))
        print("Running BUSCO analysis using Docker...")
        protein_abs = os.path.abspath(args.protein)
        protein_dir = os.path.dirname(protein_abs)
        # Determine mount directory and container BUSCO output path:
        host_busco_out = os.path.abspath(args.busco_out_path)
        if args.busco_workdir:
            mount_dir = os.path.abspath(args.busco_workdir)
            container_busco_out = os.path.basename(host_busco_out)
        else:
            mount_dir = os.path.dirname(host_busco_out)
            container_busco_out = os.path.basename(host_busco_out)
        # Ensure the host BUSCO output directory exists.
        full_host_busco_out = os.path.join(mount_dir, container_busco_out)
        os.makedirs(full_host_busco_out, exist_ok=True)
        docker_cmd = [
            "docker", "run", "--rm",
            "-v", f"{os.path.abspath(protein_dir)}:/data",
            "-v", f"{mount_dir}:/output",
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
    # --- End BUSCO analysis ---

    # --- Begin OMArk analysis via conda (if omark parameters are provided) ---
    if (args.omark_dir and args.protein and args.threads is not None):
        print("Running OMArk analysis using conda run...")
        # Get the absolute path of run_omark.sh (assumed to be in the same directory as this script)
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
    # --- End OMArk analysis ---

    print(f"Done. Additional checks in '{qc_extra_file}', histogram in '{qc_plot_file}'.")

# -----------------------------------------------------------------------
# InterPro Processing
# -----------------------------------------------------------------------
def run_interpro(args):
    import re
    proteins = set()
    genes = set()
    with open(args.protein, "r") as f:
        for line in f:
            if line.startswith(">"):
                protein_id = line.strip().split()[0][1:]
                proteins.add(protein_id)
                gene_id = re.sub(r"\.p\d+$", "", protein_id)
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
            gene_name = re.sub(r"\.p\d+$", "", protein_name)
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
        f"Proteins with IPR: {num_ipr_prot} ({ipr_prot_percent:.2f}%)\n"
        f"Genes with IPR:    {num_ipr_gene} ({ipr_gene_percent:.2f}%)\n\n"
        "[GO Annotation]\n"
        f"Proteins with GO:  {num_go_prot} ({go_prot_percent:.2f}%)\n"
        f"Genes with GO:     {num_go_gene} ({go_gene_percent:.2f}%)\n\n"
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

# -----------------------------------------------------------------------
# KAAS Processing
# -----------------------------------------------------------------------
def run_kaas(args):
    file_path = args.query
    out_dir = args.outdir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

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
                gene_name = protein_name.split(".", 1)[0]
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
    summary_lines.append(f"KAAS annotated Protein seq:      {f2_count} ({ratio_lines:.2f}%)")
    summary_lines.append(f"Total unique genes (first col):  {total_genes}")
    summary_lines.append(f"Annotated unique genes:          {annotated_genes} ({ratio_genes:.2f}%)")

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

# -----------------------------------------------------------------------
# Diamond BLASTP to NR Processing ("diamond-nr" mode)
# -----------------------------------------------------------------------
def run_nr(args):
    proteins = set()
    genes = set()
    with open(args.protein, "r") as f:
        for line in f:
            if line.startswith(">"):
                protein_id = line.strip().split()[0][1:]
                proteins.add(protein_id)
                gene_id = re.sub(r"\.p\d+$", "", protein_id)
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
            gene_id = re.sub(r"\.p\d+$", "", qseqid)
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
    output_lines.append(f"Protein-level annotation percentage: {pct_protein_hits:.2f}%")
    output_lines.append("")
    output_lines.append(f"Total genes in FASTA: {num_genes}")
    output_lines.append(f"Genes with Diamond BLASTP hits: {num_hits_gene}")
    output_lines.append(f"Gene-level annotation percentage: {pct_gene_hits:.2f}%")
    output_lines.append("")
    output_lines.append("[Uncharacterized / Hypothetical / Low-quality / Unnamed] (Protein level)")
    output_lines.append("Counts among all proteins, plus two percentages:")
    output_lines.append("  1) % among ALL proteins")
    output_lines.append("  2) % among proteins WITH hits")
    output_lines.append("")
    output_lines.append("Uncharacterized:")
    output_lines.append(f"  Count: {unchar_count_p}")
    output_lines.append(f"  Among ALL proteins: {unchar_pct_all_p:.2f}%")
    output_lines.append(f"  Among proteins WITH hits: {unchar_pct_hits_p:.2f}%")
    output_lines.append("")
    output_lines.append("Hypothetical protein:")
    output_lines.append(f"  Count: {hypot_count_p}")
    output_lines.append(f"  Among ALL proteins: {hypot_pct_all_p:.2f}%")
    output_lines.append(f"  Among proteins WITH hits: {hypot_pct_hits_p:.2f}%")
    output_lines.append("")
    output_lines.append("Low quality protein:")
    output_lines.append(f"  Count: {lowq_count_p}")
    output_lines.append(f"  Among ALL proteins: {lowq_pct_all_p:.2f}%")
    output_lines.append(f"  Among proteins WITH hits: {lowq_pct_hits_p:.2f}%")
    output_lines.append("")
    output_lines.append("Unnamed protein product:")
    output_lines.append(f"  Count: {unnamed_count_p}")
    output_lines.append(f"  Among ALL proteins: {unnamed_pct_all_p:.2f}%")
    output_lines.append(f"  Among proteins WITH hits: {unnamed_pct_hits_p:.2f}%")
    output_lines.append("")
    output_lines.append("[Uncharacterized / Hypothetical / Low-quality / Unnamed] (Gene level)")
    output_lines.append("Counts among all genes, plus two percentages:")
    output_lines.append("  1) % among ALL genes")
    output_lines.append("  2) % among genes WITH hits")
    output_lines.append("")
    output_lines.append("Uncharacterized genes:")
    output_lines.append(f"  Count: {unchar_count_g}")
    output_lines.append(f"  Among ALL genes: {unchar_pct_all_g:.2f}%")
    output_lines.append(f"  Among genes WITH hits: {unchar_pct_hits_g:.2f}%")
    output_lines.append("")
    output_lines.append("Hypothetical genes:")
    output_lines.append(f"  Count: {hypot_count_g}")
    output_lines.append(f"  Among ALL genes: {hypot_pct_all_g:.2f}%")
    output_lines.append(f"  Among genes WITH hits: {hypot_pct_hits_g:.2f}%")
    output_lines.append("")
    output_lines.append("Low-quality genes:")
    output_lines.append(f"  Count: {lowq_count_g}")
    output_lines.append(f"  Among ALL genes: {lowq_pct_all_g:.2f}%")
    output_lines.append(f"  Among genes WITH hits: {lowq_pct_hits_g:.2f}%")
    output_lines.append("")
    output_lines.append("Unnamed-product genes:")
    output_lines.append(f"  Count: {unnamed_count_g}")
    output_lines.append(f"  Among ALL genes: {unnamed_pct_all_g:.2f}%")
    output_lines.append(f"  Among genes WITH hits: {unnamed_pct_hits_g:.2f}%")
    output_lines.append("")
    output_lines.append("[Informative Ratio Among Proteins With Hits]")
    informative_hits_p = hits_protein - (unchar_proteins | hypot_proteins | low_quality_proteins | unnamed_proteins)
    count_informative_p = len(informative_hits_p)
    ratio_informative_p = pct(count_informative_p, num_hits_protein)
    output_lines.append(f"  Informative proteins (hits): {count_informative_p}")
    output_lines.append(f"  Among proteins WITH hits: {ratio_informative_p:.2f}%")
    output_lines.append("")
    output_lines.append("[Informative Ratio Among Genes With Hits]")
    informative_hits_g = hits_gene - (unchar_genes | hypot_genes | low_quality_genes | unnamed_genes)
    count_informative_g = len(informative_hits_g)
    ratio_informative_g = pct(count_informative_g, num_hits_gene)
    output_lines.append(f"  Informative genes (hits): {count_informative_g}")
    output_lines.append(f"  Among genes WITH hits: {ratio_informative_g:.2f}%")
    output_lines.append("")
    informative_genes_file = f"{args.outprefix}.informative_gene_list.tsv"
    output_lines.append(f"Informative gene list written to: {informative_genes_file}")
    output_lines.append("")
    if species_counter:
        most_common_3 = species_counter.most_common(3)
        top_three_sum = sum(count for spec, count in most_common_3)
        top_three_ratio = (top_three_sum / sum(species_counter.values()) * 100) if sum(species_counter.values()) else 0
        output_lines.append(f"Top three species total hits: {top_three_sum} ({top_three_ratio:.2f}% of total species hits)")
    output_lines.append("Top 3 species among hits (and ratio among ALL proteins):")
    if species_counter:
        for i, (spec, count) in enumerate(species_counter.most_common(3), start=1):
            species_ratio_all = pct(count, num_proteins)
            output_lines.append(f"{i}. {spec} ({count} hits) => {species_ratio_all:.2f}% of all proteins")
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

# -----------------------------------------------------------------------
# Main Parser
# -----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Combined script for Protein Structural Annotation (PSA), InterPro, KAAS, and Diamond BLASTP (NR) processing."
    )
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    psa_parser = subparsers.add_parser("psa", help="Protein Structural Annotation mode (process GFF QC and gene length histogram).")
    psa_parser.add_argument("-g", "--gff", required=True, help="Input GFF file")
    psa_parser.add_argument("-f", "--fasta", required=True, help="Input genome FASTA file (.fna)")
    psa_parser.add_argument("-b", "--basename", required=True, help="Base name for output files")
    psa_parser.add_argument("-o", "--outdir", default=".", help="Output directory for PSA mode files (default: current directory)")
    psa_parser.add_argument("-p", "--protein", help="Protein sequence provided for BUSCO (ex. braker.aa path)")
    psa_parser.add_argument("--busco_lineage", help="BUSCO lineage (e.g. bunco lineage)")
    psa_parser.add_argument("-t", "--threads", type=int, help="Number of threads for BUSCO and OMArk")
    psa_parser.add_argument("--busco_out_path", help="Full host path for BUSCO results output (this directory will be mounted under /output in the container)")
    psa_parser.add_argument("--busco_docker_image", default="busco/busco:latest", help="Docker image for BUSCO analysis (default: busco/busco:latest)")
    psa_parser.add_argument("--busco_workdir", help="Optional: Docker working directory to mount as /output. If provided, it is used as the container's working directory.")
    psa_parser.add_argument("--omark_dir", help="Directory for OMArk run output (for -d option)")
    psa_parser.set_defaults(func=run_psa)

    ip = subparsers.add_parser("interpro", help="Process an InterProScan TSV file.")
    ip.add_argument("-p", "--protein", required=True, help="Input protein FASTA file")
    ip.add_argument("-i", "--interpro_tsv", required=True, help="InterProScan TSV file")
    ip.add_argument("-o", "--outprefix", default="interpro_output", help="Prefix for output annotation files.")
    ip.set_defaults(func=run_interpro)

    ka = subparsers.add_parser("kaas", help="Process a KAAS result file and generate a summary file. Usage: python eva_annotation.py kaas query.ko [outdir]")
    ka.add_argument("query", help="KAAS result file (e.g., query.ko)")
    ka.add_argument("outdir", nargs="?", default=".", help="Output directory (default: current directory)")
    ka.set_defaults(func=run_kaas)

    dn = subparsers.add_parser("diamond-nr", help="Evaluate Diamond BLASTP TSV output (NR database mode).")
    dn.add_argument("-p", "--protein", required=True, help="Path to the input protein FASTA file.")
    dn.add_argument("-d", "--diamond_tsv", required=True, help="Path to the Diamond BLASTP TSV output file.")
    dn.add_argument("-o", "--outprefix", default="output", help="Prefix for output files (default: 'output').")
    dn.set_defaults(func=run_nr)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
