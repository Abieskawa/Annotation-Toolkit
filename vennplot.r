#!/usr/bin/env Rscript
# Set a CRAN mirror so that install.packages works in non‐interactive mode.
options(repos = c(CRAN = "https://cloud.r-project.org"))

# List of required packages
required_pkgs <- c("ggplot2", "ggpolypath", "XVector", "Biostrings", "venn")

# Install any missing packages
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # XVector and Biostrings are Bioconductor packages
    if (pkg %in% c("XVector", "Biostrings")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(ggplot2)        # for saving the plot and ggplot2 functions
  library(ggpolypath)     # needed by venn for ggplot2-based plotting
  library(XVector)        # for handling FASTA files
  library(Biostrings)     # for reading AAStringSet from FASTA
  library(venn)           # for generating Venn diagrams
})

args <- commandArgs(trailingOnly = TRUE)

# Show help if no arguments or a help flag is used
if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  cat("
Usage:
  vennplot.r -i NAME FILE -i NAME FILE ... [-o OUTPUTFILE]

Description:
  This script generates a Venn diagram based on gene sets.
  You can input multiple named gene lists using -i flags.
  Optionally, specify a FASTA file using -a (reads sequence names from a FASTA file)
  and an output image file with -o.
  Optionally, specify a GFF3 file using -g (maps protein IDs to gene IDs and can define the 'all' gene universe).

Options:
  -i NAME FILE     Add a gene list with a name.
                   Example: -i GO go_genes.tsv
  -a FILE          Add a gene list with the default name 'all' from a protein/gene FASTA file.
                   But make sure only protein-coding gene inside fasta sequence.
  -g FILE          GFF3 file. If provided, tool-specific IDs will be mapped to genes using GFF
                   (CDS protein_id -> Parent chain -> gene). If -a is not provided, 'all' will be set
                   to all genes in the GFF3.
  -o FILE          Specify output filename (default: intersect.png).
  -h, --help       Show this help message and exit.

Example:
  ./vennplot.r -i GO go.tsv -i IPR ipr.tsv -i diamond diamond.tsv -g genome.gff -o venn.png
\n")
  quit(status = 0)
}

# -------------------- GFF mapper helpers (NEW) --------------------
parse_gff_attrs <- function(attr_str) {
  attrs <- list()
  if (is.na(attr_str) || attr_str == "" || attr_str == ".") return(attrs)
  parts <- strsplit(attr_str, ";", fixed = TRUE)[[1]]
  for (p in parts) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) == 2) {
      key <- kv[1]
      val <- kv[2]
      attrs[[key]] <- val
    }
  }
  attrs
}

# Helper function to trim gene names (remove everything from the first dot),
# and normalize common GFF/Gene prefixes so lists can intersect (e.g. gene-LOC... vs LOC...)
trim_gene <- function(x) {
  x <- sub("\\..*$", "", x)
  x <- sub("^gene-", "", x)
  x <- sub("^rna-", "", x)
  x <- sub("^cds-", "", x)
  x <- sub("^id-", "", x)
  x
}

gene_label_from_gene_attrs <- function(attrs) {
  # Priority: gene -> locus_tag -> Name -> ID (normalized)
  if (!is.null(attrs[["gene"]]) && nzchar(attrs[["gene"]])) return(trim_gene(attrs[["gene"]]))
  if (!is.null(attrs[["locus_tag"]]) && nzchar(attrs[["locus_tag"]])) return(trim_gene(attrs[["locus_tag"]]))
  if (!is.null(attrs[["Name"]]) && nzchar(attrs[["Name"]])) return(trim_gene(attrs[["Name"]]))
  if (!is.null(attrs[["ID"]]) && nzchar(attrs[["ID"]])) return(trim_gene(attrs[["ID"]]))
  NA_character_
}

# Build:
# - protein_to_gene: named character vector, names are protein IDs (qseqid) and values are gene labels (LOC...)
# - gff_all_genes: character vector of all gene labels in GFF
build_gff_mapper <- function(gff_file) {
  con <- file(gff_file, open = "r")
  on.exit(close(con), add = TRUE)

  # child_id -> parent_id (store first parent only; OK for typical NCBI Parent chains)
  parent_of <- new.env(hash = TRUE, parent = emptyenv())
  type_of   <- new.env(hash = TRUE, parent = emptyenv())
  attrs_of  <- new.env(hash = TRUE, parent = emptyenv())

  gene_ids <- character(0)
  gene_label <- new.env(hash = TRUE, parent = emptyenv())  # gene_feature_id -> label
  gene_label_rev <- new.env(hash = TRUE, parent = emptyenv()) # label -> gene_feature_id

  # Store CDS->protein_id mapping candidates so we can resolve after we have parent tables
  cds_records <- list()

  while (TRUE) {
    ln <- readLines(con, n = 50000, warn = FALSE)
    if (length(ln) == 0) break
    # Skip comments
    ln <- ln[!grepl("^#", ln)]
    if (length(ln) == 0) next

    fields <- strsplit(ln, "\t", fixed = FALSE)
    for (f in fields) {
      if (length(f) < 9) next
      ftype <- f[3]
      attr_str <- f[9]
      attrs <- parse_gff_attrs(attr_str)

      fid <- if (!is.null(attrs[["ID"]])) attrs[["ID"]] else NA_character_
      if (!is.na(fid) && nzchar(fid)) {
        type_of[[fid]] <- ftype
        attrs_of[[fid]] <- attr_str
      }

      # parent table
      if (!is.na(fid) && nzchar(fid) && !is.null(attrs[["Parent"]]) && nzchar(attrs[["Parent"]])) {
        # Take first parent if multiple
        par <- strsplit(attrs[["Parent"]], ",", fixed = TRUE)[[1]][1]
        parent_of[[fid]] <- par
      }

      # gene feature label table
      if (ftype == "gene" && !is.na(fid) && nzchar(fid)) {
        lab <- gene_label_from_gene_attrs(attrs)
        if (!is.na(lab) && nzchar(lab)) {
          gene_ids <- c(gene_ids, fid)
          gene_label[[fid]] <- lab
          gene_label_rev[[lab]] <- fid
        }
      }

      # capture CDS record
      if (ftype == "CDS") {
        pid <- if (!is.null(attrs[["protein_id"]]) && nzchar(attrs[["protein_id"]])) attrs[["protein_id"]] else NA_character_
        # special NCBI cases may not have protein_id; use gene/locus_tag/Name as proxy "protein ID"
        proxy <- NA_character_
        if (is.na(pid) || !nzchar(pid)) {
          if (!is.null(attrs[["gene"]]) && nzchar(attrs[["gene"]])) proxy <- trim_gene(attrs[["gene"]])
          else if (!is.null(attrs[["locus_tag"]]) && nzchar(attrs[["locus_tag"]])) proxy <- trim_gene(attrs[["locus_tag"]])
          else if (!is.null(attrs[["Name"]]) && nzchar(attrs[["Name"]])) proxy <- trim_gene(attrs[["Name"]])
        }
        par <- if (!is.null(attrs[["Parent"]]) && nzchar(attrs[["Parent"]])) strsplit(attrs[["Parent"]], ",", fixed = TRUE)[[1]][1] else NA_character_
        cds_records[[length(cds_records) + 1]] <- list(protein_id = pid, proxy_id = proxy, parent = par, cds_gene = if (!is.null(attrs[["gene"]])) trim_gene(attrs[["gene"]]) else NA_character_)
      }
    }
  }

  # Walk up parent chain to find a gene feature ID
  resolve_to_gene_id <- function(start_id) {
    cur <- start_id
    steps <- 0
    while (!is.na(cur) && nzchar(cur) && steps < 20) {
      if (!is.null(type_of[[cur]]) && identical(type_of[[cur]], "gene")) return(cur)
      if (!is.null(parent_of[[cur]])) {
        cur <- parent_of[[cur]]
      } else {
        return(NA_character_)
      }
      steps <- steps + 1
    }
    NA_character_
  }

  protein_to_gene <- character(0)

  # Resolve CDS-derived mappings
  for (rec in cds_records) {
    par <- rec$parent
    gid <- resolve_to_gene_id(par)

    # If parent chain didn't reach gene, try using CDS gene label to find gene feature
    if (is.na(gid) || !nzchar(gid)) {
      if (!is.na(rec$cds_gene) && nzchar(rec$cds_gene) && !is.null(gene_label_rev[[rec$cds_gene]])) {
        gid <- gene_label_rev[[rec$cds_gene]]
      }
    }

    if (!is.na(gid) && nzchar(gid) && !is.null(gene_label[[gid]])) {
      glab <- gene_label[[gid]]
      if (!is.na(rec$protein_id) && nzchar(rec$protein_id)) {
        protein_to_gene[rec$protein_id] <- glab
      } else if (!is.na(rec$proxy_id) && nzchar(rec$proxy_id)) {
        # Proxy mapping: LOC... (as "protein") -> gene label
        protein_to_gene[rec$proxy_id] <- glab
      }
    }
  }

  list(
    protein_to_gene = protein_to_gene,
    gff_all_genes = unique(unname(as.list(gene_label)))
  )
}
# ------------------ end GFF mapper helpers ------------------

# Parse command-line arguments
input_list <- list()
i <- 1
outputfile <- "intersect.png"
gff_file <- NULL
while (i <= length(args)) {
  flag <- args[i]
  if (flag == "-i") {
    name <- args[i + 1]
    file <- args[i + 2]
    input_list[[name]] <- file
    i <- i + 3
  } else if (flag == "-a") {
    # Use default name "all" for FASTA input
    name <- "all"
    file <- args[i + 1]
    input_list[[name]] <- file
    i <- i + 2
  } else if (flag == "-g") {
    gff_file <- args[i + 1]
    i <- i + 2
  } else if (flag == "-o") {
    outputfile <- args[i + 1]
    i <- i + 2
  } else {
    stop("Unknown option: ", flag)
  }
}

# Build GFF mapper if provided
gff_mapper <- NULL
if (!is.null(gff_file)) {
  gff_mapper <- build_gff_mapper(gff_file)
}

# Initialize gene_sets list
gene_sets <- list()

# First, process the "all" set (FASTA file) if provided.
if ("all" %in% names(input_list)) {
  aa_set <- readAAStringSet(input_list[["all"]])
  # Trim names (e.g. "T_jap_gene_atp6.p" becomes "T_jap_gene_atp6")
  all_genes <- trim_gene(names(aa_set))
  # If GFF mapper is provided, map FASTA IDs to genes where possible
  if (!is.null(gff_mapper)) {
    mapped <- gff_mapper$protein_to_gene[all_genes]
    use <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, all_genes)
    all_genes <- use
  }
  gene_sets[["all"]] <- all_genes
} else if (!is.null(gff_mapper)) {
  # If FASTA isn't provided, but GFF is, define "all" as all genes in the GFF
  all_genes <- gff_mapper$gff_all_genes
  gene_sets[["all"]] <- all_genes
} else {
  all_genes <- NULL
}

# Set this flag to TRUE if you want to filter tool-specific gene lists by the "all" set.
# Set to FALSE if the tool-specific names (e.g., "T_jap_g1") are not expected to match the FASTA names.
filter_by_all <- FALSE

# Process each tool-specific gene list.
for (tool in setdiff(names(input_list), "all")) {
  # Read the file and take only its first column (to avoid extra tokens like headers)
  df <- read.delim(input_list[[tool]], header = FALSE, stringsAsFactors = FALSE)
  genes <- df[[1]]
  # Apply the same trimming to the gene list
  trimmed_genes <- trim_gene(genes)

  # If GFF mapper exists, map protein IDs to gene labels (LOC...)
  if (!is.null(gff_mapper)) {
    mapped <- gff_mapper$protein_to_gene[trimmed_genes]
    trimmed_genes <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, trimmed_genes)
  }

  # Optionally filter out entries not present in the "all" set
  if (!is.null(all_genes) && filter_by_all) {
    trimmed_genes <- trimmed_genes[trimmed_genes %in% all_genes]
  }
  gene_sets[[tool]] <- trimmed_genes
}

# Debug output: show the sizes of each gene set
cat("Gene set sizes:\n")
print(sapply(gene_sets, length))

# ------- Plotting using the venn package -------
# Now, generate the Venn diagram.
# Use ilabels = "counts" so that it displays intersection counts (instead of default concatenated labels).
p <- venn(gene_sets,
          ilabels = "counts",     # show computed counts for each intersection
          zcolor = c("#b90489", "#04c6f7f1", "#00a651", "#ffcc00", "#ff8400"),
          col = c("red", "black", "black", "black", "black"),
          lty = 1,
          ilcs = 1, sncs = 1,                # line type for borders
          ggplot = TRUE)          # produce the plot as a ggplot object

# Debug output: inspect the computed intersections in the returned object
cat("Computed intersections:\n")
print(attr(p, "intersections"))

# Save the plot using ggsave (p is a ggplot object)
ggsave(filename = outputfile, plot = p)
cat("Plot saved to", outputfile, "\n")
