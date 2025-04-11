#!/usr/bin/env Rscript
# Set a CRAN mirror so that install.packages works in non‐interactive mode.
options(repos = c(CRAN = "https://cloud.r-project.org"))

# List of required packages
required_pkgs <- c("ggplot2", "ggpolypath", "XVector", "Biostrings", "venn")

# Install any missing packages
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
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

Options:
  -i NAME FILE     Add a gene list with a name.
                   Example: -i GO go_genes.tsv
  -a FILE          Add a gene list with the default name 'all' from a FASTA file.
  -o FILE          Specify output filename (default: intersect.png).
  -h, --help       Show this help message and exit.

Example:
  ./vennplot.r -i GO go.tsv -i IPR ipr.tsv -a sequences.fasta -o venn.png
\n")
  quit(status = 0)
}

# Parse command-line arguments
input_list <- list()
i <- 1
outputfile <- "intersect.png"
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
  } else if (flag == "-o") {
    outputfile <- args[i + 1]
    i <- i + 2
  } else {
    stop("Unknown option: ", flag)
  }
}

# Helper function to trim gene names (remove everything from the first dot)
trim_gene <- function(x) sub("\\..*$", "", x)

# Initialize gene_sets list
gene_sets <- list()

# First, process the "all" set (FASTA file) if provided.
if ("all" %in% names(input_list)) {
  aa_set <- readAAStringSet(input_list[["all"]])
  # Trim names (e.g. "T_jap_gene_atp6.p" becomes "T_jap_gene_atp6")
  all_genes <- trim_gene(names(aa_set))
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
          lty = 1,                # line type for borders
          ggplot = TRUE)          # produce the plot as a ggplot object

# Debug output: inspect the computed intersections in the returned object
cat("Computed intersections:\n")
print(attr(p, "intersections"))

# Save the plot using ggsave (p is a ggplot object)
ggsave(filename = outputfile, plot = p)
cat("Plot saved to", outputfile, "\n")
