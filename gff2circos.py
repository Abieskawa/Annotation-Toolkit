#!/usr/bin/env python3
import argparse
import math
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description="Count GFF features per window (by start) for Circos; requires unique IDs."
    )
    parser.add_argument("--input", required=True, help="Input GFF file")
    parser.add_argument("--output", required=True, help="Output file (Circos format)")
    parser.add_argument("--window", type=int, default=100000, help="Window size in bp (default: 100000)")
    parser.add_argument("--biotype", help="Optional: Only include features with gene_biotype equal to this value")
    parser.add_argument("--exclude_biotype", help="Optional: Exclude features with gene_biotype equal to this value")
    parser.add_argument("--genome_fasta", help="Optional: Genome FASTA to clamp last bin to chromosome length")
    return parser.parse_args()

def read_fasta_lengths(fasta_file):
    """
    Parse a simple multi-FASTA, returning a dict {chr_name: length}.
    The 'chr_name' is the text after '>' up to the first whitespace.
    Lines are concatenated for each chromosome until the next '>'.
    """
    lengths = {}
    current_chr = None
    current_len = 0

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # If we were tracking a previous chromosome, store its length
                if current_chr is not None:
                    lengths[current_chr] = current_len
                # Start a new chromosome
                # Typically the ID is the first token after '>'
                # e.g. ">scaffold1 some description..."
                current_chr = line[1:].split()[0]
                current_len = 0
            else:
                # Sequence line - add its length
                current_len += len(line)
        # Don't forget the last chromosome
        if current_chr is not None:
            lengths[current_chr] = current_len

    return lengths

def read_gff(filename, include_biotype=None, exclude_biotype=None):
    chrom_starts = defaultdict(list)
    seen_ids = set()
    duplicates = set()
    
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            try:
                start = int(fields[3])
            except ValueError:
                continue
            attributes = fields[8]
            # Build a dictionary for the attributes
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key.strip()] = value.strip()
            if "ID" not in attr_dict:
                continue
            feature_id = attr_dict["ID"]
            if feature_id in seen_ids:
                duplicates.add(feature_id)
            else:
                seen_ids.add(feature_id)
            
            biotype_val = attr_dict.get("gene_biotype", None)
            # If a biotype filter is provided, only include matching features.
            if include_biotype and biotype_val != include_biotype:
                continue
            # If an exclusion filter is provided, skip matching features.
            if exclude_biotype and biotype_val == exclude_biotype:
                continue

            chrom_starts[chrom].append(start)
    
    if duplicates:
        print("Error: Multiple entries with same IDs found: " + ", ".join(sorted(duplicates)))
        exit(1)
    return chrom_starts

def count_features(chrom_starts, window_size, chr_lengths=None):
    """
    If chr_lengths is provided, the last bin for each chromosome is clamped
    to the chromosome's length.
    """
    results = {}
    for chrom, starts in chrom_starts.items():
        if not starts:
            continue
        max_coord = max(starts)
        # If no reference length is available, use the max coordinate to define windows
        # else clamp to the chromosome length if known
        chr_len = chr_lengths[chrom] if (chr_lengths and chrom in chr_lengths) else max_coord
        
        # e.g. if chr_len=50000 and window=10000 => 5 windows
        n_windows = math.ceil(chr_len / window_size)
        counts = [0] * n_windows

        # Count how many starts land in each bin
        for s in starts:
            idx = (s - 1) // window_size
            counts[idx] += 1

        # Build the bin list
        window_list = []
        for i, count in enumerate(counts):
            w_start = i * window_size + 1
            w_end   = (i + 1) * window_size
            # If we know the chromosome length, clamp the final bin
            if w_end > chr_len:
                w_end = chr_len
            window_list.append((w_start, w_end, count))

        results[chrom] = window_list
    return results

def write_results(results, output_file):
    with open(output_file, "w") as out:
        for chrom in sorted(results):
            for start, end, count in results[chrom]:
                out.write(f"{chrom}\t{start}\t{end}\t{count}\n")

def main():
    args = parse_args()
    # If user provided a FASTA, parse to get chromosome lengths
    chr_lengths = None
    if args.genome_fasta:
        chr_lengths = read_fasta_lengths(args.genome_fasta)
        print(f"Parsed chromosome lengths from {args.genome_fasta}")

    chrom_starts = read_gff(args.input, args.biotype, args.exclude_biotype)
    window_results = count_features(chrom_starts, args.window, chr_lengths)
    write_results(window_results, args.output)
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
