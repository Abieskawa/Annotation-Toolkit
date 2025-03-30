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
    # Create a mutually exclusive group that requires one of --biotype or --feature to be set.
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--biotype", help="Only include features with gene_biotype equal to this value")
    group.add_argument("--feature", help="Only include features with third column (feature type) equal to this value")
    parser.add_argument("--exclude_biotype", help="Optional: Exclude features with gene_biotype equal to this value")
    parser.add_argument("--genome_fasta", help="Optional: Genome FASTA to clamp last bin to chromosome length")
    parser.add_argument("--seqids", help="Optional: File containing list of sequence IDs (one per line) to include")
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
                if current_chr is not None:
                    lengths[current_chr] = current_len
                current_chr = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        if current_chr is not None:
            lengths[current_chr] = current_len

    return lengths

def read_gff(filename, include_biotype=None, exclude_biotype=None, include_feature=None, seqids=None):
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
            if seqids and chrom not in seqids:
                continue
            # Check the feature filter: if --feature was provided, then only include lines with matching third column.
            if include_feature and fields[2] != include_feature:
                continue
            try:
                start = int(fields[3])
            except ValueError:
                continue
            attributes = fields[8]
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
            if include_biotype and biotype_val != include_biotype:
                continue
            if exclude_biotype and biotype_val == exclude_biotype:
                continue

            chrom_starts[chrom].append(start)
    
    if duplicates:
        print("Error: Multiple entries with same IDs found: " + ", ".join(sorted(duplicates)))
        exit(1)
    return chrom_starts

def count_features(chrom_starts, window_size, chr_lengths=None):
    results = {}
    for chrom, starts in chrom_starts.items():
        if not starts:
            continue
        max_coord = max(starts)
        chr_len = chr_lengths[chrom] if (chr_lengths and chrom in chr_lengths) else max_coord
        n_windows = math.ceil(chr_len / window_size)
        counts = [0] * n_windows

        for s in starts:
            idx = (s - 1) // window_size
            counts[idx] += 1

        window_list = []
        for i, count in enumerate(counts):
            w_start = i * window_size + 1
            w_end   = (i + 1) * window_size
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
    
    chr_lengths = None
    if args.genome_fasta:
        chr_lengths = read_fasta_lengths(args.genome_fasta)
        print(f"Parsed chromosome lengths from {args.genome_fasta}")

    allowed_seqids = None
    if args.seqids:
        with open(args.seqids, "r") as f:
            allowed_seqids = {line.strip() for line in f if line.strip()}

    chrom_starts = read_gff(
        args.input,
        include_biotype=args.biotype,
        exclude_biotype=args.exclude_biotype,
        include_feature=args.feature,
        seqids=allowed_seqids
    )
    window_results = count_features(chrom_starts, args.window, chr_lengths)
    write_results(window_results, args.output)
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
