#!/usr/bin/env python3
import os
import sys
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate masking statistics for a soft-masked genome."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input soft-masked genome file."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="04_trf_mask",
        help="Output directory where the results will be saved (default: 04_trf_mask)."
    )
    parser.add_argument(
        "-f",
        "--outfile",
        default="TRF_stats.tsv",
        help="Name of the output TSV file (default: TRF_stats.tsv)."
    )
    return parser.parse_args()

def read_and_concatenate(input_file):
    """
    Reads the input file, skips lines starting with '>', strips newlines,
    and concatenates the sequence into a single string.
    """
    sequence_parts = []
    try:
        with open(input_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                sequence_parts.append(line.strip())
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)
    return "".join(sequence_parts)

def calculate_stats(sequence):
    """
    Calculates the following statistics:
      - Total length of the sequence.
      - Count of masked bases (lowercase letters).
      - Count of 'N' and 'X' characters.
      - Effective length (total length minus NX bases).
      - Masking ratio (masked bases / effective length).
    """
    total = len(sequence)
    masked = sum(1 for c in sequence if c.islower())
    # Count uppercase 'N' and 'X'
    nx = sequence.count("N") + sequence.count("X")
    effective = total - nx
    ratio = masked / effective if effective > 0 else 0
    return total, masked, nx, effective, ratio

def write_stats_to_file(output_file, stats):
    header = "Total_Length\tMasked_Bases\tNX_Bases\tEffective_Length\tMasking_Ratio"
    total, masked, nx, effective, ratio = stats
    data_line = f"{total}\t{masked}\t{nx}\t{effective}\t{ratio:.4f}"
    try:
        with open(output_file, "w") as f:
            f.write(header + "\n")
            f.write(data_line + "\n")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

def print_formatted_table(file_path):
    """
    Reads the TSV file and prints a formatted table similar to the
    output of the 'column -t' command.
    """
    try:
        with open(file_path, "r") as f:
            lines = [line.rstrip().split("\t") for line in f]
    except Exception as e:
        print(f"Error reading output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Determine the maximum width for each column
    num_cols = len(lines[0])
    col_widths = [0] * num_cols
    for row in lines:
        for i, item in enumerate(row):
            col_widths[i] = max(col_widths[i], len(item))
    
    # Print each row with padded columns
    for row in lines:
        formatted_row = "  ".join(item.ljust(col_widths[i]) for i, item in enumerate(row))
        print(formatted_row)

def main():
    args = parse_arguments()
    input_file = args.input
    output_dir = args.outdir
    output_file_name = args.outfile

    # Check that the input file exists
    if not os.path.isfile(input_file):
        print(f"Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    # Create the output directory if it does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Read and process the input file
    sequence = read_and_concatenate(input_file)
    stats = calculate_stats(sequence)

    # Prepare the full path to the output file
    output_file = os.path.join(output_dir, output_file_name)

    # Write statistics to the output TSV file
    write_stats_to_file(output_file, stats)

    # Display the formatted results
    print_formatted_table(output_file)

if __name__ == "__main__":
    main()
