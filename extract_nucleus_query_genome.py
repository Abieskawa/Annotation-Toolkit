#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO

def main():
    # Set up command-line argument parsing.
    parser = argparse.ArgumentParser(
        description="Extract longest N chromosomes and special sequences (e.g. mitochondria, chloroplast) from a FASTA file."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FASTA file containing sequences."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output FASTA file to write the extracted sequences."
    )
    parser.add_argument(
        "-n", "--num", type=int, required=True,
        help="Number of longest chromosomes to extract."
    )
    parser.add_argument(
        "-s", "--special", nargs="+",
        help="List of keywords to locate special sequences (e.g. 'mitochondria' 'chloroplast')."
    )
    args = parser.parse_args()

    # Parse all records from the input FASTA file.
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        print("No sequences found in the input file.")
        return

    # Sort records by sequence length in descending order and extract the top N.
    sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=True)
    longest_chromosomes = sorted_records[:args.num]

    # Search for each special sequence by keyword in the record description using an exact word match.
    special_records = []
    if args.special:
        for keyword in args.special:
            # Build a regex pattern that matches the exact keyword as a whole word (case-insensitive)
            pattern = r'\b' + re.escape(keyword.lower()) + r'\b'
            found = False
            for record in records:
                if re.search(pattern, record.description.lower()):
                    if record not in special_records:
                        special_records.append(record)
                    found = True
                    break
            if not found:
                print(f"Warning: Special sequence not found for keyword: {keyword}")

    # Combine the longest chromosomes and special sequences.
    extracted_records = longest_chromosomes + special_records

    # Write the extracted sequences to the output FASTA file.
    SeqIO.write(extracted_records, args.output, "fasta")
    print(f"Extraction complete. {len(extracted_records)} sequences written to {args.output}.")

if __name__ == "__main__":
    main()
