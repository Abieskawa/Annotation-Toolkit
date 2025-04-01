#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def main():
    # Set up command-line argument parsing.
    parser = argparse.ArgumentParser(
        description="Extract longest N nuclear chromosomes and special sequences (e.g. mitochondria, chloroplast, scaffold names) from a FASTA file, outputting to separate files."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FASTA file containing sequences."
    )
    parser.add_argument(
        "-on", "--output_nuclear", required=True,
        help="Output FASTA file for nuclear sequences (longest chromosomes)."
    )
    parser.add_argument(
        "-os", "--output_special", required=False,
        help="Output FASTA file for special sequences (e.g. mitochondria, chloroplast, scaffold names)."
    )
    parser.add_argument(
        "-n", "--num", type=int, required=True,
        help="Number of longest nuclear chromosomes to extract."
    )
    parser.add_argument(
        "-s", "--special", nargs="+", required=False,
        help="List of keywords to locate special sequences (e.g. 'mitochondria', 'chloroplast', scaffold names)."
    )
    args = parser.parse_args()

    # Parse all records from the input FASTA file.
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        print("No sequences found in the input file.")
        return

    # Extract nuclear sequences: sort records by sequence length and take the top N.
    sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=True)
    nuclear_records = sorted_records[:args.num]
    SeqIO.write(nuclear_records, args.output_nuclear, "fasta")
    print(f"Nuclear extraction complete. {len(nuclear_records)} sequences written to {args.output_nuclear}.")

    # Special extraction: Only run if both special keywords and output file are provided.
    if args.special and args.output_special:
        special_records = []
        for keyword in args.special:
            found = False
            for record in records:
                # Use a simple substring search (case-insensitive).
                if keyword.lower() in record.description.lower():
                    if record not in special_records:
                        special_records.append(record)
                    found = True
                    break
            if not found:
                print(f"Warning: Special sequence not found for keyword: {keyword}")
        SeqIO.write(special_records, args.output_special, "fasta")
        print(f"Special extraction complete. {len(special_records)} sequences written to {args.output_special}.")
    elif args.special or args.output_special:
        print("Both special keywords (-s) and output file (-os) must be provided to perform special extraction. Skipping special extraction.")

if __name__ == "__main__":
    main()
