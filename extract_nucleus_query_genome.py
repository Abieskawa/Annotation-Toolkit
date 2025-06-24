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
        "-on", "--output_nuclear", required=False,
        help="[Optional] Output FASTA file for nuclear sequences (longest chromosomes)."
    )
    parser.add_argument(
        "-os", "--output_special", required=False,
        help="[Optional]Output FASTA file for special sequences (e.g. mitochondria, chloroplast, scaffold names)."
    )
    parser.add_argument(
        "-n", "--num", type=int, required=False,
        help="[Optional] Number of longest nuclear chromosomes to extract."
    )
    parser.add_argument(
        "-s", "--special", nargs="+", required=False,
        help="[Optional]List of keywords to locate special sequences (e.g. 'mitochondria', 'chloroplast', scaffold names)."
    )
    parser.add_argument(
        "-e", "--exclude", nargs="+", required=False,
        help="[Optional]List of keywords to exclude sequences."
    )
    
    args = parser.parse_args()

    # Parse all records from the input FASTA file.
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        print("No sequences found in the input file.")
        return

    # Apply exclusion filter FIRST if specified
    if args.exclude:
        filtered_records = []
        excluded_count = 0
        for record in records:
            exclude_this = False
            for exclude_keyword in args.exclude:
                if exclude_keyword.lower() in record.description.lower():
                    exclude_this = True
                    excluded_count += 1
                    break
            if not exclude_this:
                filtered_records.append(record)
        records = filtered_records
        print(f"Excluded {excluded_count} sequences based on keywords: {args.exclude}")

    # Extract nuclear sequences: Only if both -on and -n are provided
    if args.output_nuclear and args.num:
        sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=True)
        nuclear_records = sorted_records[:args.num]
        SeqIO.write(nuclear_records, args.output_nuclear, "fasta")
        print(f"Nuclear extraction complete. {len(nuclear_records)} sequences written to {args.output_nuclear}.")
    elif args.output_nuclear or args.num:
        print("Both nuclear output file (-on) and number (-n) must be provided to perform nuclear extraction. Skipping nuclear extraction.")

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