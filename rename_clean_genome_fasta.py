#!/usr/bin/env python3
import argparse

def parse_fasta(path):
    records = []
    with open(path) as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    records.append((header, ''.join(seq)))
                header, seq = line, []
            else:
                seq.append(line)
        if header:
            records.append((header, ''.join(seq)))
    return records

def write_fasta(records, path):
    with open(path, 'w') as f:
        for h, seq in records:
            f.write(h + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def process_records(records, sort_length, remove_list, strip_n):
    # Optionally remove leading and trailing 'N' characters from each sequence.
    if strip_n:
        records = [(h, s.strip('N')) for h, s in records]

    # If remove_list is provided, remove all records with matching scaffold names.
    if remove_list:
        filtered = []
        for h, s in records:
            scaffold_name = h[1:].split("~", 1)[0]
            if scaffold_name in remove_list:
                continue
            filtered.append((h, s))
        records = filtered

    if sort_length:
        records = sorted(records, key=lambda x: len(x[1]), reverse=True)

    # Rename each header to simply >scaffoldN
    new_records = [(f">scaffold{idx}", s) for idx, (h, s) in enumerate(records, start=1)]
    return new_records

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers to >scaffoldN, optionally sort by length, remove all occurrences of specified scaffolds, and strip beginning and ending N's."
    )
    parser.add_argument("-i", "--infile", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--outfile", required=True, help="Output FASTA file")
    parser.add_argument("--sort-length", action="store_true", help="Sort sequences by length (longest first)")
    parser.add_argument("--remove", nargs="*", default=[], help="Scaffold names to remove all occurrences from (matching text before '~')")
    parser.add_argument("--strip-n", action="store_true", help="Remove leading and trailing 'N' characters from sequences")
    args = parser.parse_args()

    records = parse_fasta(args.infile)
    new_records = process_records(records, args.sort_length, args.remove, args.strip_n)
    write_fasta(new_records, args.outfile)
