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
    # Optionally strip leading and trailing 'N' characters
    if strip_n:
        records = [(h, s.strip('N')) for h, s in records]

    # Optionally remove records with a matching scaffold name.
    # Here we assume the scaffold name is the first token (without the '>').
    if remove_list:
        filtered = []
        for h, s in records:
            scaffold_name = h[1:].split()[0]
            if scaffold_name in remove_list:
                continue
            filtered.append((h, s))
        records = filtered

    if sort_length:
        records = sorted(records, key=lambda x: len(x[1]), reverse=True)

    # Rename headers to >chrN and create a mapping (original -> new).
    mapping = {}
    new_records = []
    for idx, (h, s) in enumerate(records, start=1):
        original_name = h[1:].split()[0]  # take the first word after '>'
        new_name = f"chr{idx}"
        mapping[original_name] = new_name
        new_records.append((f">{new_name}", s))
    return new_records, mapping

def write_mapping(mapping, path):
    with open(path, 'w') as f:
        f.write("Original\tNew\n")
        for orig, new in mapping.items():
            f.write(f"{orig}\t{new}\n")

def modify_gff(mapping, infile, outfile):
    """Replace the first column of a GFF file using the provided mapping.
    If a record’s first column isn’t found in the mapping, it is left unchanged."""
    with open(infile) as fin, open(outfile, 'w') as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            else:
                parts = line.strip().split("\t")
                if parts:
                    if parts[0] in mapping:
                        parts[0] = mapping[parts[0]]
                fout.write("\t".join(parts) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers to >chrN, optionally sort by length, remove specified scaffolds, strip N's, generate a mapping table, and modify a GFF file accordingly."
    )
    parser.add_argument("-i", "--infile", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--outfile", required=True, help="Output FASTA file")
    parser.add_argument("--sort-length", action="store_true", help="Sort sequences by length (longest first)")
    parser.add_argument("--remove", nargs="*", default=[], help="Scaffold names to remove (exact match of the first token in header)")
    parser.add_argument("--strip-n", action="store_true", help="Strip leading and trailing 'N' characters from sequences")
    parser.add_argument("--mapping-out", help="Output file for the mapping table (tab-separated).")
    parser.add_argument("--gff-in", help="Input GFF file to modify")
    parser.add_argument("--gff-out", help="Output modified GFF file")
    args = parser.parse_args()

    records = parse_fasta(args.infile)
    new_records, mapping = process_records(records, args.sort_length, args.remove, args.strip_n)
    write_fasta(new_records, args.outfile)
    
    if args.mapping_out:
        write_mapping(mapping, args.mapping_out)
    
    if args.gff_in and args.gff_out:
        modify_gff(mapping, args.gff_in, args.gff_out)
