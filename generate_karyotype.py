#!/usr/bin/env python3
import argparse

def parse_fasta(fasta_file):
    """Reads a FASTA file and returns a dictionary {seq_name: sequence}."""
    sequences = {}
    with open(fasta_file) as f:
        seq_name, seq_lines = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_name is not None:
                    sequences[seq_name] = ''.join(seq_lines)
                # Use the first word after '>' as the sequence name
                seq_name = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_name is not None:
            sequences[seq_name] = ''.join(seq_lines)
    return sequences

def generate_karyotype(sequences, top_n, species_prefix):
    """
    Returns a list of Circos karyotype lines in the format:
        chr - ID LABEL 0 END COLOR
    If a header is purely numeric, or starts with "chr" or "scaffold" followed by digits,
    it extracts the numeric part. For headers starting with "scaffold":
      - The ID remains the original header (e.g., "scaffold23")
      - The LABEL and COLOR become "chr" plus the numeric part (e.g., "chr23").
    For numeric headers or those starting with "chr", the ID is set to species_prefix + number.
    Otherwise, the header is used for both ID and LABEL.
    """
    seq_lengths = [(name, len(seq)) for name, seq in sequences.items()]
    seq_lengths.sort(key=lambda x: x[1], reverse=True)
    top_chromosomes = seq_lengths[:top_n]
    
    lines = []
    for name, length in top_chromosomes:
        num_part = None
        is_scaffold = False
        if name.isdigit():
            num_part = name
        elif name.lower().startswith("chr") and name[3:].isdigit():
            num_part = name[3:]
        elif name.lower().startswith("scaffold") and name[8:].isdigit():
            num_part = name[8:]
            is_scaffold = True
            
        if num_part is not None:
            if is_scaffold:
                chrom_id = name  # leave the original scaffold name as ID
            else:
                chrom_id = f"{species_prefix}{num_part}"
            label = f"chr{num_part}"
            color = label  # set color same as label
        else:
            chrom_id = name
            label = name
            color = label
        lines.append(f"chr - {chrom_id} {label} 0 {length} {color}")
    return lines

def main():
    parser = argparse.ArgumentParser(
        description="Generate a Circos karyotype file from a FASTA file (with a color column)."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("-n", "--num", type=int, default=24, help="Number of sequences to include")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    parser.add_argument("--prefix", default="hs",
                        help="Species prefix for numeric chromosome names (default: hs)")
    args = parser.parse_args()
    
    sequences = parse_fasta(args.fasta)
    karyotype_lines = generate_karyotype(sequences, args.num, args.prefix)
    
    if args.output:
        with open(args.output, "w") as out_f:
            for line in karyotype_lines:
                out_f.write(line + "\n")
    else:
        for line in karyotype_lines:
            print(line)

if __name__ == "__main__":
    main()
