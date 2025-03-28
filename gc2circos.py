#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def compute_gc(seq):
    """Calculate the GC percentage for a given sequence."""
    return (seq.count('G') + seq.count('C')) / len(seq) * 100 if seq else 0

def load_targets(target_file):
    """Load target scaffold names from a TSV file (one per line)."""
    targets = set()
    with open(target_file) as f:
        for line in f:
            line = line.strip()
            if line:
                targets.add(line)
    return targets

def main():
    parser = argparse.ArgumentParser(description="Calculate GC ratio in windows for a given FASTA file.")
    parser.add_argument('-i', '--input', required=True, help="Input FASTA file.")
    parser.add_argument('-o', '--output', required=True, help="Output text file.")
    parser.add_argument('-w', '--window', required=True, type=int, help="Window size (integer).")
    parser.add_argument('-t', '--target', help="Optional TSV file listing target scaffolds (one per line).")
    
    args = parser.parse_args()
    
    target_scaffolds = None
    if args.target:
        target_scaffolds = load_targets(args.target)
    
    with open(args.output, 'w') as out:
        for record in SeqIO.parse(args.input, "fasta"):
            if target_scaffolds is not None and record.id not in target_scaffolds:
                continue
            seq = str(record.seq).upper()
            for i in range(0, len(seq), args.window):
                start = i
                end = min(i + args.window, len(seq))
                win_seq = seq[start:end]
                gc_percent = compute_gc(win_seq)
                # Output: scaffold, 1-based start, end, GC percentage with two decimals
                out.write(f"{record.id}\t{start+1}\t{end}\t{gc_percent:.2f}\n")

if __name__ == '__main__':
    main()
