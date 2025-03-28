#!/usr/bin/env python3
import argparse

def load_targets(target_file):
    """Read a TSV file with one scaffold per line."""
    with open(target_file) as f:
        return {line.strip() for line in f if line.strip()}

def parse_fasta_lengths(fasta_file):
    """
    Parse a FASTA file and return a dictionary with scaffold names as keys
    and their lengths as values.
    """
    lengths = {}
    current_scaffold = None
    current_length = 0
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_scaffold:
                    lengths[current_scaffold] = current_length
                current_scaffold = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)
        if current_scaffold:
            lengths[current_scaffold] = current_length
    return lengths

def main():
    parser = argparse.ArgumentParser(
        description="Count unique repeat IDs in windows for specified scaffolds and repeat pattern."
    )
    parser.add_argument("input", help="Input file (RepeatMasker-like output)")
    parser.add_argument("output", help="Output file for Circos")
    parser.add_argument("--scaffold", help="Process a single scaffold (e.g. scaffold1)")
    parser.add_argument("--targets", help="TSV file with one target scaffold per line")
    parser.add_argument("--window", type=int, default=100000, help="Window size (default: 100000)")
    parser.add_argument("--pattern", default="DNA", help="Repeat pattern to filter (default: DNA)")
    parser.add_argument("--fasta", help="FASTA file to obtain scaffold lengths")
    args = parser.parse_args()

    # Build the target set: either one scaffold or a list from a file.
    targets = set()
    if args.targets:
        targets = load_targets(args.targets)
    elif args.scaffold:
        targets.add(args.scaffold)
    else:
        parser.error("Provide either --scaffold or --targets.")

    # Optionally, parse FASTA lengths.
    fasta_lengths = {}
    if args.fasta:
        fasta_lengths = parse_fasta_lengths(args.fasta)

    # windows: key = (scaffold, window_index); value = set of unique repeat IDs.
    windows = {}

    with open(args.input) as infile:
        for line in infile:
            line = line.strip()
            # Skip headers and blank lines.
            if not line or line.startswith("SW") or line.startswith("score"):
                continue
            fields = line.split()
            if len(fields) < 15:
                continue

            query = fields[4]
            if query not in targets:
                continue

            # Filter by repeat class (field index 10).
            repeat_class = fields[10]
            if not repeat_class.startswith(args.pattern):
                continue

            try:
                begin = int(fields[5])
            except ValueError:
                continue

            # Mimic BuildSummary.pl: do not strip trailing "*" from ID.
            rid = fields[14]

            win_index = (begin - 1) // args.window
            key = (query, win_index)
            windows.setdefault(key, set()).add(rid)

    # Write output: adjust the last window if scaffold length is known.
    with open(args.output, "w") as outfile:
        for (scaf, win_index) in sorted(windows.keys()):
            win_start = win_index * args.window + 1
            win_end = (win_index + 1) * args.window
            if args.fasta and scaf in fasta_lengths:
                # Do not allow window end to exceed the scaffold's length.
                if win_end > fasta_lengths[scaf]:
                    win_end = fasta_lengths[scaf]
                # Skip windows starting beyond the scaffold length.
                if win_start > fasta_lengths[scaf]:
                    continue
            count = len(windows[(scaf, win_index)])
            outfile.write(f"{scaf}\t{win_start}\t{win_end}\t{count}\n")

if __name__ == "__main__":
    main()
