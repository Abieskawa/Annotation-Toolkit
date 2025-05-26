#!/usr/bin/env python3
import argparse
import sys

def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return the header and sequence.
    This script is designed for a mitochondria FASTA input.
    If the file contains multiple sequences, only the first (mitochondria) is processed.
    """
    header = None
    sequence = ""
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line
                else:
                    sys.stderr.write("Warning: Multiple sequences detected in FASTA file. Only the first (mitochondria) sequence will be used.\n")
                    break
            else:
                if header is not None:
                    sequence += line
    return header, sequence

def write_fasta(output_file, header, sequence, line_length=100):
    """Write the sequence in FASTA format with the given header."""
    with open(output_file, "w") as f:
        f.write(header + "\n")
        for i in range(0, len(sequence), line_length):
            f.write(sequence[i:i+line_length] + "\n")

def merge_intervals(intervals):
    """Merge overlapping intervals from a list of (start, end) tuples."""
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev = merged[-1]
        # Merge if current interval overlaps or is adjacent.
        if current[0] <= prev[1] + 1:
            merged[-1] = (prev[0], max(prev[1], current[1]))
        else:
            merged.append(current)
    return merged

def fix_genes_in_gff(gff_file, L):
    """
    Read gene (and ncRNA_gene) features from the GFF file and, for those on the negative strand,
    fix problematic annotations. A negative strand gene is considered problematic if:
      1. Its length is greater than half of the sequence length, or
      2. It completely encloses another gene on the same strand, or
      3. It touches the first or last base pair of the sequence.
    For problematic genes, we fix them by shifting their coordinates so that:
         new_start = original end, and new_end = L + 1.
    In addition, we record an extra forbidden interval covering from position 1 up to the original start.
    This extra interval ensures that the sequence before the problematic gene is forbidden.
    The function prints the problematic gene's line and its original start.
    Returns a tuple (fixed_genes, extra_intervals, problematic_flag) where:
       - fixed_genes is a list of tuples (start, end) for all genes (fixed if needed),
       - extra_intervals is a list of extra forbidden intervals,
       - problematic_flag is True if any gene was fixed.
    """
    genes = []
    problematic_flag = False
    extra_intervals = []
    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2].strip()
            # Only look at "gene" entries exactly
            if feature_type != "gene":
                continue
            
            # Handle the start position - check for boundary indicators first
            start_str = parts[3].strip()
            end_str = parts[4].strip()
            
            # Handle special boundary notations
            touches_start_boundary = False
            touches_end_boundary = False
            
            try:
                if start_str.startswith('<'):
                    # Gene extends before the start of the sequence
                    s = int(start_str[1:])  # Remove '<' and parse the number
                    touches_start_boundary = True
                else:
                    s = int(start_str)
                
                if end_str.startswith('>'):
                    # Gene extends beyond the end of the sequence
                    e = int(end_str[1:])  # Remove '>' and parse the number
                    touches_end_boundary = True
                else:
                    e = int(end_str)
                    
            except ValueError:
                # Skip genes with unparseable coordinates
                continue
            
            strand = parts[6].strip()
            genes.append({
                "start": s,
                "end": e,
                "strand": strand,
                "orig_start": s,
                "orig_end": e,
                "line": line,
                "touches_start_boundary": touches_start_boundary,
                "touches_end_boundary": touches_end_boundary
            })
    
    # Process genes for problematic conditions
    fixed_genes = []
    for g in genes:
        # Condition 1: gene length greater than half of the sequence length.
        cond_length = ((g["end"] - g["start"] + 1) > (L / 2))
        
        # Condition 2: gene completely encloses another gene on the same strand.
        same_strand_genes = [other for other in genes if other["strand"] == g["strand"] and other is not g]
        cond_overlap = any((other["start"] > g["start"] and other["end"] < g["end"])
                           for other in same_strand_genes)
        
        # Condition 3: gene touches the first or last base pair of the sequence.
        cond_boundary = (g["start"] == 1 or g["end"] == L or 
                        g["touches_start_boundary"] or g["touches_end_boundary"])
        
        if cond_length or cond_overlap or cond_boundary:
            problematic_flag = True
            # Output the problematic gene information.
            reason = []
            if cond_length:
                reason.append("length greater than half of the sequence")
                fixed_genes.append((1,g["start"]))
                fixed_genes.append((g["end"], L))

            if cond_overlap:
                reason.append("completely encloses another gene")
                fixed_genes.append((1,g["start"]))
                fixed_genes.append((g["end"], L))
                
            if cond_boundary:
                if g["touches_start_boundary"]:
                    reason.append("extends before sequence start")
                    g["start"] = 1
                if g["touches_end_boundary"]:
                    reason.append("extends beyond sequence end")
                    g["end"] = L
                if g["start"] == 1:
                    reason.append("starts at position 1")
                if g["end"] == L:
                    reason.append("ends at last position")
                fixed_genes.append((g["start"], g["end"]))
            
            sys.stderr.write(f"Problematic gene detected: {g['line']}\n")
            sys.stderr.write(f"Original start: {g['orig_start']}, Original end: {g['orig_end']}\n")
            sys.stderr.write(f"Reasons: {', '.join(reason)}\n")
        else:
            fixed_genes.append((g["start"], g["end"]))
    

    return fixed_genes, extra_intervals, problematic_flag

def parse_gff(gff_file, margin, L):
    """
    Parse the GFF file and return merged forbidden intervals.
    For each gene (types "gene" or "ncRNA_gene"), using fixed coordinates if necessary,
    a forbidden interval is defined as:
         [max(1, start - margin), min(L, end + margin)]
    Also, extra forbidden intervals from problematic genes are added.
    Returns a tuple (merged_intervals, problematic_flag).
    """
    fixed_genes, extra_intervals, problematic_flag = fix_genes_in_gff(gff_file, L)
    intervals = []
    for (s, e) in fixed_genes:
        # Only create forbidden intervals for genes that are still within [1, L]
        # Fixed genes that extend beyond L should not create forbidden intervals in the original range
        if s <= L:  # Gene starts within the original sequence
            start_forbid = max(1, s - margin)
            end_forbid = min(L, e + margin)
            intervals.append((start_forbid, end_forbid))
    intervals.extend(extra_intervals)
    merged = merge_intervals(intervals)
    return merged, problematic_flag

def find_allowed_intervals(merged_intervals, L):
    """
    Compute allowed intervals in [1, L] that are not covered by any forbidden interval.
    """
    allowed = []
    current = 1
    for interval in merged_intervals:
        if current < interval[0]:
            allowed.append((current, interval[0] - 1))
        current = interval[1] + 1
    if current <= L:
        allowed.append((current, L))
    return allowed

def choose_split_point(allowed_intervals):
    """
    Choose a splitting point from allowed intervals.
    This implementation selects the midpoint of the largest allowed interval.
    """
    if not allowed_intervals:
        return None, None
    max_interval = max(allowed_intervals, key=lambda x: x[1] - x[0])
    split_point = (max_interval[0] + max_interval[1]) // 2
    return split_point, max_interval

def find_fallback_split_point(genes, L):
    """
    Find a fallback split point when no allowed intervals are found.
    Looks for the largest gap between adjacent genes.
    """
    # Sort genes by start position
    sorted_genes = sorted([(g["start"], g["end"]) for g in genes], key=lambda x: x[0])
    
    if not sorted_genes:
        return L // 2  # default to middle if no genes
    
    # Find gaps between adjacent genes
    gaps = []
    for i in range(len(sorted_genes)-1):
        current_end = sorted_genes[i][1]
        next_start = sorted_genes[i+1][0]
        if next_start > current_end:  # There's a gap
            gap_size = next_start - current_end - 1
            mid_point = current_end + gap_size // 2
            gaps.append((gap_size, mid_point))
    
    # Check gap between last gene and first gene (wrapping around)
    if sorted_genes[0][0] > 1:  # Gap at the beginning
        gap_size = sorted_genes[0][0] - 1
        mid_point = 1 + gap_size // 2
        gaps.append((gap_size, mid_point))
    
    if sorted_genes[-1][1] < L:  # Gap at the end
        gap_size = L - sorted_genes[-1][1]
        mid_point = sorted_genes[-1][1] + gap_size // 2
        gaps.append((gap_size, mid_point))
    
    if not gaps:
        # No gaps found, pick middle of sequence
        return L // 2
    
    # Find the largest gap
    largest_gap = max(gaps, key=lambda x: x[0])
    return largest_gap[1]  # Return the midpoint of the largest gap

def main():
    parser = argparse.ArgumentParser(
        description="Adjust circular mitochondria genome to avoid splitting in range of problematic genes.\n"
                    "NOTE: The input FASTA must contain only the mitochondria genome, and GFF is the rsult from mitos."
    )
    parser.add_argument("-f", "--fasta", required=True,
                        help="Input mitochondria FASTA file (should contain a single mitochondria genome)")
    parser.add_argument("-g", "--gff", required=True,
                        help="Input GFF file from first prediction (use GFF directly from Mitos)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output adjusted mitochondria FASTA file")
    parser.add_argument("-i", "--info", required=True,
                        help="Output split information text file")
    parser.add_argument("--margin", type=int, default=0,
                        help="Minimum distance (bp) from gene boundaries (default: 5)")
    args = parser.parse_args()

    # Step 1: Read the mitochondria FASTA file.
    header, sequence = parse_fasta(args.fasta)
    if header is None or not sequence:
        sys.exit("Error: No valid FASTA record found. Please check your input file.")
    L = len(sequence)
    print(f"Original mitochondria length: {L} bp")

    # Step 2: Double the sequence to emulate circularity.
    doubled_seq = sequence + sequence

    # Step 3: Parse the GFF file (with gene fixes) and determine forbidden intervals.
    merged_intervals, problematic_flag = parse_gff(args.gff, args.margin, L)
    if not problematic_flag:
        print("No adjustment necessary: no problematic genes detected. Exiting without changes.")
        sys.exit(0)
    
    # Get genes for fallback split point calculation
    genes = []
    with open(args.gff) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2].strip()
            if feature_type != "gene":
                continue
            try:
                start_str = parts[3].strip()
                if start_str.startswith('<'):
                    s = int(start_str[1:])
                else:
                    s = int(start_str)
                e = int(parts[4])
            except ValueError:
                continue
            strand = parts[6].strip()
            genes.append({
                "start": s,
                "end": e,
                "strand": strand
            })
    
    allowed_intervals = find_allowed_intervals(merged_intervals, L)
    
    # Step 4: Choose a safe split point
    if allowed_intervals:
        print(f"Found allowed intervals: {allowed_intervals}")
        split_point, chosen_interval = choose_split_point(allowed_intervals)
        print(f"Chosen split point: {split_point} (from allowed interval {chosen_interval})")
    else:
        print("No allowed intervals found. Using fallback method to find split point.")
        split_point = find_fallback_split_point(genes, L)
        print(f"Fallback split point: {split_point} (chosen from largest gap between genes)")
    
    if split_point is None:
        sys.exit("Error: Failed to determine a split point.")
    
    # Step 5: Use the chosen split point (from the first sequence) on the doubled sequence.
    adjusted_seq = doubled_seq[split_point - 1 : split_point - 1 + L]
    print(f"Adjusted sequence length: {len(adjusted_seq)} bp")
    
    new_header = header.split()[0]
    
    # Write the adjusted mitochondria genome as a FASTA file.
    write_fasta(args.output, new_header, adjusted_seq)
    
    # Write the split information into a text file.
    with open(args.info, "w") as f:
        f.write("Mitochondria Adjustment Information\n")
        f.write(f"Original length: {L} bp\n")
        f.write(f"Chosen split point (from first sequence): {split_point}\n")
        if allowed_intervals:
            f.write(f"Allowed intervals: {allowed_intervals}\n")
        else:
            f.write("No allowed intervals found. Used fallback method to find split point.\n")
        f.write("Problem detected due to genes touching sequence boundaries.\n")
        f.write("Fixed problematic gene coordinates where necessary and added extra forbidden intervals.\n")
        f.write("Adjusted genome is extracted from the doubled sequence starting at the chosen split point.\n")
        f.write("Reminder: Use the GFF file directly from Mitos for accurate gene annotation.\n")

if __name__ == '__main__':
    main()