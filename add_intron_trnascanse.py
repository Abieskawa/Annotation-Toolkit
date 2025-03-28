#!/usr/bin/env python3
import sys
import argparse

def parse_attributes(attr_str):
    """Parse semicolon‐separated attributes into a dictionary."""
    attrs = {}
    for token in attr_str.split(';'):
        token = token.strip()
        if token and '=' in token:
            key, value = token.split('=', 1)
            attrs[key] = value
    return attrs

def compute_intron_from_exons(exon_lines):
    """
    Given a list of exon lines (each with 9 fields), compute the intron coordinates.
    Exons are sorted by their start (as integers). The intron is defined as the gap 
    between the end of the first exon and the start of the second exon.
    """
    exon_coords = []
    for line in exon_lines:
        fields = line.split("\t")
        start = int(fields[3])
        end = int(fields[4])
        exon_coords.append((start, end))
    exon_coords.sort(key=lambda x: x[0])
    intron_start = exon_coords[0][1] + 1
    intron_end = exon_coords[1][0] - 1
    return str(intron_start), str(intron_end)

def flush_group(group, output):
    """
    Output the current gene group.
    
    Group is a dictionary with keys:
      'gene'     : the gene (or pseudogene) line,
      'exons'    : list of exon lines belonging to that gene,
      'gene_type': type from the gene's attribute (e.g. "tRNA" or "pseudogene")
    
    If gene_type is "tRNA" and there are exactly 2 exon lines,
    compute an intron and insert it immediately after the gene line.
    """
    if not group:
        return
    gene_line = group['gene']
    gene_attrs = parse_attributes(gene_line.split("\t")[8])
    gene_id = gene_attrs.get("ID", "unknown")
    gene_type = group['gene_type'].lower()
    # Output the gene line first.
    output.append(gene_line)
    # For tRNA groups with exactly 2 exons, add an intron line.
    if gene_type == "trna" and len(group['exons']) == 2:
        intron_start, intron_end = compute_intron_from_exons(group['exons'])
        fields = gene_line.split("\t")
        intron_attrs = f"ID={gene_id}.intron;Parent={gene_id};"
        intron_line = "\t".join([
            fields[0],       # scaffold
            fields[1],       # source
            "intron",        # feature type
            intron_start,    # intron start
            intron_end,      # intron end
            ".",             # score
            fields[6],       # strand
            ".",             # phase
            intron_attrs     # attributes
        ])
        output.append(intron_line)
    # Output exon lines (in their original order)
    for exon in group['exons']:
        output.append(exon)

def process_lines(lines):
    """
    Process input lines, grouping gene (or pseudogene) and their exons together.
    
    For tRNA groups with exactly two exon lines, computes an intron line from the 
    gap between the exons and inserts it immediately after the gene line.
    Other records are output unchanged.
    """
    output = []
    current_group = None
    for line in lines:
        line = line.rstrip("\n")
        if not line.strip() or len(line.split("\t")) < 9:
            output.append(line)
            continue
        fields = line.split("\t")
        feature = fields[2].lower()
        attrs = parse_attributes(fields[8])
        # Start a new group when a gene (or pseudogene) record is encountered.
        if feature in ("trna", "pseudogene"):
            flush_group(current_group, output)
            current_group = {
                'gene': line,
                'exons': [],
                'gene_type': attrs.get("gene_biotype", feature)
            }
        elif feature == "exon":
            # If there's an active group and this exon belongs to it, add it.
            if current_group:
                parent = parse_attributes(fields[8]).get("Parent", "")
                gene_id = parse_attributes(current_group['gene'].split("\t")[8]).get("ID", "")
                if parent == gene_id:
                    current_group['exons'].append(line)
                else:
                    flush_group(current_group, output)
                    current_group = None
                    output.append(line)
            else:
                output.append(line)
        else:
            flush_group(current_group, output)
            current_group = None
            output.append(line)
    flush_group(current_group, output)
    return output

def main():
    parser = argparse.ArgumentParser(
        description="Add intron lines to tRNA groups that have exactly two exons. "
                    "Pseudogene groups are left unchanged."
    )
    parser.add_argument("-i", "--input", help="Input file (default: STDIN)", default=None)
    parser.add_argument("-o", "--output", help="Output file (default: STDOUT)", default=None)
    args = parser.parse_args()

    # Read input either from file or STDIN.
    if args.input:
        with open(args.input, 'r') as f:
            lines = f.readlines()
    else:
        lines = sys.stdin.readlines()
    
    processed = process_lines(lines)
    
    # Write output to file or STDOUT.
    if args.output:
        with open(args.output, 'w') as f:
            f.write("\n".join(processed) + "\n")
    else:
        for line in processed:
            print(line)

if __name__ == "__main__":
    main()
