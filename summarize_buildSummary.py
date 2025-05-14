#!/usr/bin/env python3
import re
import sys

def parse_total_length(lines):
    """
    Look for a line like "Total Length: 1085333649 bp" and return the integer.
    """
    total_length = None
    for line in lines:
        m = re.search(r"Total Length:\s*([\d,]+)\s*bp", line)
        if m:
            total_length = int(m.group(1).replace(',', ''))
            break
    return total_length

def parse_table(lines):
    """
    Extract the table lines after the header "Class   Count   bpMasked   %masked".
    We assume that the table starts after this header line.
    """
    table_lines = []
    header_found = False
    for line in lines:
        if not header_found:
            if re.search(r"^Class\s+Count\s+bpMasked\s+%masked", line):
                header_found = True
            continue
        # Skip dashed lines or empty lines.
        if re.match(r"^[-=]+", line.strip()) or not line.strip():
            continue
        table_lines.append(line.rstrip())
    return table_lines

def process_table(table_lines, total_length):
    """
    Process the table lines into a summary dictionary.
    
    For each top-level category (lines with no indentation), record its top-level values.
    For each indented child line, add its values to the parent's children sum.
    
    The final value is the sum of top-level count and children count, and similarly for bpMasked.
    """
    summary = {}
    current_category = None

    for line in table_lines:
        # Check if the line is top-level (no leading whitespace)
        if line and not line[0].isspace():
            parts = re.split(r"\s{2,}", line.strip())
            if len(parts) >= 4:
                cat = parts[0]
                try:
                    top_count = int(parts[1])
                except ValueError:
                    top_count = 0
                try:
                    top_bp = int(parts[2])
                except ValueError:
                    top_bp = 0
                summary[cat] = {
                    "top_count": top_count,
                    "top_bp": top_bp,
                    "children_count": 0,
                    "children_bp": 0
                }
                current_category = cat
            else:
                current_category = None
        else:
            # This is an indented child line.
            if current_category is not None:
                parts = re.split(r"\s{2,}", line.strip())
                if len(parts) >= 4:
                    try:
                        child_count = int(parts[1])
                    except ValueError:
                        child_count = 0
                    try:
                        child_bp = int(parts[2])
                    except ValueError:
                        child_bp = 0
                    summary[current_category]["children_count"] += child_count
                    summary[current_category]["children_bp"] += child_bp

    final_summary = {}
    for cat, data in summary.items():
        total_count = data["top_count"] + data["children_count"]
        total_bp = data["top_bp"] + data["children_bp"]
        ratio = (total_bp / total_length * 100) if total_length and total_length > 0 else 0
        final_summary[cat] = {"count": total_count, "bp": total_bp, "ratio": ratio}
    return final_summary

def main():
    if len(sys.argv) < 2:
        print("Usage: {} <summary_repeat_table.txt from buildSummary.pl>".format(sys.argv[0]))
        sys.exit(1)
    filename = sys.argv[1]
    with open(filename, "r") as f:
        lines = f.readlines()

    total_length = parse_total_length(lines)
    if total_length is None:
        print("Total Length not found in file. Please ensure the header contains it.")
        sys.exit(1)

    table_lines = parse_table(lines)
    summary = process_table(table_lines, total_length)

    print("Summary of Top-Level Categories (Total Length: {} bp):".format(total_length))
    print("{:<15} {:>10} {:>15} {:>10}".format("Category", "Count", "bpMasked", "%masked"))
    print("-" * 55)
    for cat, data in summary.items():
        print("{:<15} {:>10} {:>15} {:>9.2f}%".format(
            cat, data["count"], data["bp"], data["ratio"]
        ))

if __name__ == "__main__":
    main()
