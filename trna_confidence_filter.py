#!/usr/bin/env python3
"""
tRNA GFF Confidence Filter

This script reads EukConfidenceFilter output and adjusts tRNAscan-SE GFF files
based on confidence levels. Entries marked as "high confidence set" remain as tRNA,
while others are converted to pseudogenes.

Usage:
    python trna_confidence_filter.py -i input.gff -o output_prefix -c confidence.out [-d directory]
    
Example:
    python trna_confidence_filter.py -i genome.gff -o filtered_genome -c euk_confidence.out
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path


class tRNAConfidenceFilter:
    def __init__(self):
        self.high_confidence_entries = set()
        self.confidence_data = {}
        
    def parse_confidence_file(self, confidence_file):
        """
        Parse EukConfidenceFilter output file to extract high confidence entries
        """
        print(f"Reading confidence file: {confidence_file}")
        
        try:
            with open(confidence_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: Confidence file '{confidence_file}' not found.")
            sys.exit(1)
            
        # Skip header lines (usually first 3-4 lines)
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('--------'):
                data_start = i + 1
                break
                
        entries_processed = 0
        high_confidence_count = 0
        
        for line in lines[data_start:]:
            line = line.strip()
            if not line:
                continue
                
            # Split by tab - EukConfidenceFilter uses tab delimiters
            parts = line.split('\t')
            if len(parts) < 13:
                continue
                
            try:
                scaffold = parts[0].strip()
                trna_num = int(parts[1].strip())
                original_start = int(parts[2].strip())
                original_end = int(parts[3].strip())
                note = parts[13].strip() if len(parts) > 13 else ""
                
                # Normalize coordinates for matching (GFF always has start <= end)
                # Confidence file may have start > end for minus strand
                normalized_start = min(original_start, original_end)
                normalized_end = max(original_start, original_end)
                trna_key = (scaffold, normalized_start, normalized_end)
                
                # Store confidence data
                self.confidence_data[trna_key] = {
                    'scaffold': scaffold,
                    'trna_num': trna_num,
                    'start': normalized_start,
                    'end': normalized_end,
                    'original_start': original_start,
                    'original_end': original_end,
                    'note': note,
                    'is_high_confidence': 'high confidence set' in note
                }
                
                # Track high confidence entries
                if 'high confidence set' in note:
                    self.high_confidence_entries.add(trna_key)
                    high_confidence_count += 1
                    
                entries_processed += 1
                
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line: {line.strip()}")
                continue
                
        print(f"Processed {entries_processed} confidence entries")
        print(f"Found {high_confidence_count} high confidence tRNAs")
        
        return entries_processed > 0
    
    def process_gff_file(self, input_gff, output_prefix, output_dir=None):
        """
        Process GFF file and adjust tRNA/pseudogene annotations based on confidence
        """
        input_path = Path(input_gff)
        
        if not input_path.exists():
            print(f"Error: Input GFF file '{input_gff}' not found.")
            return False
            
        # Determine output file path
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = output_dir / f"{output_prefix}.gff"
        else:
            output_file = Path(f"{output_prefix}.gff")
            
        print(f"Processing GFF file: {input_gff}")
        print(f"Output file: {output_file}")
        
        modifications_made = 0
        total_trna_entries = 0
        
        try:
            with open(input_gff, 'r') as infile, open(output_file, 'w') as outfile:
                for line_num, line in enumerate(infile, 1):
                    original_line = line
                    line = line.strip()
                    
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        outfile.write(original_line)
                        continue
                        
                    # Parse GFF line
                    parts = line.split('\t')
                    if len(parts) != 9:
                        outfile.write(original_line)
                        continue
                        
                    scaffold, source, feature_type, start, end, score, strand, phase, attributes = parts
                    
                    # Only process tRNA and pseudogene entries
                    if feature_type not in ['tRNA', 'pseudogene']:
                        outfile.write(original_line)
                        continue
                        
                    if feature_type in ['tRNA', 'pseudogene']:
                        total_trna_entries += 1
                        
                    try:
                        start_pos = int(start)
                        end_pos = int(end)
                        # GFF coordinates are always start <= end, normalize for matching
                        trna_key = (scaffold, start_pos, end_pos)
                        
                        # Check if this entry has confidence data
                        if trna_key in self.confidence_data:
                            confidence_info = self.confidence_data[trna_key]
                            
                            if confidence_info['is_high_confidence']:
                                # High confidence: should be tRNA
                                if feature_type != 'tRNA':
                                    parts[2] = 'tRNA'
                                    # Update gene_biotype in attributes
                                    parts[8] = self.update_gene_biotype(attributes, 'tRNA')
                                    modifications_made += 1
                                    print(f"  Line {line_num}: Changed {feature_type} to tRNA for {scaffold}:{start}-{end} (high confidence set)")
                            else:
                                # Not high confidence: should be pseudogene
                                if feature_type != 'pseudogene':
                                    parts[2] = 'pseudogene'
                                    # Update gene_biotype in attributes
                                    parts[8] = self.update_gene_biotype(attributes, 'pseudogene')
                                    modifications_made += 1
                                    reason = confidence_info['note']
                                    print(f"  Line {line_num}: Changed {feature_type} to pseudogene for {scaffold}:{start}-{end} (reason: {reason})")
                        else:
                            # No confidence data found
                            print(f"  Warning: No confidence data found.")
                            
                        # Write the (possibly modified) line
                        outfile.write('\t'.join(parts) + '\n')
                        
                    except ValueError as e:
                        print(f"Warning: Could not parse coordinates on line {line_num}: {e}")
                        outfile.write(original_line)
                        continue
                        
        except Exception as e:
            print(f"Error processing GFF file: {e}")
            return False
            
        print(f"\nProcessing complete!")
        print(f"Total tRNA/pseudogene entries processed: {total_trna_entries}")
        print(f"Modifications made: {modifications_made}")
        print(f"Output written to: {output_file}")
        
        return True
    
    def update_gene_biotype(self, attributes, new_biotype):
        """
        Update the gene_biotype attribute in GFF attributes string
        """
        # Split attributes by semicolon
        attr_parts = attributes.split(';')
        updated_parts = []
        biotype_found = False
        
        for part in attr_parts:
            part = part.strip()
            if part.startswith('gene_biotype='):
                # Replace the existing gene_biotype with new value
                updated_parts.append(f'gene_biotype={new_biotype}')
                biotype_found = True
            elif part:  # Only add non-empty parts
                updated_parts.append(part)
                
        # If gene_biotype wasn't found, add it
        if not biotype_found:
            updated_parts.append(f'gene_biotype={new_biotype}')
            
        return ';'.join(updated_parts)


def main():
    parser = argparse.ArgumentParser(
        description='Filter tRNA GFF annotations based on EukConfidenceFilter confidence levels',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python trna_confidence_filter.py -i genome.gff -c confidence.out -o filtered_genome
  
  # Specify output directory
  python trna_confidence_filter.py -i genome.gff -c confidence.out -o filtered_genome -d results/
  
  # Process multiple files in a directory
  python trna_confidence_filter.py -i *.gff -c confidence.out -o filtered -d output/
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input GFF file from tRNAscan-SE')
    parser.add_argument('-c', '--confidence', required=True,
                        help='EukConfidenceFilter output file (.out)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file prefix (without .gff extension)')
    parser.add_argument('-d', '--directory', default=None,
                        help='Output directory (optional)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Initialize the filter
    filter_tool = tRNAConfidenceFilter()
    
    # Parse confidence file
    if not filter_tool.parse_confidence_file(args.confidence):
        print("Error: Failed to parse confidence file")
        sys.exit(1)
    
    # Process GFF file
    if not filter_tool.process_gff_file(args.input, args.output, args.directory):
        print("Error: Failed to process GFF file")
        sys.exit(1)
    
    print("\n" + "="*50)
    print("SUCCESS: tRNA confidence filtering completed!")
    print("="*50)


if __name__ == '__main__':
    main()