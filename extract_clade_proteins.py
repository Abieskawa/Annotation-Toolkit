#!/usr/bin/env python3
import argparse
import gzip
import re
import sys

def get_taxid(clade, levels_file):
    """
    Given a clade name, find its NCBI tax id from the levels file.
    
    Expected file (tab-delimited):
      Column 1: Level NCBI tax id
      Column 2: Scientific name (clade)
      Column 5: Total non-redundant count of species underneath
    
    Returns:
      (taxid, expected_species_count) or (None, None) if not found.
    """
    with gzip.open(levels_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            if fields[1] == clade:
                return fields[0], fields[4] if len(fields) >= 5 else None
    return None, None

def get_species_ids(taxid, level2species_file):
    """
    Extract species IDs from the level2species file that belong to the given taxid.
    
    Expected file (tab-delimited):
      Column 1: Top-most level NCBI tax id (ignored here)
      Column 2: Ortho DB organism id (species id)
      Column 4: Ordered list of intermediate levels (e.g., "{2,1052212}")
    
    The taxid is looked for in the list from column 4.
    """
    species_ids = []
    with gzip.open(level2species_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            levels = fields[3].strip('{}').split(',')
            if taxid in levels:
                species_ids.append(fields[1])
    return species_ids

def report_species(species_ids, species_map_file=None):
    """
    Report the species found. If a species mapping file is provided,
    also include the scientific name.
    
    The species mapping file (tab-delimited) is expected to have:
      Column 2: Ortho DB organism id
      Column 3: Scientific name
    """
    species_info = {}
    if species_map_file:
        open_func = gzip.open if species_map_file.endswith('.gz') else open
        with open_func(species_map_file, 'rt') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                species_info[fields[1]] = fields[2]
    print("Species found under the clade:", file=sys.stderr)
    for s in species_ids:
        if species_map_file:
            name = species_info.get(s, "N/A")
            print(f"{s}\t{name}", file=sys.stderr)
        else:
            print(s, file=sys.stderr)

def extract_fasta(species_ids, fasta_file, output_file=None):
    """
    Extract FASTA records whose header contains any of the species IDs.
    
    Assumes FASTA records are two lines each: a header (starting with '>') and a sequence.
    """
    regex = r'\b(?:' + '|'.join(re.escape(s) for s in species_ids) + r')\b'
    pattern = re.compile(regex)
    
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    outfh = sys.stdout if output_file is None else open(output_file, 'w')
    
    with open_func(fasta_file, 'rt') as f:
        header = None
        sequence = None
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header and pattern.search(header):
                    outfh.write(header + '\n')
                    outfh.write(sequence + '\n')
                header = line
                sequence = ""
            else:
                sequence = line
        if header and pattern.search(header):
            outfh.write(header + '\n')
            outfh.write(sequence + '\n')
    
    if output_file is not None:
        outfh.close()

def main():
    parser = argparse.ArgumentParser(
        description="Extract protein sequences for a given clade from OrthoDB files and report the species found."
    )
    parser.add_argument('--clade', required=True, help="Clade name (e.g., Apicomplexa)")
    parser.add_argument('--levels', required=True, help="Path to levels file (e.g., odb12v0_levels.tab.gz)")
    parser.add_argument('--level2species', required=True, help="Path to level2species file (e.g., odb12v0_level2species.tab.gz)")
    parser.add_argument('--fasta', required=True, help="Path to FASTA file (e.g., odb12v0_aa_fasta.gz)")
    parser.add_argument('--output', default=None, help="Output file for extracted FASTA (default: stdout)")
    parser.add_argument('--report-species', action="store_true", help="Report species underneath the clade (to stderr)")
    parser.add_argument('--speciesmap', default=None, help="Path to species mapping file (e.g., odb12v0_species.tab or odb12v0_species.tab.gz) to include scientific names")
    args = parser.parse_args()

    # Retrieve taxonomic id and expected species count for the clade.
    taxid, expected_species_count = get_taxid(args.clade, args.levels)
    if taxid is None:
        print(f"Error: Clade '{args.clade}' not found in levels file.", file=sys.stderr)
        sys.exit(1)
    
    # Extract species IDs.
    species_ids = get_species_ids(taxid, args.level2species)
    if expected_species_count is not None:
        try:
            expected = int(expected_species_count)
            if len(species_ids) != expected:
                print(f"Warning: Expected {expected} species, but found {len(species_ids)}.", file=sys.stderr)
        except ValueError:
            pass
    print(f"Found {len(species_ids)} species IDs for clade '{args.clade}' (taxid: {taxid}).", file=sys.stderr)
    
    # Optionally report species (with scientific names if a species map is provided).
    if args.report_species:
        report_species(species_ids, args.speciesmap)

    # Extract FASTA records matching the species IDs.
    extract_fasta(species_ids, args.fasta, args.output)

if __name__ == "__main__":
    main()
