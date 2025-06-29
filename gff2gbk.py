#!/usr/bin/env python3
"""
RefSeq-style GFF + FASTA to GenBank converter
Produces publication-quality GenBank files following RefSeq standards
Usage: python refseq_gff2gbk.py fasta_file gff_file output_file organism [division]
"""

import sys
import os
from datetime import datetime
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Standard mitochondrial gene mappings following RefSeq conventions
MITO_GENE_MAP = {
    # NADH dehydrogenase genes
    'ND1': {'gene': 'ND1', 'product': 'NADH dehydrogenase subunit 1'},
    'ND2': {'gene': 'ND2', 'product': 'NADH dehydrogenase subunit 2'},
    'ND3': {'gene': 'ND3', 'product': 'NADH dehydrogenase subunit 3'},
    'ND4': {'gene': 'ND4', 'product': 'NADH dehydrogenase subunit 4'},
    'ND4L': {'gene': 'ND4L', 'product': 'NADH dehydrogenase subunit 4L'},
    'ND5': {'gene': 'ND5', 'product': 'NADH dehydrogenase subunit 5'},
    'ND6': {'gene': 'ND6', 'product': 'NADH dehydrogenase subunit 6'},
    
    # Cytochrome oxidase genes
    'COX1': {'gene': 'COX1', 'product': 'cytochrome c oxidase subunit I'},
    'COX2': {'gene': 'COX2', 'product': 'cytochrome c oxidase subunit II'},
    'COX3': {'gene': 'COX3', 'product': 'cytochrome c oxidase subunit III'},
    'COI': {'gene': 'COX1', 'product': 'cytochrome c oxidase subunit I'},
    'COII': {'gene': 'COX2', 'product': 'cytochrome c oxidase subunit II'},
    'COIII': {'gene': 'COX3', 'product': 'cytochrome c oxidase subunit III'},
    
    # Cytochrome b
    'CYTB': {'gene': 'CYTB', 'product': 'cytochrome b'},
    'COB': {'gene': 'CYTB', 'product': 'cytochrome b'},
    
    # ATP synthase
    'ATP6': {'gene': 'ATP6', 'product': 'ATP synthase F0 subunit 6'},
    'ATP8': {'gene': 'ATP8', 'product': 'ATP synthase F0 subunit 8'},
    
    # Ribosomal RNAs - FIXED: gene should be same as note
    'RRN12': {'gene': '12S ribosomal RNA', 'product': 's-rRNA', 'note': '12S ribosomal RNA'},
    'RRN16': {'gene': '16S ribosomal RNA', 'product': 'l-rRNA', 'note': '16S ribosomal RNA'},
    'RRNS': {'gene': '12S ribosomal RNA', 'product': 's-rRNA', 'note': '12S ribosomal RNA'},
    'RRNL': {'gene': '16S ribosomal RNA', 'product': 'l-rRNA', 'note': '16S ribosomal RNA'},
    '12S': {'gene': '12S ribosomal RNA', 'product': 's-rRNA', 'note': '12S ribosomal RNA'},
    '16S': {'gene': '16S ribosomal RNA', 'product': 'l-rRNA', 'note': '16S ribosomal RNA'},
    'S-RRNA': {'gene': '12S ribosomal RNA', 'product': 's-rRNA', 'note': '12S ribosomal RNA'},
    'L-RRNA': {'gene': '16S ribosomal RNA', 'product': 'l-rRNA', 'note': '16S ribosomal RNA'},
}

# tRNA amino acid mappings - MODIFIED to match RefSeq format
TRNA_MAP = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

def parse_gene_info(feature_type, qualifiers):
    """Parse gene information following RefSeq standards"""
    
    # Get name from various qualifier fields
    name = None
    for key in ['gene', 'Name', 'ID', 'product']:
        if key in qualifiers:
            name = qualifiers[key]
            if isinstance(name, list):
                name = name[0]
            break
    
    if not name:
        name = feature_type
    
    name = name.upper().strip()
    
    # Handle standard mitochondrial genes
    if name in MITO_GENE_MAP:
        return MITO_GENE_MAP[name]
    
    # Handle tRNA genes - MODIFIED LOGIC: gene and product should be same
    if 'TRN' in name or 'TRNA' in name:
        # Extract amino acid and codon from various formats
        if '(' in name:
            # Format: trnL2(taa) or TRN-A(UGC)
            aa_part = name.split('(')[0]
            codon_part = name.split('(')[1].split(')')[0].upper()
            
            # Extract amino acid letter
            if 'TRN' in aa_part:
                aa = aa_part.replace('TRN', '').replace('TRNA', '').replace('-', '').strip()
                # Remove numbers (like L2 -> L)
                aa = ''.join([c for c in aa if c.isalpha()])
                if len(aa) >= 1:
                    aa = aa[0]  # Take first letter
                    if aa in TRNA_MAP:
                        trna_name = f'tRNA-{TRNA_MAP[aa]}'
                        return {'gene': trna_name, 'product': trna_name}
        elif '-' in name:
            # Format: TRNA-ALA
            parts = name.split('-')
            if len(parts) > 1:
                aa = parts[1][:3]
                if aa in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                         'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                    trna_name = f'tRNA-{aa.capitalize()}'
                    return {'gene': trna_name, 'product': trna_name}
    
    # Default fallback
    return {'gene': name, 'product': name}

def create_refseq_features(gff_file, seq_record):
    """Create RefSeq-style features from GFF"""
    
    features = []
    
    with open(gff_file, 'r') as handle:
        for line_num, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 8:
                continue
            
            try:
                seqid, source, feature_type, start, end, score, strand, phase = parts[:8]
                attributes = parts[8] if len(parts) > 8 else ""
                
                # Only process features matching our sequence
                if seqid != seq_record.id:
                    continue
                
                # Parse attributes
                attr_dict = {}
                if attributes:
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key.strip()] = value.strip().strip('"')
                
                # Create location
                start_pos = int(start) - 1
                end_pos = int(end)
                strand_val = -1 if strand == '-' else 1
                location = FeatureLocation(start_pos, end_pos, strand=strand_val)
                
                # Get gene information
                gene_info = parse_gene_info(feature_type, attr_dict)
                
                # Create appropriate qualifiers based on feature type
                qualifiers = {}
                
                if feature_type.lower() == 'gene':
                    if gene_info['gene']:
                        qualifiers['gene'] = gene_info['gene']
                    
                elif feature_type.lower() == 'cds':
                    if gene_info['gene']:
                        qualifiers['gene'] = gene_info['gene']
                    qualifiers['codon_start'] = '1'
                    qualifiers['transl_table'] = '2'  # Mitochondrial genetic code
                    qualifiers['product'] = gene_info['product']
                    
                    # Add note for special cases
                    if 'partial' in str(attr_dict.get('Note', '')).lower():
                        if start_pos == 0:  # Starts at beginning
                            qualifiers['note'] = 'start codon not determined'
                        qualifiers['transl_except'] = f'(pos:{end_pos},aa:TERM)'
                        qualifiers['note'] = qualifiers.get('note', '') + '; TAA stop codon is completed by the addition of 3\' A residues to the mRNA'
                    
                elif feature_type.lower() == 'trna':
                    # FIXED: both gene and product should be same for tRNA
                    if gene_info['gene']:
                        qualifiers['gene'] = gene_info['gene']
                    qualifiers['product'] = gene_info['product']
                    
                    # Add codon recognition if available - MODIFIED
                    if 'codon_recognized' in attr_dict:
                        qualifiers['codon_recognized'] = attr_dict['codon_recognized']
                    elif '(' in str(attr_dict.get('Name', '')):
                        # Extract from name like trnL(taa)
                        name_val = attr_dict['Name']
                        if '(' in name_val:
                            codon = name_val.split('(')[1].split(')')[0].upper()
                            if len(codon) == 3:
                                qualifiers['codon_recognized'] = codon
                    
                elif feature_type.lower() == 'rrna':
                    # FIXED: gene should be filled with note content
                    if gene_info['gene']:
                        qualifiers['gene'] = gene_info['gene']
                    qualifiers['product'] = gene_info['product']
                    if 'note' in gene_info:
                        qualifiers['note'] = gene_info['note']
                
                else:
                    # For other feature types, preserve original qualifiers
                    qualifiers.update(attr_dict)
                    if gene_info['gene']:
                        qualifiers['gene'] = gene_info['gene']
                    if 'product' not in qualifiers:
                        qualifiers['product'] = gene_info['product']
                
                # Create feature
                feature = SeqFeature(
                    location=location,
                    type=feature_type,
                    qualifiers=qualifiers
                )
                features.append(feature)
                
            except Exception as e:
                print(f"Warning: Error processing line {line_num}: {e}")
                continue
    
    return features

def create_refseq_genbank(fasta_file, gff_file, output_file, organism, division="VRT"):
    """Create RefSeq-style GenBank file"""
    
    # Load sequence
    with open(fasta_file, 'r') as handle:
        sequences = list(SeqIO.parse(handle, 'fasta'))
    
    if not sequences:
        print("Error: No sequences found in FASTA file")
        return
    
    if len(sequences) > 1:
        print(f"Found {len(sequences)} sequences:")
        for i, seq in enumerate(sequences):
            print(f"  {i+1}. {seq.id} ({len(seq.seq):,} bp)")
        
        choice = input("Enter sequence number to use (or press Enter for first): ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(sequences):
            seq_record = sequences[int(choice) - 1]
        else:
            seq_record = sequences[0]
        print(f"Using: {seq_record.id}")
    else:
        seq_record = sequences[0]
    
    print(f"Sequence: {seq_record.id} ({len(seq_record.seq):,} bp)")
    
    # Create features
    features = create_refseq_features(gff_file, seq_record)
    print(f"Created {len(features)} features")
    
    # Detect genome type
    is_mitochondrial = (
        len(seq_record.seq) < 25000 and
        any('ND' in str(f.qualifiers.get('gene', '')) or 
            'COX' in str(f.qualifiers.get('gene', '')) for f in features)
    )
    
    # Create GenBank record
    gb_record = SeqRecord(
        seq_record.seq,
        id=seq_record.id,
        name=seq_record.id[:10],  # GenBank name limit
        description=f"{organism} mitochondrion, complete genome" if is_mitochondrial else f"{organism}",
        features=features
    )
    
    # Set RefSeq-style annotations
    gb_record.annotations = {
        "molecule_type": "DNA",
        "organism": organism,
        "topology": "circular" if is_mitochondrial else "linear",
        "data_file_division": division,
        "date": datetime.now().strftime("%d-%b-%Y").upper(),
        "accessions": [seq_record.id],
        "sequence_version": 1,
        "keywords": ["RefSeq"] if not is_mitochondrial else [],
        "source": f"mitochondrion {organism}" if is_mitochondrial else organism,
    }
    
    # Add taxonomy (basic eukaryotic hierarchy)
    if is_mitochondrial:
        if division == "VRT":  # Vertebrate
            gb_record.annotations["taxonomy"] = [
                "Eukaryota", "Metazoa", "Chordata", "Craniata", "Vertebrata", 
                "Euteleostomi", "Actinopterygii"
            ]
        elif division == "INV":  # Invertebrate
            gb_record.annotations["taxonomy"] = [
                "Eukaryota", "Metazoa"
            ]
        else:
            gb_record.annotations["taxonomy"] = ["Eukaryota"]
    else:
        gb_record.annotations["taxonomy"] = ["Eukaryota"]
    
    # Add source feature if not present
    has_source = any(f.type == "source" for f in gb_record.features)
    if not has_source:
        source_qualifiers = {
            "organism": organism,
            "mol_type": "genomic DNA"
        }
        
        if is_mitochondrial:
            source_qualifiers["organelle"] = "mitochondrion"
        
        source_feature = SeqFeature(
            location=FeatureLocation(0, len(gb_record.seq)),
            type="source",
            qualifiers=source_qualifiers
        )
        gb_record.features.insert(0, source_feature)
    
    # Write GenBank file
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as handle:
            SeqIO.write(gb_record, handle, 'genbank')
        
        print(f"\nRefSeq-style GenBank file created: {output_file}")
        print(f"Features: {len(gb_record.features)}")
        print(f"Topology: {gb_record.annotations['topology']}")
        print(f"Division: {gb_record.annotations['data_file_division']}")
        
        # Feature summary
        feature_counts = {}
        for f in gb_record.features:
            ftype = f.type
            feature_counts[ftype] = feature_counts.get(ftype, 0) + 1
        
        print("\nFeature summary:")
        for ftype, count in sorted(feature_counts.items()):
            print(f"  {ftype}: {count}")
        
    except Exception as e:
        print(f"Error writing GenBank file: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) < 4:
        print("Usage: python refseq_gff2gbk.py <fasta> <gff> <output> <organism> [division]")
        print("\nDivision codes:")
        print("  VRT = Vertebrate (default for mitochondria)")
        print("  INV = Invertebrate")
        print("  PLN = Plant") 
        print("  BCT = Bacterial")
        print("  UNC = Unclassified")
        print("\nExample:")
        print("  python refseq_gff2gbk.py seq.fa anno.gff result.gb 'Trichiurus japonicus' VRT")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    gff_file = sys.argv[2]
    output_file = sys.argv[3]
    organism = sys.argv[4]
    division = sys.argv[5] if len(sys.argv) > 5 else "VRT"
    
    # Check files exist
    for f in [fasta_file, gff_file]:
        if not os.path.exists(f):
            print(f"Error: File not found: {f}")
            sys.exit(1)
    
    print(f"Creating RefSeq-style GenBank file for: {organism}")
    print(f"Division: {division}")
    print("-" * 50)
    
    create_refseq_genbank(fasta_file, gff_file, output_file, organism, division)

if __name__ == "__main__":
    main()