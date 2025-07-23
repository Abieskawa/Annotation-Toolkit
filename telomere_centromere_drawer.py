#!/usr/bin/env python3
"""
Telomere Drawer - Extracted from quarTeT
A standalone script for drawing chromosome ideograms with telomere information
Enhanced with centromere support
"""

import sys
import os
import argparse
import subprocess
import re
import math
from collections import defaultdict


def readFastaAsDict(fastafile):
    """Read FASTA file and return as dictionary"""
    fastaDict = {}
    with open(fastafile, 'r') as fil:
        allline = fil.read()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict


def parse_centromere_file(centrofile):
    """
    Parse centromere file in format: chromosome,start-end[,color]
    Example: scaffold1,37651276-38156789,ff0000
    Color is optional (defaults to red)
    """
    print(f'[Info] Parsing centromere file: {centrofile}')
    centromeres = {}
    
    with open(centrofile, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            try:
                parts = line.split(',')
                if len(parts) < 2:
                    print(f'[Warning] Line {line_num}: Invalid format, expected "chromosome,start-end[,color]": {line}')
                    continue
                
                chrom = parts[0].strip()
                coords = parts[1].strip()
                
                # Parse coordinates
                if '-' not in coords:
                    print(f'[Warning] Line {line_num}: Invalid coordinates format, expected "start-end": {coords}')
                    continue
                    
                start_str, end_str = coords.split('-', 1)
                start = int(start_str.strip())
                end = int(end_str.strip())
                
                if start >= end:
                    print(f'[Warning] Line {line_num}: Invalid range, start >= end: {start}-{end}')
                    continue
                
                # Parse color (optional, default to red)
                color = 'ff0000'  # Default red
                if len(parts) >= 3 and parts[2].strip():
                    color = parts[2].strip()
                    # Remove # if present
                    if color.startswith('#'):
                        color = color[1:]
                    # Validate hex color
                    if not all(c in '0123456789abcdefABCDEF' for c in color) or len(color) != 6:
                        print(f'[Warning] Line {line_num}: Invalid hex color "{color}", using default red')
                        color = 'ff0000'
                
                centromeres[chrom] = (start, end, color)
                print(f'[Info] {chrom}: centromere at {start}-{end} (color: #{color}, length: {end-start+1} bp)')
                
            except ValueError as e:
                print(f'[Warning] Line {line_num}: Could not parse coordinates: {line} - {e}')
                continue
            except Exception as e:
                print(f'[Warning] Line {line_num}: Error parsing line: {line} - {e}')
                continue
    
    print(f'[Info] Successfully parsed {len(centromeres)} centromeres')
    return centromeres


def parse_tidk_telomere_file(tidk_file, min_repeat_threshold=100):
    """
    Parse tidk telomeric_repeat_windows.tsv file and convert to telomere info
    Uses quarTeT's exact logic for telomere detection
    """
    print(f'[Info] Parsing tidk telomere file: {tidk_file}')
    print(f'[Info] Using quarTeT telomere detection logic with threshold: {min_repeat_threshold}')
    
    # Read tidk file - same as quarTeT
    block = defaultdict(list)
    with open(tidk_file, 'r') as f:
        for line in f:
            if line.startswith('id'):  # Skip header
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrid, window, forward_repeats, reverse_repeats = parts[:4]
                block[chrid].append([int(forward_repeats), int(reverse_repeats)])
    
    # Process telomere information - EXACT quarTeT logic
    telodict = {}
    for chrid in block:
        if not block[chrid]:  # Skip empty chromosomes
            continue
        
        # quarTeT logic: determine terminal region size
        sidesize = math.floor(len(block[chrid])/2) if len(block[chrid]) < 30 else 15
        print(f'[Info] {chrid}: {len(block[chrid])} windows, analyzing {sidesize} windows at each end')
        
        # LEFT telomere analysis
        lefttelo = block[chrid][0:sidesize]
        totalforwardrepeatnum = 0
        totalreverserepeatnum = 0
        for forwardrepeatnum, reverserepeatnum in lefttelo:
            totalforwardrepeatnum += forwardrepeatnum
            totalreverserepeatnum += reverserepeatnum
        
        if max(totalforwardrepeatnum, totalreverserepeatnum) >= min_repeat_threshold:
            if totalforwardrepeatnum >= totalreverserepeatnum:
                telodict[f'{chrid}.L'] = [totalforwardrepeatnum, '+']
                print(f'[Info] {chrid} LEFT telomere: {totalforwardrepeatnum} repeats (+)')
            else:
                telodict[f'{chrid}.L'] = [totalreverserepeatnum, '-']
                print(f'[Info] {chrid} LEFT telomere: {totalreverserepeatnum} repeats (-)')
        else:
            print(f'[Info] {chrid} LEFT: {max(totalforwardrepeatnum, totalreverserepeatnum)} repeats < {min_repeat_threshold} (no telomere)')
        
        # RIGHT telomere analysis
        righttelo = block[chrid][len(block[chrid])-sidesize:]
        totalforwardrepeatnum = 0
        totalreverserepeatnum = 0
        for forwardrepeatnum, reverserepeatnum in righttelo:
            totalforwardrepeatnum += forwardrepeatnum
            totalreverserepeatnum += reverserepeatnum
        
        if max(totalforwardrepeatnum, totalreverserepeatnum) >= min_repeat_threshold:
            if totalforwardrepeatnum >= totalreverserepeatnum:
                telodict[f'{chrid}.R'] = [totalforwardrepeatnum, '+']
                print(f'[Info] {chrid} RIGHT telomere: {totalforwardrepeatnum} repeats (+)')
            else:
                telodict[f'{chrid}.R'] = [totalreverserepeatnum, '-']
                print(f'[Info] {chrid} RIGHT telomere: {totalreverserepeatnum} repeats (-)')
        else:
            print(f'[Info] {chrid} RIGHT: {max(totalforwardrepeatnum, totalreverserepeatnum)} repeats < {min_repeat_threshold} (no telomere)')
    
    return telodict, block


def convert_tidk_to_telo_info(tidk_file, fasta_file, output_file, min_repeat_threshold=100):
    """
    Convert tidk telomeric_repeat_windows.tsv to quarTeT .telo.info format
    Uses exact quarTeT telomere detection logic
    """
    telodict, block = parse_tidk_telomere_file(tidk_file, min_repeat_threshold)
    fasta = readFastaAsDict(fasta_file)
    
    with open(output_file, 'w') as t:
        both, side, no = 0, 0, 0
        status = {}
        
        for chrid in block:
            if f'{chrid}.L' in telodict and f'{chrid}.R' in telodict:
                both += 1
                status[chrid] = 'both'
            elif f'{chrid}.L' in telodict and f'{chrid}.R' not in telodict:
                side += 1
                status[chrid] = 'left'
            elif f'{chrid}.R' in telodict and f'{chrid}.L' not in telodict:
                side += 1
                status[chrid] = 'right'
            else:
                no += 1
                status[chrid] = 'no'
        
        t.write(f'# Telomere data converted from tidk format using quarTeT logic\n')
        t.write(f'# Minimum repeat threshold: {min_repeat_threshold}\n')
        t.write(f'# Both telomere found: {both}\n')
        t.write(f'# Only one telomere found: {side}\n')
        t.write(f'# No telomere found: {no}\n')
        t.write('# Chrid\tchrlen\tstatus\tleftnum\tleftdirection\trightnum\trightdirection\n')
        
        for chrid in block:
            if chrid not in fasta:
                print(f'[Warning] Chromosome {chrid} found in tidk file but not in FASTA')
                continue
                
            leftid = f'{chrid}.L'
            if leftid in telodict:
                leftinfo = f'{telodict[leftid][0]}\t{telodict[leftid][1]}'
            else:
                leftinfo = '0\t'
                
            rightid = f'{chrid}.R'
            if rightid in telodict:
                rightinfo = f'{telodict[rightid][0]}\t{telodict[rightid][1]}'
            else:
                rightinfo = '0\t'
                
            t.write(f'{chrid}\t{len(fasta[chrid])}\t{status[chrid]}\t{leftinfo}\t{rightinfo}\n')
    
    print(f'[Info] Converted tidk format to telo.info format: {output_file}')
    return output_file


def create_simple_agp_from_fasta(infile, outagp):
    """Create a simple AGP file from FASTA without gap information"""
    os.makedirs(os.path.dirname(outagp), exist_ok=True)
    
    fastadict = readFastaAsDict(infile)
    with open(outagp, 'w') as w:
        for sid, seq in fastadict.items():
            w.write(f'{sid}\t1\t{len(seq)}\t1\tW\t.\t.\t.\t.\n')


def agpGap(infile, outagp):
    """Generate AGP file from FASTA file, identifying gaps (N stretches)"""
    os.makedirs(os.path.dirname(outagp), exist_ok=True)
    
    fastadict = readFastaAsDict(infile)
    with open(outagp, 'w') as w:
        for sid, seq in fastadict.items():
            gapsitelist = [r.span() for r in re.finditer(r'N+', seq)]
            if gapsitelist == []:
                w.write(f'{sid}\t1\t{len(seq)}\t1\tW\t.\t.\t.\t.\n')
            else:
                sitelist = [(1, gapsitelist[0][0]), (gapsitelist[0][0]+1, gapsitelist[0][1])]
                for i in range(len(gapsitelist)):
                    if i + 1 == len(gapsitelist):
                        sitelist += [(gapsitelist[i][1]+1, len(seq))]
                    else:
                        sitelist += [(gapsitelist[i][1]+1, gapsitelist[i+1][0]), (gapsitelist[i+1][0]+1, gapsitelist[i+1][1])]
                count = 1
                for site in sitelist:
                    ty = 'W' if count % 2 == 1 else 'U'
                    w.write(f'{sid}\t{site[0]}\t{site[1]}\t{count}\t{ty}\t.\t.\t.\t.\n')
                    count += 1


def drawTelomereIdeogram(agpfile, outprefix, centrofile=None, telofile=None, show_gaps=False, show_position_markers=True):
    """
    Draw chromosome ideogram with telomere and colored centromere information
    
    Parameters:
    - agpfile: AGP file path
    - outprefix: Output file prefix
    - centrofile: Centromere file path (optional)
    - telofile: Telomere file path (optional)
    - show_gaps: Whether to show gaps in the ideogram (default: False)
    - show_position_markers: Whether to show 5 Mb position markers (default: True)
    """
    # Create temp directory
    subprocess.run('mkdir -p tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    # Use default files if they exist
    if centrofile is None and os.path.exists(f'{outprefix}.best.candidate'):
        centrofile = f'{outprefix}.best.candidate'
    if telofile is None and os.path.exists(f'{outprefix}.telo.info'):
        telofile = f'{outprefix}.telo.info'
    
    # Parse AGP file
    with open(agpfile, 'r') as agp:
        agpblock = defaultdict(list)
        chrnumber = 0
        for line in agp:
            if chrnumber == 0 or line.split()[0] != chrid:
                chrnumber += 1
            chrid, seqstart, seqend, num, attr, tigid, tigstart, tigend, strand = line.split()
            agpblock[chrid].append([chrnumber, seqstart, seqend, num, attr, tigid, tigstart, tigend, strand])
    
    # Parse centromere file if provided
    centromeres = {}
    if centrofile is not None:
        try:
            centromeres = parse_centromere_file(centrofile)
        except Exception as e:
            print(f'[Warning] Error parsing centromere file: {e}')
            centromeres = {}
    
    # Write chromosome information file
    with open(f'tmp/{outprefix}.chr.txt', 'w') as chrfile:
        chrfile.write('Chr\tStart\tEnd\n')
        for chrid in agpblock:
            chrfile.write(f'{agpblock[chrid][-1][0]}\t1\t{agpblock[chrid][-1][2]}\n')
    
    # Prepare data for visualization
    has_heatmap = False
    has_markers = False
    
    # Write heatmap data for centromeres (to actually color the full chromosome regions)
    if centromeres:
        with open(f'tmp/{outprefix}.heatmap.txt', 'w') as heatmapfile:
            heatmapfile.write('Chr\tStart\tEnd\tValue\n')
            for chrid, (start, end, color) in centromeres.items():
                if chrid in agpblock:
                    chr_num = agpblock[chrid][-1][0]
                    # Use discrete value for each centromere (avoid gradient)
                    heatmapfile.write(f'{chr_num}\t{start}\t{end}\t1\n')
                    has_heatmap = True
                    print(f'[Info] Added centromere heatmap region for {chrid} at {start}-{end} (color: #{color})')
                else:
                    print(f'[Warning] Centromere found for {chrid} but chromosome not in AGP file')
        
        # Write color mapping file for discrete colors
        if has_heatmap:
            unique_colors = list(set(color for _, _, color in centromeres.values()))
            with open(f'tmp/{outprefix}.color_config.txt', 'w') as colorfile:
                # Create color mapping - RIdeogram can use custom colors
                for i, color in enumerate(unique_colors):
                    colorfile.write(f'#{color}\n')
                print(f'[Info] Created color mapping for {len(unique_colors)} unique centromere colors')
    
    # Write marker data for telomeres and other features
    marker_data = []
    
    # Add 5 Mb position markers
    if show_position_markers:
        for chrid in agpblock:
            chr_num = agpblock[chrid][-1][0]
            chr_size = int(agpblock[chrid][-1][2])
            
            # 5 Mb intervals
            interval = 5000000
            positions = list(range(interval, chr_size + 1, interval))
            
            for pos in positions:
                if pos <= chr_size:
                    marker_data.append(f'5Mb\tcircle\t{chr_num}\t{pos}\t{pos}\t666666\n')
                    has_markers = True
    
    # Add telomeres as markers
    if telofile is not None:
        with open(telofile, 'r') as te:
            for line in te:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 7:
                    chrid, chrlen, status, leftnum, leftdirection, rightnum, rightdirection = parts[:7]
                    
                    if chrid not in agpblock:
                        continue
                    
                    chr_num = agpblock[chrid][-1][0]
                    
                    # Left telomere
                    if int(leftnum) != 0:
                        marker_data.append(f'telomere\ttriangle\t{chr_num}\t1\t10000\t0000ff\n')
                        has_markers = True
                    
                    # Right telomere  
                    if int(rightnum) != 0:
                        marker_data.append(f'telomere\ttriangle\t{chr_num}\t{int(chrlen)-10000}\t{chrlen}\t0000ff\n')
                        has_markers = True
    
    # Add gaps as markers if requested
    if show_gaps:
        for chrid in agpblock:
            for line in agpblock[chrid]:
                if line[4] != 'W':  # Non-sequence (gap)
                    marker_data.append(f'gap\tbox\t{line[0]}\t{line[1]}\t{line[2]}\tff7f00\n')
                    has_markers = True
    
    # Write marker file if needed
    if has_markers:
        with open(f'tmp/{outprefix}.marker.txt', 'w') as markerfile:
            markerfile.write('Type\tShape\tChr\tStart\tEnd\tcolor\n')
            markerfile.writelines(marker_data)
    
    # Generate R script with proper region coloring
    with open(f'tmp/{outprefix}.genomedrawer.r', 'w') as r:
        if has_heatmap and has_markers:
            # Both heatmap regions and markers
            Rscript = f'''library(RIdeogram)
chr <- read.table("tmp/{outprefix}.chr.txt", sep = "\\t", header = T, stringsAsFactors = F)
heatmap_data <- read.table("tmp/{outprefix}.heatmap.txt", sep = "\\t", header = T, stringsAsFactors = F)
marker <- read.table("tmp/{outprefix}.marker.txt", sep = "\\t", header = T, stringsAsFactors = F)

print("Chromosome data:")
print(chr)
print("Heatmap data for centromeres:")
print(heatmap_data)
print("Marker data:")
print(head(marker))

# Create custom color vector for centromeres'''

            # Add custom colors for each chromosome's centromere
            centromere_colors = []
            for chrid, (start, end, color) in centromeres.items():
                if chrid in agpblock:
                    centromere_colors.append(f'#{color}')
            
            if centromere_colors:
                Rscript += f'''
custom_colors <- c({", ".join([f'"{c}"' for c in centromere_colors])})
print("Custom colors:")
print(custom_colors)

ideogram(karyotype = chr, overlaid = heatmap_data, label = marker, label_type = "marker", 
         colorset1 = custom_colors, output = "{outprefix}.svg")'''
            else:
                Rscript += f'''
ideogram(karyotype = chr, overlaid = heatmap_data, label = marker, label_type = "marker", output = "{outprefix}.svg")'''
            
            Rscript += f'''
convertSVG("{outprefix}.svg", file = "{outprefix}", device = "png")'''

        elif has_heatmap:
            # Only heatmap regions
            Rscript = f'''library(RIdeogram)
chr <- read.table("tmp/{outprefix}.chr.txt", sep = "\\t", header = T, stringsAsFactors = F)
heatmap_data <- read.table("tmp/{outprefix}.heatmap.txt", sep = "\\t", header = T, stringsAsFactors = F)

print("Chromosome data:")
print(chr)
print("Heatmap data for centromeres:")
print(heatmap_data)

# Create custom color vector'''
            
            # Add custom colors
            centromere_colors = []
            for chrid, (start, end, color) in centromeres.items():
                if chrid in agpblock:
                    centromere_colors.append(f'#{color}')
            
            if centromere_colors:
                Rscript += f'''
custom_colors <- c({", ".join([f'"{c}"' for c in centromere_colors])})
print("Custom colors:")
print(custom_colors)

ideogram(karyotype = chr, overlaid = heatmap_data, colorset1 = custom_colors, output = "{outprefix}.svg")'''
            else:
                Rscript += f'''
ideogram(karyotype = chr, overlaid = heatmap_data, output = "{outprefix}.svg")'''
            
            Rscript += f'''
convertSVG("{outprefix}.svg", file = "{outprefix}", device = "png")'''
            
        elif has_markers:
            # Only markers
            Rscript = f'''library(RIdeogram)
chr <- read.table("tmp/{outprefix}.chr.txt", sep = "\\t", header = T, stringsAsFactors = F)
marker <- read.table("tmp/{outprefix}.marker.txt", sep = "\\t", header = T, stringsAsFactors = F)
ideogram(karyotype = chr, label = marker, label_type = "marker", output = "{outprefix}.svg")
convertSVG("{outprefix}.svg", file = "{outprefix}", device = "png")'''
        else:
            # Basic ideogram only
            Rscript = f'''library(RIdeogram)
chr <- read.table("tmp/{outprefix}.chr.txt", sep = "\\t", header = T, stringsAsFactors = F)
ideogram(karyotype = chr, output = "{outprefix}.svg")
convertSVG("{outprefix}.svg", file = "{outprefix}", device = "png")'''
        r.write(Rscript)
    
    # Run R script
    cmdr = subprocess.run(f'Rscript tmp/{outprefix}.genomedrawer.r', 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    if cmdr.returncode != 0:
        print(f'[Warning] Unexpected error occurred in Rscript figure drawing:')
        print(f'cmd: {cmdr.args}')
        print(f'returncode: {cmdr.returncode}')
        print('stdout:')
        print(cmdr.stdout.decode("utf-8"))
        print('stderr:')
        print(cmdr.stderr.decode("utf-8"))
        return False
    else:
        print(f'[Output] Chromosome plot written to: {outprefix}.png')
        return True


def main():
    parser = argparse.ArgumentParser(
        description='Draw chromosome ideograms with telomere and colored centromere information',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage with FASTA file
  python3 telomere_drawer.py -f genome.fasta -o output

  # With telomere information
  python3 telomere_drawer.py -f genome.fasta -t telomere.info -o output

  # With centromere and telomere information
  python3 telomere_drawer.py -f genome.fasta -c centromere.txt -t telomere.info -o output

  # With colored centromeres
  python3 telomere_drawer.py -f genome.fasta -c colored_centromeres.txt -o output

  # Show gaps in the ideogram
  python3 telomere_drawer.py -f genome.fasta -t telomere.info -o output --show-gaps

  # Hide 5 Mb position markers
  python3 telomere_drawer.py -f genome.fasta -t telomere.info -o output --hide-position-markers

Centromere file format:
  chromosome,start-end[,color]
  
  Without colors (defaults to red):
    scaffold1,37651276-38156789
    scaffold2,32935702-33158580
  
  With custom colors (hex without #):
    scaffold1,37651276-38156789,ff0000
    scaffold2,32935702-33158580,00ff00
    scaffold3,29623332-29948992,0000ff
    scaffold4,41961254-42174310,ff6600
        ''')
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-f', '--fasta', 
                           help='Input FASTA file (will generate AGP automatically)')
    input_group.add_argument('-a', '--agp', 
                           help='Input AGP file')
    
    # Optional files
    parser.add_argument('-c', '--centromere', 
                       help='Centromere file in format: chromosome,start-end[,color] (color as hex without #, optional)')
    parser.add_argument('-t', '--telomere', 
                       help='Telomere information file (quarTeT .telo.info format)')
    parser.add_argument('--tidk-telomere', 
                       help='Telomere file in tidk format (telomeric_repeat_windows.tsv)')
    parser.add_argument('--min-telomere-repeats', type=int, default=100,
                       help='Minimum repeat count to consider as telomere for tidk format (default: 100, same as quarTeT)')
    
    # Output and display options
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for generated files')
    parser.add_argument('--show-gaps', action='store_true',
                       help='Show gaps in the ideogram (default: False)')
    parser.add_argument('--hide-position-markers', action='store_true',
                       help='Hide 5 Mb position markers (default: False)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Handle display logic
    show_position_markers = not args.hide_position_markers
    
    # Determine AGP file
    if args.fasta:
        if not os.path.exists(args.fasta):
            print(f'[Error] FASTA file not found: {args.fasta}')
            sys.exit(1)
        
        agp_file = f'tmp/{args.output}.genome.agp'
        
        if args.show_gaps:
            print(f'[Info] Generating AGP file with gap information from FASTA: {args.fasta}')
            agpGap(args.fasta, agp_file)
        else:
            print(f'[Info] Generating simple AGP file (no gap detection) from FASTA: {args.fasta}')
            create_simple_agp_from_fasta(args.fasta, agp_file)
    else:
        if not os.path.exists(args.agp):
            print(f'[Error] AGP file not found: {args.agp}')
            sys.exit(1)
        agp_file = args.agp
    
    # Check optional files
    if args.centromere and not os.path.exists(args.centromere):
        print(f'[Warning] Centromere file not found: {args.centromere}')
        args.centromere = None
    
    # Handle telomere files
    telomere_file = None
    if args.telomere and args.tidk_telomere:
        print('[Error] Cannot specify both --telomere and --tidk-telomere')
        sys.exit(1)
    
    if args.telomere:
        if not os.path.exists(args.telomere):
            print(f'[Warning] Telomere file not found: {args.telomere}')
        else:
            # Check if it's tidk format by looking at the first line
            with open(args.telomere, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('id\t') or 'forward_repeat_number' in first_line:
                    print('[Info] Detected tidk format in telomere file, converting...')
                    # Convert tidk format to telo.info format
                    converted_file = f'tmp/{args.output}.telo.info'
                    fasta_file = args.fasta if args.fasta else None
                    if fasta_file is None:
                        print('[Error] FASTA file required when using tidk format for chromosome length information')
                        sys.exit(1)
                    telomere_file = convert_tidk_to_telo_info(args.telomere, fasta_file, converted_file, args.min_telomere_repeats)
                else:
                    telomere_file = args.telomere
    
    if args.tidk_telomere:
        if not os.path.exists(args.tidk_telomere):
            print(f'[Warning] Tidk telomere file not found: {args.tidk_telomere}')
        else:
            # Convert tidk format to telo.info format
            converted_file = f'tmp/{args.output}.telo.info'
            fasta_file = args.fasta if args.fasta else None
            if fasta_file is None:
                print('[Error] FASTA file required when using --tidk-telomere for chromosome length information')
                sys.exit(1)
            telomere_file = convert_tidk_to_telo_info(args.tidk_telomere, fasta_file, converted_file, args.min_telomere_repeats)
    
    # Draw ideogram
    print(f'[Info] Drawing chromosome ideogram...')
    print(f'[Info] Parameters: AGP={agp_file}, Centromere={args.centromere}, Telomere={telomere_file}, Show gaps={args.show_gaps}, Show position markers={show_position_markers}')
    
    success = drawTelomereIdeogram(
        agpfile=agp_file,
        outprefix=args.output,
        centrofile=args.centromere,
        telofile=telomere_file,
        show_gaps=args.show_gaps,
        show_position_markers=show_position_markers
    )
    
    if success:
        print('[Info] Ideogram drawing completed successfully!')
    else:
        print('[Error] Failed to generate ideogram')
        sys.exit(1)


if __name__ == '__main__':
    main()