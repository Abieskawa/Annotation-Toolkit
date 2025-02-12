#!/bin/bash
#----------------------------------------------------------------------
# Minimal TRF Masking Pipeline with Automatic Script Detection
# and a symlink approach for a clean base name.
#
# This script uses:
#   - splitMfasta.pl (from https://github.com/Gaius-Augustus/Augustus)
#   - trf
#   - parseTrfOutput.py (from https://github.com/gatech-genemark/BRAKER2-exp)
#   - bedtools
#   - GNU parallel
#
# It expects an input genome file (which may include extensions such as .fa, .fasta, .masked, etc.).
# The script creates a symlink with a clean base name so that output files are named as:
#   <BASE_NAME>.split.<N>.fa, <BASE_NAME>.combined.masked, etc.
#
# Usage:
#   ./trf_mask.sh -i <path/to/input> -o <output_dir> [-t <threads>]
#
# Example:
#   nohup bash -c "time (bash trf_mask.sh -i /home/abieskawa/output/asian_hard_clam/03_tetools_SoftMask/MTW_genome.fasta.masked \
#       -o /home/abieskawa/output/asian_hard_clam/04_TRF_SoftMask -t 128 > repeatmask_trf.log 2>&1)" \
#       2>run_trf_time.log >/dev/null &
#----------------------------------------------------------------------

usage() {
    echo "Usage: $0 -i <input_file> -o <output_dir> [-t <threads>]"
    echo ""
    echo "  -i, --input     Path to the input genome file (can include extensions)"
    echo "  -o, --output    Output directory (all results will be written here)"
    echo "  -t, --threads   (Optional) Number of threads for GNU parallel"
    echo "  -h, --help      Show this help message"
    exit 1
}

#----------------------------------------------------------------------
# Parse command-line arguments
#----------------------------------------------------------------------
if [ "$#" -lt 4 ]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        -o|--output)
            OUTDIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

if [ -z "$INPUT" ] || [ -z "$OUTDIR" ]; then
    usage
fi

#----------------------------------------------------------------------
# Compute the absolute path of the input file BEFORE changing directories.
#----------------------------------------------------------------------
ABS_INPUT=$(realpath "$INPUT") || { echo "Error: Unable to resolve absolute path for $INPUT"; exit 1; }

#--------------------------------------------------------------------------------
# Create a clean base name by removing common extensions.
# For example: "MTW_genome.fasta.masked" becomes "MTW_genome"
#--------------------------------------------------------------------------------
BASE=$(basename "$INPUT")                    
BASE_NAME=$(echo "$BASE" | sed -E 's/\.(fa|fna|fasta|masked)//g')  

#--------------------------------------------------------------------------------
# Automatically detect the directory of this script.
# Place 'splitMfasta.pl' and 'parseTrfOutput.py' in the same directory.
#--------------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SPLIT_MFASTA="${SCRIPT_DIR}/splitMfasta.pl"
PARSE_TRF_OUTPUT="${SCRIPT_DIR}/parseTrfOutput.py"

#--------------------------------------------------------------------------------
# Check that required scripts and tools are available.
#--------------------------------------------------------------------------------
for f in "$SPLIT_MFASTA" "$PARSE_TRF_OUTPUT"; do
    [ -x "$f" ] || { echo "Error: $f not found or not executable"; exit 1; }
done

command -v trf >/dev/null     || { echo "Error: trf not in PATH"; exit 1; }
command -v bedtools >/dev/null|| { echo "Error: bedtools not in PATH"; exit 1; }
command -v parallel >/dev/null|| { echo "Error: GNU parallel not found"; exit 1; }

#--------------------------------------------------------------------------------
# Print tool checks (optional)
#--------------------------------------------------------------------------------
echo "===== Tool Check ====="
echo "splitMfasta.pl: $("$SPLIT_MFASTA" 2>&1 | head -n 5 | grep -m1 -i usage || echo 'OK')"
echo "TRF: $(trf -v | head -n1)"
if "$PARSE_TRF_OUTPUT" 2>&1 | grep -q "usage"; then
    echo "parseTrfOutput.py: Available"
else
    echo "parseTrfOutput.py: Unknown response, but file is executable"
fi
echo "bedtools: $(bedtools --version)"
echo "parallel: $(parallel --version 2>&1 | head -n1)"
echo "========================"

#--------------------------------------------------------------------------------
# Prepare output directory and change into it.
#--------------------------------------------------------------------------------
mkdir -p "$OUTDIR" || { echo "Cannot create output directory $OUTDIR"; exit 1; }
cd "$OUTDIR" || { echo "Cannot change to output directory $OUTDIR"; exit 1; }

#--------------------------------------------------------------------------------
# Create a symlink in OUTDIR with a clean base name.
# The symlink will be named "<BASE_NAME>.fa" and point to the real input file.
#--------------------------------------------------------------------------------
CLEAN_FASTA="${BASE_NAME}.fa"
if [ ! -f "$CLEAN_FASTA" ]; then
    ln -sf "$ABS_INPUT" "$CLEAN_FASTA"
fi

echo "Using symlink: $(ls -l $CLEAN_FASTA)"

#--------------------------------------------------------------------------------
# Main Pipeline Steps
#--------------------------------------------------------------------------------

echo "Splitting $CLEAN_FASTA with splitMfasta.pl..."
# Using --outputpath=./ ensures that the split files are created in the current directory.
"$SPLIT_MFASTA" --minsize=25000000 --outputpath=./ "$CLEAN_FASTA"

echo "Running TRF on each split FASTA file..."
ls ${BASE_NAME}.split.*.fa | parallel ${THREADS:+-j $THREADS} 'trf {} 2 7 7 80 10 50 500 -d -m -h'

echo "Parsing TRF output..."
ls ${BASE_NAME}.split.*.fa.2.7.7.80.10.50.500.dat | \
    parallel ${THREADS:+-j $THREADS} "${PARSE_TRF_OUTPUT} {} --minCopies 1 --statistics {}.STATS > {}.raw.gff 2> {}.parsedLog"

echo "Sorting parsed output..."
ls ${BASE_NAME}.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff | \
    parallel ${THREADS:+-j $THREADS} 'sort -k1,1 -k4,4n -k5,5n {} > {}.sorted 2> {}.sortLog'

echo "Merging GFF files..."
FILES=${BASE_NAME}.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff.sorted
for f in $FILES; do
    bedtools merge -i "$f" | \
      awk 'BEGIN{OFS="\t"} {print $1,"trf","repeat",$2+1,$3,".",".",".","."}' > "$f.merged.gff" 2> "$f.bedtools_merge.log"
done

echo "Masking FASTA chunks..."
ls ${BASE_NAME}.split.*.fa | \
    parallel ${THREADS:+-j $THREADS} 'bedtools maskfasta -fi {} -bed {}.2.7.7.80.10.50.500.dat.raw.gff.sorted.merged.gff -fo {}.combined.masked -soft &> {}.bedools_mask.log'

echo "Concatenating masked FASTA chunks..."
cat ${BASE_NAME}.split.*.fa.combined.masked > ${BASE_NAME}.combined.masked

echo "Pipeline completed at $(date)"
