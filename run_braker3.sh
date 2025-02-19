#!/bin/bash

# Save the original command line for logging
original_command_line="$0 $*"

usage() {
    echo "Usage: $0 -t threads -g genome -s species -w working_dir -p protein -b bam_dir -l busco_lineage"
    exit 1
}

# Parse command-line options
while getopts ":t:g:s:w:p:b:l:" opt; do
  case $opt in
    t) threads="$OPTARG" ;;
    g) genome="$OPTARG" ;;
    s) species="$OPTARG" ;;
    w) wd="$OPTARG" ;;
    p) protein="$OPTARG" ;;
    b) bam_dir="$OPTARG" ;;
    l) busco_lineage="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Ensure that all required parameters have been provided
if [[ -z "$threads" || -z "$genome" || -z "$species" || -z "$wd" || -z "$protein" || -z "$bam_dir" || -z "$busco_lineage" ]]; then
    echo "Error: Missing required parameter(s)." >&2
    usage
fi

# Create the working directory if it doesn't exist
mkdir -p "$wd"

# Print the command line used (this output will be captured if you redirect stdout)
echo "Command line: $original_command_line"

# Determine if bam_dir is a file (single BAM) or a directory (multiple BAMs)
if [[ -f "$bam_dir" ]]; then
    # bam_dir is a single BAM file
    bam_files="$bam_dir"
elif [[ -d "$bam_dir" ]]; then
    # bam_dir is a directory; gather all BAM files in that directory
    bam_files=$(find "$bam_dir" -maxdepth 1 -type f -name "*.bam")
else
    echo "Error: $bam_dir is neither a file nor a directory." >&2
    exit 1
fi

# Check that we found at least one BAM file
if [[ -z "$bam_files" ]]; then
    echo "Error: No BAM files found in $bam_dir" >&2
    exit 1
fi

# (Optional) Write the list of BAM files to a text file
bam_list_file="$wd/bam_list.txt"
echo "$bam_files" | tr ' ' '\n' > "$bam_list_file"

# Run braker.pl with the gathered BAM file(s)
braker.pl --genome="$genome" \
          --species="$species" \
          --prot_seq="$protein" \
          --bam="$bam_files" \
          --threads="$threads" \
          --workingdir="$wd" \
          --busco_lineage="$busco_lineage"

echo "Finished."
