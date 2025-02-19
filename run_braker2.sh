#!/bin/bash

# Save the original command line for logging
original_command_line="$0 $*"

usage() {
    echo "Usage: $0 -t threads -g genome -s species -w working_dir -p protein -l busco_lineage"
    exit 1
}

# Parse command-line options
while getopts ":t:g:s:w:p:l:" opt; do
  case $opt in
    t) threads="$OPTARG" ;;
    g) genome="$OPTARG" ;;
    s) species="$OPTARG" ;;
    w) wd="$OPTARG" ;;
    p) protein="$OPTARG" ;;
    l) busco_lineage="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Ensure that all required parameters have been provided
if [[ -z "$threads" || -z "$genome" || -z "$species" || -z "$wd" || -z "$protein" || -z "$busco_lineage" ]]; then
    echo "Error: Missing required parameter(s)." >&2
    usage
fi

# Create the working directory if it doesn't exist
mkdir -p "$wd"

# Print the command line used (this output will be captured if you redirect stdout)
echo "Command line: $original_command_line"

# Run braker.pl without BAM evidence
braker.pl --genome="$genome" \
          --species="$species" \
          --prot_seq="$protein" \
          --threads="$threads" \
          --workingdir="$wd" \
          --busco_lineage="$busco_lineage"

echo "Finished."
