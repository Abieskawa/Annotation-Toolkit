#!/bin/bash
# run_cutadapt.sh
#
# This script uses cutadapt to process paired‐end FASTQ files.
#
# It supports front‐trimming for each sample via the -u (read1) and -U (read2)
# options. The front‐trim values can be supplied either as a single number (applied
# globally) or via a file (specified with -f) that contains either 2 or 3 columns:
#
#   2 columns: <sample_filename_prefix> <trimValue>
#       (Both read1 and read2 will be trimmed by that many bases.)
#
#   3 columns: <sample_filename_prefix> <trimValue_R1> <trimValue_R2>
#       (For example, if your FASTQ files are named SRR26245751_1.fastq and
#        SRR26245751_2.fastq then the file should contain keys “SRR26245751_1” and “SRR26245751_2”.)
#
# Adapter sequences can be specified via -a1 and -a2. If not provided, the script will
# look for a file named "adapter_list.txt" in the same directory as this script; if that
# file isn’t found, default adapters are used:
#
#   For read1 (-a):
#     AGATCGGAAGAG
#     AAAAAAAAAAAA
#     GGGGGGGGGGGG
#
#   For read2 (-A):
#     AGATCGGAAGAG
#     AAAAAAAAAAAA
#     GGGGGGGGGGGG
#
# The option --detect_adapter_for_pe is parsed (and logged) but not used by cutadapt.
#
# Example (global front trimming):
#   nohup bash run_cutadapt.sh -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA \\
#         -t 128 -l 5 -q 20 -Q 20 -f 10 -a1 YOUR_ADAPTER1 -a2 YOUR_ADAPTER2 > run_cutadapt.log 2>&1 &
#
# Example (per-sample trimming using a trim file):
#   nohup bash run_cutadapt.sh -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA \\
#         -t 128 -l 5 -q 20 -Q 20 -f sample_trim.txt > run_cutadapt.log 2>&1 &

############################################
# Defaults for optional parameters
############################################
threads=16
min_length=5         # Minimum read length (-m option)
quality1=20          # Quality cutoff for read1 (-q option)
quality2=20          # Quality cutoff for read2 (-Q option)
trim_front=""        # Either a single number (global) or a filename with per-sample trim values

# Adapter parameters. If provided, these override the defaults.
adapter1=""
adapter2=""

# (Note: --detect_adapter_for_pe is parsed but not used by cutadapt.)
detect_adapter_for_pe=0

# Required parameters
wd=""
output_dir=""

############################################
# Function to display usage information
############################################
usage() {
  cat <<EOF
Usage: $0 -d <working_directory> -o <output_directory> -t <threads> -l <min_length> \\
          -q <quality_R1> -Q <quality_R2> [-f <trim_front_value_or_file>] \\
          [-a1 <adapter_read1>] [-a2 <adapter_read2>] [--detect_adapter_for_pe]

Example (global front trimming):
  $0 -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA -t 16 -l 5 -q 20 -Q 20 -f 10 -a1 YOUR_ADAPTER1 -a2 YOUR_ADAPTER2

Example (per-sample trimming):
  $0 -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA -t 16 -l 5 -q 20 -Q 20 -f sample_trim.txt

If you supply a trim file to -f, it can have either 2 or 3 columns per line:
  2 columns: <sample_filename_prefix> <trimValue>
    => The same trim value is used for both read1 and read2.
  3 columns: <sample_filename_prefix> <trimValue_R1> <trimValue_R2>
    => For example, if your FASTQ files are named SRR26245751_1.fastq and
       SRR26245751_2.fastq, then the file should contain keys “SRR26245751_1” and “SRR26245751_2”.

Options:
  -d      Working directory (required)
  -o      Output directory (required; must already exist)
  -t      Number of threads (default: $threads)
  -l      Minimum read length (default: $min_length)
  -q      Quality cutoff for read1 (default: $quality1)
  -Q      Quality cutoff for read2 (default: $quality2)
  -f      (Optional) Either a single number or a filename with sample-specific trimming values
  -a1     (Optional) Adapter sequence for read1
  -a2     (Optional) Adapter sequence for read2
  --detect_adapter_for_pe  (Parsed but not used)
  -? or --help   Display this help message
EOF
  exit 1
}

############################################
# If no arguments provided, show usage
############################################
if [ $# -eq 0 ]; then
  usage
fi

############################################
# Save the original command line for logging
############################################
orig_cmd="$0 $@"
echo "Script invoked with: $orig_cmd"

############################################
# Parse command-line options
############################################
while [ $# -gt 0 ]; do
  case "$1" in
    -d)
      wd="$2"
      shift 2
      ;;
    -o)
      output_dir="$2"
      shift 2
      ;;
    -t)
      threads="$2"
      shift 2
      ;;
    -l)
      min_length="$2"
      shift 2
      ;;
    -q)
      quality1="$2"
      shift 2
      ;;
    -Q)
      quality2="$2"
      shift 2
      ;;
    -f)
      trim_front="$2"
      shift 2
      ;;
    -a1)
      adapter1="$2"
      shift 2
      ;;
    -a2)
      adapter2="$2"
      shift 2
      ;;
    --detect_adapter_for_pe)
      detect_adapter_for_pe=1
      shift
      ;;
    -\? | --help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

############################################
# Validate required parameters
############################################
if [[ -z "$wd" ]]; then
  echo "Error: Working directory (-d) not specified."
  usage
fi
if [[ -z "$output_dir" ]]; then
  echo "Error: Output directory (-o) not specified."
  usage
fi

############################################
# Convert directories to absolute paths
############################################
initial_dir=$(pwd)
if [[ ! -d "$wd" ]]; then
  echo "Error: Working directory '$wd' does not exist!"
  exit 1
fi
wd=$(cd "$wd" && pwd)

# Make output_dir absolute if it isn't already
if [[ "${output_dir:0:1}" != "/" ]]; then
  output_dir="$initial_dir/$output_dir"
fi

############################################
# Process the -f parameter (front trim)
############################################
# If -f is provided, it can be a file (with 2 or 3 columns) or a global number.
declare global_trim_R1=""
declare global_trim_R2=""
declare -A per_sample_trim_R1
declare -A per_sample_trim_R2

if [[ -n "$trim_front" ]]; then
  if [[ -f "$trim_front" ]]; then
    # Parse the file line by line.
    while IFS= read -r line || [[ -n "$line" ]]; do
      # Remove carriage return if present (handles DOS/Windows line endings)
      line=$(echo "$line" | tr -d '\r')
      [[ -z "$line" || "$line" =~ ^# ]] && continue
      fields=( $line )
      sample="${fields[0]}"
      if [[ ${#fields[@]} -eq 2 ]]; then
        per_sample_trim_R1["$sample"]="${fields[1]}"
        per_sample_trim_R2["$sample"]="${fields[1]}"
      elif [[ ${#fields[@]} -eq 3 ]]; then
        per_sample_trim_R1["$sample"]="${fields[1]}"
        per_sample_trim_R2["$sample"]="${fields[2]}"
      else
        echo "Warning: Unexpected number of columns in trim file line: $line"
      fi
    done < "$trim_front"
  else
    # Otherwise, treat -f as a global trim value.
    global_trim_R1="$trim_front"
    global_trim_R2="$trim_front"
  fi
fi

############################################
# Change to the working directory
############################################
cd "$wd" || { echo "Error: Failed to change directory to '$wd'"; exit 1; }

############################################
# Find unique sample basenames from FASTQ files
############################################
# This script expects paired-end files named as: <basename>_1.fastq and <basename>_2.fastq.
# We extract the basename as the part before the first underscore.
uniq_basenames=$(find . -maxdepth 1 -type f -name "*.fastq" \
  | sed 's|^\./||' \
  | cut -d '_' -f1 \
  | sort | uniq)

############################################
# Process each sample
############################################
for basename in $uniq_basenames; do
  read1="${basename}_1.fastq"
  read2="${basename}_2.fastq"
  
  if [[ ! -f "$read1" || ! -f "$read2" ]]; then
    echo "Skipping sample '$basename': one or both paired-end files are missing."
    continue
  fi
  
  echo "Processing sample: $basename"
  
  ############################################
  # Determine front trim values for each read.
  # (We append _1 and _2 to the basename to look up keys in the trim file.)
  ############################################
  sample_key_r1="${basename}_1"
  sample_key_r2="${basename}_2"
  sample_trim_r1=""
  sample_trim_r2=""
  if [[ -n "${per_sample_trim_R1[$sample_key_r1]}" ]]; then
    sample_trim_r1="${per_sample_trim_R1[$sample_key_r1]}"
  elif [[ -n "$global_trim_R1" ]]; then
    sample_trim_r1="$global_trim_R1"
  fi
  if [[ -n "${per_sample_trim_R2[$sample_key_r2]}" ]]; then
    sample_trim_r2="${per_sample_trim_R2[$sample_key_r2]}"
  elif [[ -n "$global_trim_R2" ]]; then
    sample_trim_r2="$global_trim_R2"
  fi
  
  ############################################
  # Define adapter options for cutadapt
  ############################################
  # Determine the directory where this script resides.
  # Using realpath here so that if the script is called via a relative path,
  # we still get the absolute directory.
  script_dir=$(cd "$(dirname "$(realpath "$0")")" && pwd)
  adapter_fasta="$script_dir/adapter_list.txt"
  
  if [[ -n "$adapter1" ]]; then
    adapter_opts1=( -a "$adapter1" )
  elif [[ -f "$adapter_fasta" ]]; then
    adapter_opts1=( -a "file:$adapter_fasta" )
  else
    adapter_opts1=( -a "AGATCGGAAGAG" -a "AAAAAAAAAAAA" -a "GGGGGGGGGGGG" )
  fi
  
  if [[ -n "$adapter2" ]]; then
    adapter_opts2=( -A "$adapter2" )
  elif [[ -f "$adapter_fasta" ]]; then
    adapter_opts2=( -A "file:$adapter_fasta" )
  else
    adapter_opts2=( -A "AGATCGGAAGAG" -A "AAAAAAAAAAAA" -A "GGGGGGGGGGGG" )
  fi
  
  ############################################
  # Build the cutadapt command
  ############################################
  cmd=(cutadapt)
  cmd+=( "${adapter_opts1[@]}" "${adapter_opts2[@]}" )
  cmd+=( -j "$threads" -q "$quality1" -Q "$quality2" -m "$min_length" )
  
  # Add front trim options if trim values were determined.
  if [[ -n "$sample_trim_r1" ]]; then
    cmd+=( -u "$sample_trim_r1" )
  fi
  if [[ -n "$sample_trim_r2" ]]; then
    cmd+=( -U "$sample_trim_r2" )
  fi
  
  # Define output filenames (placed in the output directory)
  out1="${output_dir}/${basename}_1.cleaned.fastq"
  out2="${output_dir}/${basename}_2.cleaned.fastq"
  
  cmd+=( -o "$out1" -p "$out2" "$read1" "$read2" )
  
  echo "Running command: ${cmd[*]}"
  "${cmd[@]}"
done

echo "All samples processed at $(date)"
