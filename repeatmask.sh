#!/bin/bash

# Short help function
show_help() {
  echo "Usage: $0 [options]
Options:
  -i, --input         Input genome file
  -o, --output        Output directory
  -t, --threads       Number of threads
  -l, --library       Path to repeat library (FamDB dir or .h5)
  -s, --search        Search target
  -f, --file-tag      File tag for output
      --log           Log file (default: repeatmask.log)
  -h, --help          Show this help and exit

Example:
  $0 -i Tilapia_TFS_strain_genome_v1.fasta \\
     -o /output/marine_tilapia \\
     -t 128 \\
     -l /output/Lib_fish/famdb \\
     -s \"oreochromis niloticus\" \\
     -f \"marine_tilapia\" \\
     --log mylog.log
"
  exit
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)         original_genome="$2"; shift ;;
    -o|--output)        output_dir="$2"; shift ;;
    -t|--threads)       threads="$2"; shift ;;
    -l|--library)       configured_library="$2"; shift ;;
    -s|--search)        search_target="$2"; shift ;;
    -f|--file-tag)      file_tag="$2"; shift ;;
    --log)              log_file="$2"; shift ;;
    -h|--help)          show_help ;;
    *) echo "Unknown parameter: $1"; show_help ;;
  esac
  shift
done

# Required args
if [[ -z "$original_genome" || -z "$output_dir" || -z "$threads" || -z "$configured_library" || -z "$search_target" || -z "$file_tag" ]]; then
  echo "Error: Missing required arguments."
  show_help
fi

# Resolve input to an absolute path (so later cd's don't break it)
case "$original_genome" in
  /*) original_genome_abs="$original_genome" ;;
   *) original_genome_abs="$(pwd)/$original_genome" ;;
esac
[[ -r "$original_genome_abs" ]] || { echo "Input not readable: $original_genome_abs" >&2; exit 1; }

# Defaults/paths
log_file="${log_file:-repeatmask.log}"
mkdir -p "$output_dir"

genome_base="$(basename "$original_genome_abs")"
softmask_dir="$output_dir/03_SoftMask"

# Idempotency guard
if [[ -d "$softmask_dir" ]] && compgen -G "$softmask_dir/${genome_base}.*" > /dev/null; then
  echo "Error: Genome '$genome_base' appears to have been processed already in '$softmask_dir'."
  echo "Exiting to avoid overwriting records."
  exit 1
fi

# Log path absolute
if [[ "$log_file" != /* ]]; then
  log_file="$output_dir/$log_file"
fi

# Start logging
{
  echo "Command: $0 $*"
  echo "Working dir: $(pwd); Host: $(hostname); Date: $(date)"
  echo "Log file path: $log_file"
  echo "----------------------------------------"
} >> "$log_file"

# Create subdirs
mkdir -p "$output_dir"/{01_RepeatMasker,01_RepeatModeler,02_RepeatMasker,03_SoftMask}

##################################
# 1) Prepare additional RepeatMasker library
##################################
cd "$output_dir/01_RepeatMasker" || exit
echo "Step 1: Creating additional RepeatMasker library in $(pwd)" >> "$log_file"

famdb.py -i "$configured_library" families -f embl -ad "$search_target" \
  > "${file_tag}_ad.embl" 2>>"$log_file"

/opt/RepeatMasker/util/buildRMLibFromEMBL.pl "${file_tag}_ad.embl" \
  > "${file_tag}_ad.fa" 2>>"$log_file"

##################################
# 2) RepeatModeler
##################################
cd "$output_dir/01_RepeatModeler" || exit
echo "Step 2: Running RepeatModeler in $(pwd)" >> "$log_file"

BuildDatabase -name GDB -engine ncbi "$original_genome_abs" >>"$log_file" 2>&1
RepeatModeler -engine ncbi -threads "$threads" -database GDB -LTRStruct >>"$log_file" 2>&1

##################################
# 3) Final RepeatMasker run
##################################
cd "$output_dir/02_RepeatMasker" || exit
echo "Step 3: Final RepeatMasker run in $(pwd)" >> "$log_file"

cat ../01_RepeatModeler/GDB-families.fa ../01_RepeatMasker/"${file_tag}_ad.fa" > repeat_db.fa
# -dir . forces RepeatMasker outputs into the current directory.
RepeatMasker -xsmall -gff -html -lib repeat_db.fa -pa "$threads" -dir . "$original_genome_abs" >>"$log_file" 2>&1

##################################
# 4) Move soft-masked results
##################################
echo "Step 4: Moving soft-masked files for $genome_base to $softmask_dir" >> "$log_file"
mkdir -p "$softmask_dir"

for ext in cat.gz masked out out.gff out.html tbl; do
  f="${genome_base}.${ext}"
  [[ -e "$f" ]] && mv "$f" "$softmask_dir"/ >>"$log_file" 2>&1
done

echo "All steps completed. See $log_file for details." >> "$log_file"
