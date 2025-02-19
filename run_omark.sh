#!/bin/bash
usage(){ echo "Usage: $0 -i <input_file> -base <basename> -d <output_dir> -t <threads>"; exit 1; }
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -i) input_file="$2"; shift 2;;
        -base) basename="$2"; shift 2;;
        -t) threads="$2"; shift 2;;
        -d) output_dir="$2"; shift 2;;
        *) echo "Unknown option: $1"; usage;;
    esac
done
[[ -z "$input_file" ]] && { echo "Error: Input file required."; usage; }
[[ ! -f "$input_file" ]] && { echo "Error: $input_file not found."; exit 1; }
threads=${threads:-1}
output_dir=${output_dir:-omark_omamer_output}
[[ -z "$basename" ]] && { basename=$(basename "$input_file"); basename="${basename%.*}"; }
mkdir -p "$output_dir"
echo "Input: $input_file; Base: $basename; Output: $output_dir; Threads: $threads"
awk '/^>/ { seq_id=substr($1,2); split(seq_id,parts,"."); tag=parts[1]; tag_list[tag]=(tag_list[tag]?tag_list[tag]"; "seq_id:seq_id) } END { for(tag in tag_list) print tag_list[tag] }' "$input_file" > "$output_dir/$basename.isoform"
omamer search --db /home/abieskawa/output/LUCA.h5 --query "$input_file" --out "$output_dir/$basename.omamer" -t "$threads"
omark -f "$output_dir/$basename.omamer" -d /home/abieskawa/output/LUCA.h5 -o "$output_dir" -i "$output_dir/$basename.isoform"
