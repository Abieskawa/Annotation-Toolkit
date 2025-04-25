#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import re
import shutil
from Bio.Seq import Seq

def usage():
    print("Usage: {} -g <input_genome> -r <input_gff> -p <file_prefix> -a <annotator> [-v] [--nameprefix <header_prefix>] [--pep_mitos <mitos_fasta>] [--pep_braker <braker_fasta>] [-d <outputdir>]".format(sys.argv[0]))
    sys.exit(1)

def run_command(cmd, shell=False):
    """Run a command and exit if it fails. Capture and display stderr output."""
    print("Running: " + (" ".join(cmd) if not shell else cmd))
    
    try:
        if shell:
            result = subprocess.run(cmd, shell=True, check=True, 
                                   stderr=subprocess.PIPE, text=True)
        else:
            result = subprocess.run(cmd, check=True, 
                                   stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print("Error running command: ", cmd, file=sys.stderr)
        print("Error output:", e.stderr, file=sys.stderr)
        sys.exit(e.returncode)

def extract_bed_from_gff(gff_file, features, bed_files, adjust_start=True, verbose=False):
    """
    Build a hierarchy of gene/transcript/exon/intron from the GFF file and output BED entries
    for the specified features in one pass.
    
    Features can include: "gene", "intron", "transcript".
    
    For transcript-level BED, the transcript boundaries are taken directly from the transcript record,
    without merging exon intervals.
    """
    # Ensure features is a list
    if isinstance(features, str):
        features = [features]
        
    genes = {}         # key: gene id, value: (chrom, start, end, strand, attr)
    transcripts = {}   # key: transcript id, value: dict with keys:
                       #   "transcript": (chrom, start, end, strand, attr),
                       #   "parent_gene": gene id,
                       #   "exons": list of (start, end),
                       #   "introns": list of (start, end)
    warned_genes = set()  # Track genes for which a warning has already been printed

    with open(gff_file) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attr = parts
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            # Parse attributes (assumes key=value; format)
            attributes = {}
            for a in attr.split(";"):
                a = a.strip()
                if not a:
                    continue
                if "=" in a:
                    key, value = a.split("=", 1)
                    attributes[key] = value

            ftype_low = ftype.lower()
            # Gene-level records.
            if ftype_low in ["gene"]:
                gene_id = attributes.get("ID") or attributes.get("gene_id")
                if gene_id:
                    genes[gene_id] = (chrom, start, end, strand, attr)
            # Exon-level records (including CDS, UTR, start_codon, stop_codon).
            elif ftype_low in ["exon", "cds", "utr", "start_codon", "stop_codon"]:
                parent = attributes.get("Parent")
                if parent:
                    if parent in transcripts:
                        transcripts[parent]["exons"].append((start, end))
                    elif parent in genes:
                        if parent not in warned_genes:
                            print(f"WARNING: Exon-level record directly attached to gene {parent} encountered. Ignoring.", file=sys.stderr)
                            warned_genes.add(parent)
                    else:
                        transcripts.setdefault(parent, {"transcript": None, "exons": [], "introns": []})
                        transcripts[parent]["exons"].append((start, end))
            # Intron-level records.
            elif ftype_low == "intron":
                parent = attributes.get("Parent")
                if parent:
                    transcripts.setdefault(parent, {"transcript": None, "exons": [], "introns": []})
                    transcripts[parent]["introns"].append((start, end))
            # Transcript-level records.
            elif ("Parent" in attributes or "gene_id" in attributes) and attributes.get("ID"):
                transcript_id = attributes.get("ID")
                parent_gene = attributes.get("Parent") or attributes.get("gene_id")
                transcripts.setdefault(transcript_id, {"transcript": None, "exons": [], "introns": []})
                transcripts[transcript_id]["transcript"] = (chrom, start, end, strand, attr)
                transcripts[transcript_id]["parent_gene"] = parent_gene

    # Filter transcripts: only keep those with a valid parent gene.
    valid_transcripts = {}
    for tid, t in transcripts.items():
        if "parent_gene" in t:
            if t["parent_gene"] in genes:
                valid_transcripts[tid] = t
            else:
                print(f"WARNING: Transcript {tid} has parent_gene {t['parent_gene']} that is not found in genes. Skipping.", file=sys.stderr)
        else:
            print(f"WARNING: Transcript {tid} does not have a parent_gene attribute. Skipping.", file=sys.stderr)

    if verbose:
        for gene_id, gene_data in genes.items():
            total_exons = 0
            for t in valid_transcripts.values():
                if t.get("parent_gene") == gene_id:
                    total_exons += len(t["exons"])
            if total_exons == 0:
                print(f"DEBUG: Gene {gene_id} has a gene record but no exon-level data.")

    # For each requested feature, write to its corresponding BED file.
    for feat in features:
        feat_lower = feat.lower()
        with open(bed_files[feat_lower], "w") as fout:
            if feat_lower == "transcript":
                for transcript_id, t in valid_transcripts.items():
                    if t["transcript"] is None:
                        continue
                    chrom, tstart, tend, strand, t_attr = t["transcript"]
                    bed_start = tstart - 1 if adjust_start else tstart
                    line_out = f"{chrom}\t{bed_start}\t{tend}\t{transcript_id}\t.\t{strand}"
                    fout.write(line_out + "\n")
            elif feat_lower == "intron":
                for transcript_id, t in valid_transcripts.items():
                    if t.get("introns") and len(t["introns"]) > 0:
                        sorted_introns = sorted(t["introns"], key=lambda x: x[0])
                        for i, (istart, iend) in enumerate(sorted_introns):
                            if t["transcript"]:
                                chrom = t["transcript"][0]
                                strand = t["transcript"][3]
                            else:
                                parent_gene = t.get("parent_gene")
                                chrom, _, _, strand, _ = genes[parent_gene]
                            bed_start = istart - 1 if adjust_start else istart
                            line_out = f"{chrom}\t{bed_start}\t{iend}\t{transcript_id}_intron{i+1}\t.\t{strand}"
                            fout.write(line_out + "\n")
            elif feat_lower == "gene":
                for gene_id, gene_data in genes.items():
                    chrom, gstart, gend, strand, gene_attr = gene_data
                    bed_start = gstart - 1 if adjust_start else gstart
                    line_out = f"{chrom}\t{bed_start}\t{gend}\t{gene_id}\t.\t{strand}"
                    fout.write(line_out + "\n")
    return

def process_pep_file(pep_in, pep_out, header_prefix, pep_type):
    records = []
    with open(pep_in) as fin:
        header = None
        seq_lines = []
        for line in fin:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq_lines)))
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header:
            records.append((header, "".join(seq_lines)))
    
    updated_records = []
    if pep_type == "mitos":
        expected_prefix = f"{header_prefix}_gene_"
        for header, seq in records:
            seq = seq.rstrip("*")
            if header.startswith(">" + expected_prefix):
                updated_records.append((header, seq))
            else:
                # Assume header format: ">MT; 487-1453; +; nad1"
                parts = header[1:].split(";")
                gene = parts[-1].strip().lower()
                new_header = f">{expected_prefix}{gene}.p"
                updated_records.append((new_header, seq))
        print("Processed mitos peptide input.")
    elif pep_type == "braker":
        expected_prefix = f"{header_prefix}_"
        for header, seq in records:
            seq = seq.rstrip("*")
        
            # Get the ID part without the ">"
            original_id = header[1:]
        
            # Replace .t with .p instead of appending .p
            if ".t" in original_id:
                modified_id = original_id.replace(".t", ".p")
            else:
                # If no .t pattern found, just append .p
                modified_id = original_id + ".p"
            
            # Check if prefix needs to be added
            if not modified_id.startswith(expected_prefix):
                new_header = f">{expected_prefix}{modified_id}"
            else:
                new_header = f">{modified_id}"
            
            updated_records.append((new_header, seq))
        print("Processed braker peptide input.")

    else:
        print("Unknown peptide type specified.", file=sys.stderr)
        sys.exit(1)
    with open(pep_out, "w") as fout:
        for rec in updated_records:
            fout.write(f"{rec[0]}\n{rec[1]}\n")
    print(f"Peptide FASTA processed and written to: {pep_out}")

def translate_geseq_cds(cds_fa: str, gene_fa: str, pep_out: str, table):
    """
    Translate CDS FASTA produced by gffread for GeSeq annotations.
    Uses gene names from gene FASTA file and adds .p to create protein IDs.
    Only sequences whose core ID matches an entry in gene_mapping are processed.
    Reports and trims sequences with partial codons.
    """
    import re
    from Bio.Seq import Seq
    import sys

    # 1) extract all gene names from the gene FASTA file
    gene_names = []
    with open(gene_fa) as f:
        for line in f:
            if line.startswith('>'):
                gene_names.append(line.strip()[1:])  # drop '>'

    # 2) build a mapping from core gene name to full gene name
    #    e.g. core "ND4" → full "M_tai_gene-blatx_ND4_1"
    gene_mapping = {}
    for gene in gene_names:
        m = re.search(r'gene-blatx_([^_]+)', gene)
        if m:
            gene_mapping[m.group(1).upper()] = gene

    # 3) translate only those CDS records that map
    with open(cds_fa) as fh, open(pep_out, "w") as out:
        hdr, seq = None, []
        for ln in fh:
            if ln.startswith(">"):
                # process the previous record
                if hdr is not None:
                    # core ID before the first dot, e.g. "M_tai_g9894"
                    core_id = hdr.split('.', 1)[0]
                    matched = next(
                        (full for core, full in gene_mapping.items()
                         if core in core_id.upper()),
                        None
                    )
                    if matched:
                        full_seq = "".join(seq)
                        rem = len(full_seq) % 3
                        if rem != 0:
                            print(f"WARNING: Partial codon detected in sequence {hdr}", file=sys.stderr)
                            print(f" Sequence length: {len(full_seq)} (not divisible by 3)", file=sys.stderr)
                            print(f" Trimming {rem} nucleotide(s) from the end", file=sys.stderr)
                            full_seq = full_seq[:-rem]

                        prot = str(Seq(full_seq).translate(table=table,
                                                           to_stop=True,
                                                           cds=False))
                        out.write(f">{matched}.p\n{prot}\n")
                    # unmatched records are skipped

                # start a new record
                hdr = ln[1:].strip()
                seq = []
            else:
                seq.append(ln.strip())

        # process the last record
        if hdr is not None:
            core_id = hdr.split('.', 1)[0]
            matched = next(
                (full for core, full in gene_mapping.items()
                 if core in core_id.upper()),
                None
            )
            if matched:
                full_seq = "".join(seq)
                rem = len(full_seq) % 3
                if rem != 0:
                    print(f"WARNING: Partial codon detected in sequence {hdr}", file=sys.stderr)
                    print(f" Sequence length: {len(full_seq)} (not divisible by 3)", file=sys.stderr)
                    print(f" Trimming {rem} nucleotide(s) from the end", file=sys.stderr)
                    full_seq = full_seq[:-rem]

                prot = str(Seq(full_seq).translate(table=table,
                                                   to_stop=True,
                                                   cds=False))
                out.write(f">{matched}.p\n{prot}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline using a GFF input (converted to GTF) to generate BED/FASTA outputs at gene, intron, and transcript levels"
    )
    parser.add_argument("-g", "--genome", required=True, help="Input genome file")
    parser.add_argument("-r", "--gff", required=True, help="Input GFF file")
    parser.add_argument("-p", "--prefix", required=True, help="File prefix for output files (e.g. T_jap)")
    parser.add_argument("-a", "--annotator", required=True, help="Annotator name (e.g. braker3)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose debug output")
    parser.add_argument("--nameprefix", help="Header prefix for FASTA entries. Defaults to file prefix if not provided.", default=None)
    parser.add_argument("--pep_mitos", help="Optional mitos peptide FASTA file provided by the user", default=None)
    parser.add_argument("--pep_braker", help="Optional braker peptide FASTA file provided by the user", default=None)
    parser.add_argument("--geseq", action="store_true",
                help="Translate GeSeq‑style CDS records in the generated CDS FASTA")
    parser.add_argument("--table", type=int, required="--geseq" in sys.argv,
                help="Genetic‑code table ID/name for Bio.Seq.translate; required with --geseq")
    parser.add_argument("-d", "--outputdir", help="Specify the output directory name (default: MOLAS_input)", default="MOLAS_input")
    args = parser.parse_args()

    input_genome = args.genome
    input_gff = args.gff
    file_prefix = args.prefix
    annotator = args.annotator
    verbose = args.verbose
    header_prefix = args.nameprefix if args.nameprefix else file_prefix

    # Define output file names.
    output_intron_lowercase_genome = f"{file_prefix}_genome_intron-lower.fa"
    intron_bed = f"{file_prefix}_{annotator}_intron.bed"
    gene_bed = f"{file_prefix}_{annotator}_gene.bed"
    trans_bed = f"{file_prefix}_{annotator}_trans.bed"
    gene_fa = f"{file_prefix}_{annotator}_gene.fa"
    trans_fa = f"{file_prefix}_{annotator}_trans.fa"
    cds_fa_temp = f"{file_prefix}_{annotator}_CDS_70.fa"
    pep_fa_final = f"{file_prefix}_{annotator}_pep.fa"

    # -----------------------------------------------------
    # Generate a GTF file from the input GFF using gffread.
    gtf_generated = "merged_fix_tran.gtf"
    cmd_gffread_gtf = ["gffread", "-T", input_gff, "-o", gtf_generated]
    run_command(cmd_gffread_gtf)
    print(f"GTF file generated: {gtf_generated}")
    # We'll use the generated GTF for any downstream processing if needed.
    # For BED extractions, we work with the original GFF.
    # -----------------------------------------------------

    # Generate BED files for intron, gene, and transcript in one pass.
    bed_files = {
        "intron": intron_bed,
        "gene": gene_bed,
        "transcript": trans_bed
    }
    features = ["intron", "gene", "transcript"]
    extract_bed_from_gff(input_gff, features, bed_files, adjust_start=True, verbose=verbose)
    print(f"Intron BED file written: {intron_bed}")
    print(f"Gene BED file written: {gene_bed}")
    print(f"Transcript BED file written: {trans_bed}")

    # Use bedtools maskfasta to soft-mask the genome using the intron BED file.
    cmd_mask = [
        "bedtools", "maskfasta", "-soft",
        "-fi", input_genome,
        "-bed", intron_bed,
        "-fo", output_intron_lowercase_genome
    ]
    run_command(cmd_mask)
    print(f"Soft-masked genome written: {output_intron_lowercase_genome}")
    # -----------------------------------------------------

    # Use bedtools getfasta to extract gene and transcript sequences.
    cmd_getfasta_gene = [
        "bedtools", "getfasta", "-fi", output_intron_lowercase_genome,
        "-bed", gene_bed, "-nameOnly", "-s", "-fo", gene_fa
    ]
    run_command(cmd_getfasta_gene)
    print(f"Gene FASTA written: {gene_fa}")
    cmd_getfasta_trans = [
        "bedtools", "getfasta", "-fi", output_intron_lowercase_genome,
        "-bed", trans_bed, "-nameOnly", "-s", "-fo", trans_fa
    ]
    run_command(cmd_getfasta_trans)
    print(f"Transcript FASTA written: {trans_fa}")

    # Remove any "(+)" and "(-)" from the FASTA headers.
    for fasta_file in [gene_fa, trans_fa]:
        with open(fasta_file, "r") as f:
            content = f.read()
        content = re.sub(r'\(\+\)', '', content)
        content = re.sub(r'\(-\)', '', content)
        with open(fasta_file, "w") as f:
            f.write(content)
        print(f"Cleaned {fasta_file} of stranded annotations.")

    # Run gffread (using the generated GTF) to extract CDS FASTA sequences only.
    cmd_gffread = [
        "gffread", "-g", output_intron_lowercase_genome,
        "-x", cds_fa_temp,
        gtf_generated
    ]
    run_command(cmd_gffread)
    print(f"gffread produced {cds_fa_temp}")

    # Process CDS FASTA file (simulate the awk/sort/awk pipeline).
    cds_fa = cds_fa_temp.replace("_70", "")
    records = []
    with open(cds_fa_temp, "r") as fin:
        header = None
        seq_lines = []
        for line in fin:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq_lines)))
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header:
            records.append((header, "".join(seq_lines)))
    records.sort(key=lambda x: x[0])
    with open(cds_fa, "w") as fout:
        for rec in records:
            fout.write(f"{rec[0]}\n{rec[1]}\n")
    print(f"Processed CDS FASTA file written: {cds_fa}")

    # ----- Peptide FASTA Processing -----
    processed_pep_files = []
    if args.pep_mitos:
        out_pep_mitos = f"{file_prefix}_{annotator}_pep_mitos.fa"
        process_pep_file(args.pep_mitos, out_pep_mitos, header_prefix, "mitos")
        processed_pep_files.append(out_pep_mitos)
    if args.geseq:
       if Seq is None:
          sys.exit("ERROR: Biopython is not installed but --geseq was requested")
       geseq_pep = f"{file_prefix}_{annotator}_pep_geseq.fa"
       translate_geseq_cds(cds_fa, gene_fa, geseq_pep, args.table)
       processed_pep_files.append(geseq_pep)

    if args.pep_braker:
        out_pep_braker = f"{file_prefix}_{annotator}_pep_braker.fa"
        process_pep_file(args.pep_braker, out_pep_braker, header_prefix, "braker")
        processed_pep_files.append(out_pep_braker)
    if processed_pep_files:
        final_pep = pep_fa_final  # final peptide FASTA file name
        with open(final_pep, "w") as fout:
            for f in processed_pep_files:
                with open(f) as fin:
                    fout.write(fin.read())
        print(f"Final peptide FASTA written (concatenated): {final_pep}")
    else:
        print("No peptide FASTA input provided; skipping peptide processing.")

    # Clean up temporary BED files and the temporary GTF file.
    for f in [intron_bed, gene_bed, trans_bed, cds_fa_temp]:
        if os.path.exists(f):
            os.remove(f)
    # Clean up processed peptide files individually.
    for f in processed_pep_files:
        if os.path.exists(f):
            os.remove(f)
    print("Temporary BED files removed.")
    if os.path.exists(gtf_generated):
        os.remove(gtf_generated)
        print(f"Removed temporary GTF file: {gtf_generated}")

    # Move output files into user-specified output directory.
    output_dir = f"./{args.outputdir}/FASTA/"
    os.makedirs(output_dir, exist_ok=True)
    files_to_move = [
        output_intron_lowercase_genome,
        output_intron_lowercase_genome + ".fai",
        gene_fa, trans_fa, cds_fa, pep_fa_final]
    for f in files_to_move:
        if os.path.exists(f):
            shutil.move(f, os.path.join(output_dir, os.path.basename(f)))
    print(f"Moved output files to {output_dir}")

if __name__ == "__main__":
    main()
