#!/usr/bin/env python3
"""
STAR ➔ StringTie versatile RNA seq pipeline

Outputs (you may request any combination):
  --transcript-fpkm    per-sample transcript FPKM values
  --gene-fpkm          per-sample gene FPKM values
  --gene-tpm           per-sample gene TPM values
  --transcript-counts  transcript count matrix via prepDE.py per sample
  --gene-counts        gene count matrix via prepDE.py per sample
  --map-only           perform mapping and filtering only; skip StringTie and prepDE
  --merge-expression   merge per-mode tables into a single expression table
Each output mode writes into its own directory under --out-dir; per-sample count matrices are placed at --out-dir.
"""
import argparse, subprocess, glob, os, sys, math

# ---------- helpers ----------

def compute_genomeSAindex_nbases(fasta):
    total = 0
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith('>'): continue
            total += len(ln.strip())
    val = int(math.log2(total) / 2 - 1)
    return val if val < 14 else 14


def run(cmd, *, cwd=None, use_shell=False, **kw):
    if use_shell:
        print(cmd, flush=True)
        subprocess.run(cmd, shell=True, check=True, cwd=cwd,
                       executable='/bin/bash', **kw)
    else:
        print(' '.join(cmd), flush=True)
        subprocess.run(cmd, check=True, cwd=cwd, **kw)


def run_version(cmd):
    print('#', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=False)

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description='STAR → StringTie multi-mode pipeline')
    ap.add_argument('--reads-dir', required=True)
    ap.add_argument('--genome-fasta', required=True)
    ap.add_argument('--annotation', required=True,
                    help='GTF annotation file (GFF input is deprecated)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--map-only', action='store_true',
                    help='Only perform mapping and filtering; skip StringTie and prepDE')
    ap.add_argument('--transcript-fpkm',   action='store_true')
    ap.add_argument('--gene-fpkm',         action='store_true')
    ap.add_argument('--gene-tpm',          action='store_true')
    ap.add_argument('--transcript-counts', action='store_true')
    ap.add_argument('--gene-counts',       action='store_true')
    ap.add_argument('--merge-expression',  action='store_true')
    args = ap.parse_args()

    print('Command invoked:', ' '.join(sys.argv), flush=True)

    if args.annotation.lower().endswith(('.gff', '.gff3')):
        sys.exit('Error: GFF input is not suggested to run STAR; please provide a GTF (*.gtf) file.')

    if not args.map_only and not any([args.transcript_fpkm, args.gene_fpkm,
                                      args.gene_tpm, args.transcript_counts,
                                      args.gene_counts]):
        sys.exit('Error: select at least one output mode or use --map-only')

    star_index    = os.path.join(args.out_dir, 'STAR_index')
    align_dir     = os.path.join(args.out_dir, 'alignments'); os.makedirs(align_dir, exist_ok=True)
    stringtie_dir = os.path.join(args.out_dir, 'stringtie'); os.makedirs(stringtie_dir, exist_ok=True)

    mode_dirs = {}
    if args.transcript_fpkm:   mode_dirs['transcript_fpkm']   = os.path.join(args.out_dir, 'transcript_fpkm')
    if args.transcript_counts: mode_dirs['transcript_counts'] = os.path.join(args.out_dir, 'transcript_counts')
    if args.gene_fpkm:         mode_dirs['gene_fpkm']         = os.path.join(args.out_dir, 'gene_fpkm')
    if args.gene_tpm:          mode_dirs['gene_tpm']          = os.path.join(args.out_dir, 'gene_tpm')
    if args.gene_counts:       mode_dirs['gene_counts']       = os.path.join(args.out_dir, 'gene_counts')
    for d in mode_dirs.values(): os.makedirs(d, exist_ok=True)

    for tool in [['STAR','--version'], ['stringtie','--version'], ['samtools','--version'], ['sort','--version']]:
        run_version(tool)
    print(flush=True)

    idx_file = os.path.join(star_index, 'SAindex')
    if os.path.exists(idx_file):
        print(f"STAR index exists.", flush=True)
    if not os.path.exists(idx_file):
        cmd = [
            'STAR','--runMode','genomeGenerate',
            '--genomeFastaFiles', args.genome_fasta,
            '--genomeSAindexNbases', str(compute_genomeSAindex_nbases(args.genome_fasta)),
            '--runThreadN', str(args.threads),
            '--genomeDir', star_index,
            '--sjdbGTFfile', args.annotation
        ]
        print('Building STAR index with command:', ' '.join(cmd), flush=True)
        os.makedirs(star_index, exist_ok=True)
        run(cmd)
    else:
        print("STAR index exists, skipping build", flush=True)

    # process each sample (support both .fastq* and .fq* extensions)
    fastq1 = glob.glob(os.path.join(args.reads_dir, '*_1.cleaned.fastq*'))
    fq1    = glob.glob(os.path.join(args.reads_dir, '*_1.cleaned.fq*'))
    r1_files = sorted(set(fastq1 + fq1))
    samples = []
    for r1 in r1_files:
        r2 = r1.replace('_1.cleaned', '_2.cleaned')
        paired = os.path.exists(r2)
        if not paired:
            print(f"Single-end sample for {sample}, no R2 found: {r2}", flush=True)
        sample = os.path.basename(r1).split('_1.cleaned')[0]
        samples.append(sample)

        pref = os.path.join(align_dir, sample + '.')
        bam  = pref + 'Aligned.sortedByCoord.out.bam'
        if not os.path.exists(bam):
            cmd = [
                'STAR','--runMode','alignReads','--twopassMode','Basic',
                '--runThreadN', str(args.threads),
                '--genomeDir', star_index,
                '--outSAMtype','BAM','SortedByCoordinate',
                '--quantMode','TranscriptomeSAM','GeneCounts',
                '--outFileNamePrefix', pref,
                '--outSAMstrandField','intronMotif'
            ]
            if paired:
                cmd += ['--readFilesIn', r1, r2]
            else:
                cmd += ['--readFilesIn', r1]

            if r1.endswith('.gz'):
                cmd += ['--readFilesCommand','zcat']
            print(f"Running STAR alignment for {sample} with command:", ' '.join(cmd), flush=True)
            run(cmd)
        else:
            print(f"Skipping STAR for {sample}, BAM exists", flush=True)

        bam_flag2 = pref + 'flag2.bam'
        samtools_cmd = ['samtools','view','-b','-f','2','-@', str(args.threads), bam, '-o', bam_flag2]
        print('Filtering BAM with command:', ' '.join(samtools_cmd), flush=True)
        run(samtools_cmd)
        
        if args.map_only:
            continue

        # run StringTie
        raw_dir = os.path.join(stringtie_dir, sample); os.makedirs(raw_dir, exist_ok=True)
        cmd = [
            'stringtie', bam_flag2, '-G', args.annotation,
            '-o', os.path.join(raw_dir,'stringtie.gtf'),
            '-p', str(args.threads), '-e', '-l', sample,
            '-b', raw_dir, '-f','0.15','-m','200','-a','10',
            '-j','1','-c','2','-g','50','-M','0.95',
            '-A', os.path.join(raw_dir,'abundance.tab')
        ]
        run(cmd)

        if args.transcript_fpkm:
            src = os.path.join(raw_dir,'t_data.ctab')
            out = os.path.join(mode_dirs['transcript_fpkm'], sample + '.transcript_fpkm.tab')
            run(f"sort -k6,6 {src} | cut -f6,12 > {out}", use_shell=True)
        if args.gene_fpkm:
            src = os.path.join(raw_dir,'abundance.tab')
            out = os.path.join(mode_dirs['gene_fpkm'], sample + '.gene_fpkm.tab')
            run(f"sort -k1,1 {src} | cut -f1,8 > {out}", use_shell=True)
        if args.gene_tpm:
            src = os.path.join(raw_dir,'abundance.tab')
            out = os.path.join(mode_dirs['gene_tpm'], sample + '.gene_tpm.tab')
            run(f"sort -k1,1 {src} | cut -f1,9 > {out}", use_shell=True)

    # generate count matrices
    if args.transcript_counts or args.gene_counts:
        mapping_file = os.path.join(args.out_dir, 'prepDE_input.txt')
        with open(mapping_file, 'w') as fh:
            for sample in samples:
                gtf_file = os.path.abspath(os.path.join(stringtie_dir, sample, 'stringtie.gtf'))
                fh.write(f"{sample}\t{gtf_file}\n")
                print('Mapping entry:', sample, '->', gtf_file, flush=True)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        prepde_script = os.path.join(script_dir, 'prepDE.py')
        prepde_cmd = ['python3', '-W', 'ignore::SyntaxWarning', prepde_script, '-i', mapping_file]
        if args.gene_counts:
            gene_matrix = os.path.join(mode_dirs['gene_counts'], 'gene_count_matrix.csv')
            prepde_cmd += ['-g', gene_matrix]
        if args.transcript_counts:
            trans_matrix = os.path.join(mode_dirs['transcript_counts'], 'transcript_count_matrix.csv')
            prepde_cmd += ['-t', trans_matrix]
        print('Running prepDE command:', ' '.join(prepde_cmd), flush=True)
        try:
            run(prepde_cmd)
            print("Successfully generated count matrices", flush=True)
            # create sorted copies
            if args.gene_counts:
                sorted_gene_counts = os.path.join(mode_dirs['gene_counts'], 'gene_count_matrix.sorted.csv')
                run(f"head -n1 {gene_matrix} > {sorted_gene_counts}; tail -n+2 {gene_matrix} | sort -t, -k1,1 >> {sorted_gene_counts}", use_shell=True)
            if args.transcript_counts:
                sorted_trans_counts = os.path.join(mode_dirs['transcript_counts'], 'transcript_count_matrix.sorted.csv')
                run(f"head -n1 {trans_matrix} > {sorted_trans_counts}; tail -n+2 {trans_matrix} | sort -t, -k1,1 >> {sorted_trans_counts}", use_shell=True)
        except subprocess.CalledProcessError as e:
            sys.exit(e.returncode)

    # merge per-mode tables if requested
    if args.merge_expression:
        table = os.path.join(args.out_dir, 'expression.table')
        data = {}
        cols = []
        # per-sample FPKM/TPM
        for mode, suffix in [('gene_fpkm','.g.fpkm'), ('gene_tpm','.g.tpm'), ('transcript_fpkm','.t.fpkm')]:
            if getattr(args, mode):
                mode_dir = os.path.join(args.out_dir, mode)
                for sample in samples:
                    fname = os.path.join(mode_dir, f"{sample}.{mode}.tab")
                    col = f"{sample}{suffix}"
                    cols.append(col)
                    if not os.path.exists(fname):
                        continue
                    with open(fname) as fh:
                        first = fh.readline()
                        fields = first.strip().split('\t')
                        is_data = False
                        for v in fields[1:]:
                            try:
                                float(v)
                                is_data = True
                                break
                            except ValueError:
                                continue
                        if is_data:
                            lines = [first] + fh.readlines()
                        else:
                            lines = fh
                        for ln in lines:
                            id_, val = ln.strip().split('\t',1)
                            data.setdefault(id_, {})[col] = val
        # gene counts matrix
        if args.gene_counts:
            if os.path.exists(sorted_gene_counts):
                with open(sorted_gene_counts) as fh:
                    header_line = next(fh).strip()
                    hdr = header_line.split(',')[1:]
                    for sample in hdr:
                        col = f"{sample}.g.count"
                        cols.append(col)
                    for ln in fh:
                        parts = ln.strip().split(',')
                        id_ = parts[0]
                        if id_.startswith("<class"):
                            print(f"DEBUG: merging phantom ID '{id_}' from file {sorted_gene_counts}", file=sys.stderr)
                            continue
                        for sample, val in zip(hdr, parts[1:]):
                            col = f"{sample}.g.count"
                            data.setdefault(id_, {})[col] = val
        # transcript counts matrix
        if args.transcript_counts:
            if os.path.exists(sorted_trans_counts):
                with open(sorted_trans_counts) as fh:
                    hdr = next(fh).strip().split(',')[1:]
                    for sample in hdr:
                        col = f"{sample}.t.count"
                        cols.append(col)
                    for ln in fh:
                        parts = ln.strip().split(',')
                        id_ = parts[0]
                        for sample, val in zip(hdr, parts[1:]):
                            col = f"{sample}.t.count"
                            data.setdefault(id_, {})[col] = val
        # write merged table
        with open(table, 'w') as fh:
            fh.write('ID	' + '	'.join(cols) + '\n')
            for id_ in sorted(data):
                row = [data[id_].get(c, '.') for c in cols]
                fh.write(id_ + '\t' + '\t'.join(row) + '\n')
    print('Everything Finished.')

if __name__=='__main__':
    main()
