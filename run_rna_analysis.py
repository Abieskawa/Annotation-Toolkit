#!/usr/bin/env python3
"""
STAR ➜ StringTie versatile RNA‑seq pipeline

Outputs (you may request any combination):
  --transcript-fpkm    transcript_expression.tab  (t_data.ctab → sort 6, cut 6,12)
  --transcript-counts  transcript_counts.tab      (t_data.ctab → sort 1, cut 1,2)
  --gene-fpkm          gene_abundance_estimation.tab (abundance.tab → sort 1, cut 1,8)
  --gene-tpm           gene_tpm.tab                  (abundance.tab → sort 1, cut 1,9)
  --gene-counts        gene_counts.tab               (STAR ReadsPerGene.out.tab → sort 1, cut 1,2)
Each mode is written to its own subdirectory under --out-dir.
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


def run(cmd, *, use_shell=False, **kw):
    """If use_shell=True run through /bin/bash so pipes/redirects work."""
    if use_shell:
        print('+', cmd)
        subprocess.run(cmd, shell=True, check=True,
                       executable='/bin/bash', **kw)
    else:
        print('+', ' '.join(cmd))
        subprocess.run(cmd, check=True, **kw)

def bg_path(outd, sample, fname):
    """
    Find Ballgown / abundance files regardless of layout.
    1) <outd>/<sample>/<fname>             (# safest: per‑sample dir)
    2) <outd>/<sample>.ballgown/<fname>    (# legacy StringTie layout)
    3) <outd>/<fname>                      (# single‑sample run)
    """
    for c in (
        os.path.join(outd, sample, fname),
        os.path.join(outd, f"{sample}.ballgown", fname),
        os.path.join(outd, fname),
    ):
        if os.path.exists(c):
            return c
    raise FileNotFoundError(f"{fname} not found for {sample}")

def run_version(cmd):
    print('#', ' '.join(cmd))
    subprocess.run(cmd, check=False)

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description='STAR → StringTie multi‑mode pipeline')
    ap.add_argument('--reads-dir', required=True)
    ap.add_argument('--genome-fasta', required=True)
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--threads', type=int, default=4)

    ap.add_argument('--transcript-fpkm',   action='store_true')
    ap.add_argument('--transcript-counts', action='store_true')
    ap.add_argument('--gene-fpkm',         action='store_true')
    ap.add_argument('--gene-tpm',          action='store_true')
    ap.add_argument('--gene-counts',       action='store_true')
    ap.add_argument('--merge-expression',  action='store_true')
    args = ap.parse_args()

    if not any([args.transcript_fpkm, args.transcript_counts,
                args.gene_fpkm, args.gene_tpm, args.gene_counts]):
        sys.exit('Error: select at least one output mode')

    star_index = os.path.join(args.out_dir, 'STAR_index')
    align_dir  = os.path.join(args.out_dir, 'alignments'); os.makedirs(align_dir, exist_ok=True)

    mode_dirs = {}
    if args.transcript_fpkm:   mode_dirs['transcript_fpkm']   = os.path.join(args.out_dir,'transcript_fpkm')
    if args.transcript_counts: mode_dirs['transcript_counts'] = os.path.join(args.out_dir,'transcript_counts')
    if args.gene_fpkm:         mode_dirs['gene_fpkm']         = os.path.join(args.out_dir,'gene_fpkm')
    if args.gene_tpm:          mode_dirs['gene_tpm']          = os.path.join(args.out_dir,'gene_tpm')
    if args.gene_counts:       mode_dirs['gene_counts']       = os.path.join(args.out_dir,'gene_counts')
    for d in mode_dirs.values(): os.makedirs(d, exist_ok=True)

    ann = args.annotation; is_gff = ann.lower().endswith(('.gff','.gff3'))

    for tool in [['STAR','--version'],['stringtie','--version'],['sort','--version']]:
        run_version(tool)
    print()

    idx_file = os.path.join(star_index, 'SAindex')
    print(f"looking for STAR index at {idx_file!r}: exists? {os.path.exists(idx_file)}", flush=True)
    if not os.path.exists(idx_file):
        print("STAR index not found, building STAR index...")
        os.makedirs(star_index, exist_ok=True)
        idx = ['STAR','--runMode','genomeGenerate','--genomeFastaFiles',args.genome_fasta,
               '--genomeSAindexNbases',str(compute_genomeSAindex_nbases(args.genome_fasta)),
               '--runThreadN',str(args.threads),'--genomeDir',star_index,'--sjdbGTFfile',ann]
        if is_gff:
            idx += ['--sjdbGTFfeatureExon','exon','--sjdbGTFtagExonParentTranscript','Parent','--sjdbGTFtagExonParentGene','Parent']
        run(idx)
    else:
        print("STAR index exists, skipping build", flush=True)

    # locate reads
    r1_files = sorted(glob.glob(os.path.join(args.reads_dir,'*_1.cleaned.fastq*')))
    samples=[]
    for r1 in r1_files:
        r2 = r1.replace('_1.cleaned.fastq.gz','_2.cleaned.fastq.gz') if r1.endswith('.gz') else r1.replace('_1.cleaned.fastq','_2.cleaned.fastq')
        if not os.path.exists(r2):
            print('skip',r1,'no R2',file=sys.stderr); continue
        sample = os.path.basename(r1).split('_1.cleaned.fastq')[0].rstrip('.gz'); samples.append(sample)

        # STAR mapping prefix and BAM path
        pref = os.path.join(align_dir, sample + '.')
        bam  = pref + 'Aligned.sortedByCoord.out.bam'
        if not os.path.exists(bam):
            star_cmd = [
                'STAR', '--runMode', 'alignReads', '--twopassMode', 'Basic',
                '--runThreadN', str(args.threads),
                '--genomeDir', star_idx,
                '--readFilesIn', r1, r2,
                '--outSAMtype', 'BAM', 'SortedByCoordinate',
                '--quantMode', 'TranscriptomeSAM', 'GeneCounts',
                '--outFileNamePrefix', pref,
                '--outSAMstrandField', 'intronMotif'
            ]
            if r1.endswith('.gz'):
                star_cmd += ['--readFilesCommand', 'zcat']
            run(star_cmd)
        else:
            print(f"Skipping STAR for sample {sample}, BAM exists: {bam}")

        # ─── new: filter for properly paired reads ───
        bam_flag2 = pref + 'flag2.bam'
        samtools_cmd = [
            'samtools', 'view', '-b', '-f', '2',
            '-@', str(args.threads),
            bam,
            '-o', bam_flag2
        ]
        run(samtools_cmd)
        # now use the filtered BAM for downstream steps
        bam = bam_flag2

        rpg  = pref+'ReadsPerGene.out.tab'

        # per‑mode work
        for mode,outd in mode_dirs.items():
            raw_dir   = os.path.join(stringtie_dir, mode, sample)
            final_dir = os.path.join(final_root, sample)           # final tables remain under '<mode>/<sample>'
            os.makedirs(raw_dir, exist_ok=True)
            os.makedirs(final_dir, exist_ok=True)
            # StringTie abundance modes
            # run StringTie (raw outputs)
            if mode != 'gene_counts':
                abundance = os.path.join(raw_dir, 'abundance.tab')
                st_cmd = [
                    'stringtie', bam,
                    '-G', args.annotation,
                    '-o', os.path.join(raw_dir, 'stringtie.gtf'),
                    '-p', str(args.threads),
                    '-e',
                    '-b', raw_dir,
                    '-f', '0.15','-m','200','-a','10','-j','1','-c','2','-g','50','-M','0.95',
                    '-A', abundance
                ]
                run(st_cmd)

            # Determine source file/columns
            # extract tables into final_dir
            if mode == 'transcript_fpkm':
                src = os.path.join(raw_dir, sample, 't_data.ctab')    
                run(f"sort -k6,6 {src} | cut -f6,12 > {os.path.join(final_dir, sample + '.transcript_fpkm.tab')}" , 
                use_shell=True)
            elif mode == 'transcript_counts':
                src = os.path.join(raw_dir, sample, 't_data.ctab')   
                run(f"sort -k1,1 {src} | cut -f1,2 > {os.path.join(final_dir, sample + '.transcript_counts.tab')}" , 
                use_shell=True)
            elif mode == 'gene_fpkm':
                src = os.path.join(raw_dir, 'abundance.tab')
                run(f"sort -k1,1 {src} | cut -f1,8 > {os.path.join(final_dir, sample + '.gene_fpkm.tab')}" , 
                use_shell=True)
            elif mode == 'gene_tpm':
                src = os.path.join(raw_dir, 'abundance.tab')
                run(f"sort -k1,1 {src} | cut -f1,9 > {os.path.join(final_dir, sample + '.gene_tpm.tab')}" , 
                use_shell=True)
            else:  # gene_counts
                rpg = pref + 'ReadsPerGene.out.tab'
                run(f"sort -k1,1 {rpg} | cut -f1,2 > {os.path.join(final_dir, sample + '.gene_counts.tab')}" , 
                use_shell=True)

    # merge for MOLAS
    if args.merge_expression:
        table = os.path.join(args.out_dir,'expression_table.tab')
        with open(table,'w') as out: out.write('ID\t'+'\t'.join(samples)+'\n')
        collect={}
        for mode,d in mode_dirs.items():
            for s in samples:
                fn=os.path.join(d,f'{s}.{mode}.tab')
                if not os.path.exists(fn): continue
                with open(fn) as fh:
                    next(fh)
                    for ln in fh:
                        gid,val=ln.strip().split('\t',1)
                        collect.setdefault(gid,{})[s]=val
        with open(table,'a') as out:
            for gid in sorted(collect):
                out.write(gid+'\t'+'\t'.join(collect[gid].get(s,'.') for s in samples)+'\n')

if __name__=='__main__':
    main()
