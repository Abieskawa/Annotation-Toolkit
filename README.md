# Annotation-Toolkit
This repository is the home place and notes for any tool applied to annotated the genome, including the command line, evaluation and download tool. While NCBI offers annotation services, users may experience delays. Also, their annotation tools are not open source until now, so a toolkit for annotation can allow the users to adjust according to their requirement. This toolkit includes structural annotation for protein-coding, ncRNA, mitochondria genes. Codes for functional annotation for protein-coding gene are also recorded here. We recommend read and run everything with your genome first to know every detail and potential issues.

This script has been tested on beltfish, tialpia 2 strains, asian hard clam and japanese eel.

## File/Directory structure
It is not necessary to follow this naming rule, but using this naming logic can make things tidier and easier to follow.
```
file basename {ex. asian_hard_clam}
├──00_raw_genome
├──00_check_mitochondria
├──00_check_mitochondria_blastn
├──00_fcs_adaptor_mito_nuclear
├──00_fcs_adaptor_double_check_mito_nuclear
├──00_fcs_gx_mito_nuclear
├──00_fcs_gx_double_check_mito_nuclear
├──00_start_annotation
├──01_RepeatMasker
├──01_RepeatModeler
├──02_RepeatMasker
├──03_SoftMask
├──04_TRF_mask
├──05_raw_RNA
├──05_cleaned_RNA
├──05_RNA_mapping_star_2pass
├──star_index
├──06_braker3
├──06_mitoz/06_mitohifi/06_mitos
├──06_ncRNA
├──07_merge_gff
├──08_FASTA
├──09_deeptmhmm_mito_nuclear
├──10_kaas_mito_nuclear
├──11_diamond_blastp_mito_nulcear
├──12_interproscan_mito_nulcear
└──13_RNA_ana
``` 
Database file structure
```
output
├──file basename {Working dir}
├──diamond_nr
├──fcs-db
├──LUCA.h5
├──Lib_fish {Repeat DB, constructed with instruction in [tetools](https://github.com/Dfam-consortium/TETools)}
├──MITOS2-refdata
└──Rfam_db
```
# Environment


# Check the quality of received genome 
Check the assembly FCS-gx/adapter, mitochondrial contamination, and information, BUSCO, (and clean if it is necessary)

## Check if mitochondria genome exists
[*mitoZ*](https://github.com/linzhi2013/MitoZ) v3.6
```
mitoz findmitoscaf --fastafile {genome, gz supported here} --outprefix {prefix} \
                   --workdir {dir} --thread_number {threads} \
                   --requiring_taxa {Chordata,Arthropoda...} \
                   --clade {Chordata,Arthropoda...} --min_abundance 0
```
## Extract mitochondria and nucleus chromosomes (Longest N)
```
python {the path of this directory}/extract_nucleus_query_genome.py \
        -i {the genome after cleaning with fcs} -on {nuclear fasta name} \
        -os {optional, fasta name of mitochondria or chloreplast} \
        -n {the longest n seqs} -s {mitochondria,chloroplast,or the name of the target genome}
``` 

## Assmeble and annotate mitochondrial genome
NCBI genetic code for mitochondria (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) 
Ex. fish genetic code: 2, asian hard clam genetic code: 5
Other tools, such as GeSeq, mitoZ, MFannot, DeGeCI can also applied here. In my point of view, MitoZ, mitos is more recommended here, expecially for those non-model species, such as asian hard clam.\
[*mitos*](https://gitlab.com/Bernt/MITOS): use relative obsolete reference file, but it output gff and provide other complete useful data, ex. protein sequence. The script for self-configured reference data is removed by the author, so user cannot update the reference by themselves.\
[*mitoZ*](https://github.com/linzhi2013/MitoZ): It can manually add the reference sequence, but output only genebank file, nor gff or gtf. User can use all subcommend to start from Illumina WGS data to mitochondria genome and .gbk annotation file.\
[*MitoHiFi*](https://github.com/marcelauliano/MitoHiFi?tab=readme-ov-file): It can start from PacBio HiFi reads to mito genome and .gbk annotation file. Note that MitoHiFi 3.0.0_cv1 in [DockerHub](https://hub.docker.com/r/biocontainers/mitohifi/tags) is the latest container executable. However, I found that the gene prediction result of MitoHiFi may only have on the same strand, which is weird. As the result, i would recommend to use MitoZ.

Program: mitoZ v3.6\
bio gff cannot accomodate complement(strand:-) info, so we will deal with it in the stage of fix and merge gffs later.
```
#Assembler: megahit, spades, mitoassemble
cd 06_mitoz
mitoz all --outprefix {prefix of output file, ex.japanese_eel_mito} \
          --thread_number 128 --workdir {wd} --clade {clade} \
          --genetic_code {genetic code} \
          --species_name {species name} \
          --fq1 {gzipped/ungzipped R1} --fq2 {gzipped/ungzipped R2} \
          --insert_size {WGS insert len} --requiring_taxa {clade} \ 
          --memory 200 --data_size_for_mt_assembly 0 \
          --assembler {megahit, spades, mitoassemble}
bio gff {output}.gbf > {output}.gff
```

MitoHiFi 3.0.0_cv1 docker
*Note: We recommend to use mitoZ to annotate with sequence again, since there might be some relative rare start codon issue that MitoHiFi cannot idnetify.*
Note: bio gff cannot digest the stranded info in gff, so we will correct it during fix and merge stage.
```
cd 06_mitohifi
nohup mitohifi.py -r {HiFi raw reads} -f {reference fasta} \
                  -g {reference .gbk file} -t 128 -d -a animal \
                  -o {genetic code} > mitohifi.log 2>&1 & 
bio gff final_mitogenome.gb > final_mitogenome.gff
```

If the gene range cross the split region, it is recommended to split on the safe places and annotate again
```
python {the path of this directory}/check_adjust_mito.py 
        -f 06_mitoz/{prefix of output file}.result/{prefix of output file}.{prefix of output file}.megahit.mitogenome.fa.result/{prefix of output file}_{prefix of output file}.megahit.mitogenome.fa_mitoscaf.fa.gbf.fasta \
        -g 06_mitoz/{out gff} -o 00_start_annotation/{file basename}_mito_adjusted.fa \
        -i 00_start_annotation/{file basename}_mito_adjusted.info

mitoz annotate --outprefix mitoz --thread_number {cpus} --fastafiles {file basename, cp adjusted fa to the 06_mitoz_adjusted dir}_mito_adjusted.fa --species_name {species name} --genetic_code {genetic code} --clade {clade}
```

### Check nuclear or nuclear+mitochondrial genome with FCS adaptor/gx
It requires some efforts to configure, users can follow the instruction recorded in the link.  
gx manual and example: https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart
gx report info: https://github.com/ncbi/fcs/wiki/FCS-GX-output
```
#Check contamination from adpator 
{the path of this directory}/run_fcsadaptor.sh --fasta-input 00_raw_genome/{genome fasta} \
                                               --output-dir 00_fcs_adaptor_mito_nuclear --euk

#Remove the adaptor
cat 00_raw_genome/{genome fasta} | python3 {the path of this directory}/fcs.py clean genome --action-report 00_fcs_adapter_mito_nuclear/fcs_adaptor_report.txt --output 00_fcs_adaptor_mito_nuclear/{genome fasta basename}_clean_adaptor.fasta --contam-fasta-out 00_fcs_adaptor_mito_nuclear/contam.fasta

#Double check contamination from adpator after removing the adaptors 
{the path of this directory}/run_fcsadaptor.sh --fasta-input 00_fcs_adapter_mito_nuclear/{genome fasta basename}_clean_adaptor.fasta --output-dir 00_fcs_adapter_double_check_mito_nuclear --euk

#Check contamination from other species
python3 {the path of this directory}/fcs.py screen genome --fasta 00_fcs_adapter_mito_nuclear/{genome fasta basename}_clean_adaptor.fasta --out-dir 00_fcs_gx_mito_nuclear --tax-id {tax id} --gx-db ../fcs-db/gxdb

#replace all in fifth column to EXCLUDE or other tag (depends on user).
#genome fasta basename, ex. MTW_genome.fa -> MTW_genome
#tax id: hypothesized species or the related one

awk -F'\t' -v OFS='\t' '{ $5 = "EXCLUDE"; print }'  00_fcs_adaptor_mito_nuclear/{genome fasta basename}.{tax id}.fcs_gx_report.txt >  00_fcs_adapter_mito_nuclear/{genome fasta basename}.{tax id}.fcs_gx_report_modified.txt

#Clean genome
cat 00_fcs_adapter_mito_nuclear/{genome fasta basename}_clean_adaptor.fasta | python3 {the path of this directory}/fcs.py clean genome --action-report ./00_fcs_gx_mito_nuclear/{genome fasta basename}.{tax id}.fcs_gx_report_modified.txt --output 00_fcs_gx_mito_nuclear/{genome fasta basename}_clean.fasta --contam-fasta-out 00_fcs_gx_mito_nuclear/contam.fasta

#Double check contamination from other species
python3 {the path of this directory}/fcs.py screen genome --fasta 00_fcs_gx_mito_nuclear/{genome fasta basename}_clean.fasta --out-dir 00_fcs_gx_double_check_mito_nuclear --tax-id {tax id from NCBI}
```
Note
*If FCS-adapter/gx detect the contamination, users can decide whether they remove the candidates or not. If users modify the sequence, please rememeber to check BUSCO and basic statistics again*

### Check basics number of genome assembly
```
python {the path of this directory}/assembly_quality_check_v5.py {fasta file} 
```
### Check with BUSCO genome mode
The miniprot is recommended to apply here. Please start a docker of BUSCO first and run the command within. It takes less an hour (threads/cpus = 128).
```
nohup busco -i {genome fasta} -m genome -l {ex.mollusca_odb10} -c {cpu numbers} -o {output dir basename} --miniprot --out_path {output directory} > run_busco_genome_miniprot.log 2>&1 &
```

Copy genome fastas (nuclear, mitochondria, nuclear+mitochondria) to 00_start_annotation.

# Repeat soft masking
Run repeat masking with RepeatModeler, RepeatMasker in TETools and second-time TRF (followed the instruction from Braker3 2024 [paper](https://arxiv.org/abs/2403.19416))
```
#Start a docker container with TETools
docker run --name tetools-v1.89.2 -it -d -u root -v /home/abieskawa/output:/output --workdir /output/ dfam/tetools:1.89.2 bash

docker exec -it tetools-v1.89.2 bash

#Repeat masking with RepeatMasker and RepeatModeler
nohup bash -c "time {the path of this directory}/repeatmask.sh 
           -i /output/{file basename}/00_start_annotation/{genome fasta basename}_clean.fasta{or genome fasta in 00_raw_genome if there is no mitochondria or contamination issue in the raw genome} \
           -o /output/{file basename}/ -t 128 \
           -l {ex./output/RepeatMasker_DB/famdb} -s '{search target}' \
           -f '{file basename}' --log repeatmask.log 2>&1" 2> run_repeatmask_time.log > /dev/null &

#Run repeat masking again with TRF
nohup bash -c "time ({the path of this directory}/trf_mask.sh -i 03_SoftMask/{genome after repeat masking} -o 04_TRF_mask -t 128 > run_trf_mask.log 2>&1)" 2> run_trf_time.log > /dev/null &
```
The generated table can sometimes hard to interpret. As a result, users can summarize by themselves. Here is the example script to process the result.
```
#buildSummary.pl is already built inside the tetools docker.
python {the path of this directory}/generate_karyotype.py -n {longest x scaffold} {genome fasta} | awk '{print $3 "\t" $6}' | sed 's/  */\t/g' > RepeatCount.tsv
[within tetools docker] buildSummary.pl -genome RepeatCount.tsv {input genome}.out > repeat_report.txt
python {the path of this directory}/summarize_buildSummary.py 03_SoftMask/repeat_report.txt > 03_SoftMask/Repeatsummary.txt
```

# Structural Annotation 
## Run mitochondria Annotation with mitos if it is not done with MitoZ/MitoHiFi
If you have run mitoZ --all, there is no need to run annotation again.
Version: runmitos.py 2.1.9
```
runmitos.py -i {input mitochondria genome} -c {code for different species} -R ../MITOS2-refdata/ -r refseq89m -o 06_mitos/ --best --debug

#Check and adjust if any gene span the split point
python {the path of this directory}//check_adjust_mito.py -f 06_mitoz/{prefix of output file}.result/{prefix of output file}.japanese_eel_mito.megahit.mitogenome.fa.result/{prefix of output file}_{prefix of output file}.megahit.mitogenome.fa_mitoscaf.fa.gbf.fasta -g 06_mitoz/{out gff} -o 00_start_annotation/{file basename}_mito_adjusted.fa -i 00_start_annotation/{file basename}_mito_adjusted.info

#Then, run mitos annotation again
```

## Run ncRNA Annotation
```
#cmscan
awk '!/^>/ {sum += length($0)} END {print sum}' {genome fasta} -> genome size
{genome size} * 2 / 10**6 = genome size param
nohup bash -c "time cmscan --cpu 128 -Z {genome size param} --cut_ga --rfam --
nohmmonly --tblout 06_ncRNA/{file basename}.tblout --fmt 2 --clanin ../Rfam_db/Rfam.clanin --oskip ../Rfam_db/Rfam.cm 00_start_annotation/{genome fasta sequence} > 06_ncRNA/{file basename}.cmscan 2> run_cmscan.log" > run_cmscan_time.log 2>&1 &

#tRNAscan-se
nohup bash -c "time tRNAscan-SE -E -I --detail -H --gff 06_ncRNA/{file basename}_tRNAscan.gff --stats 06_ncRNA/{file basename}.stats -d --thread 128 00_annotation_start_genome/{genome fasta sequence} -o 06_ncRNA/{file basename}_tRNAscan.out -f 06_ncRNA/{file basename}_tRNAscan.struct> run_{outfile name}_tRNA.log 2>&1" > run_{file basename}_tRNAscan_time.log 2>&1 &

#(Optional)Filter with EukHighConfidenceFilter
#It may require species specific criteria
EukHighConfidenceFilter --result 06_ncRNA/{file basename}_tRNAscan.out --ss 06_ncRNA/{file basename}_tRNAscan.struct --output 06_ncRNA/ --prefix {file basename}_filter 

#Filter gff file with EukHighConfidenceFilter output file
python {the path of this directory}/trna_confidence_filter.py -i 06_ncRNA/{file basename}_tRNAscan.gff -c 06_ncRNA/{file basename}_filter.out -o {file basename}_tRNAscan_filter -d 06_ncRNA/ -v > 06_ncRNA/trna_confidence_filter.log

#Insert intron to the tRNAscan-se gff file
python {the path of this directory}/add_intron_trnascanse.py -i 06_ncRNA/{file basename}_tRNAscan.gff -o 06_ncRNA/{file basename}_tRNAscan_add_intron.gff
```

## Prepare protein evidence:
The script was modified with some informatiuon provided in https://github.com/gatech-genemark/ProtHint/issues/39
Users should download odb12v0_level2species.tab.gz, odb12v0_levels.tab.gz, odb12v0_aa_fasta.gz, odb12v0_species.tab.gz can be downloaded from OrthoDB v12 website first.
```
nohup python {the path of this directory}/extract_clade_proteins.py --clade {Mollusca} --levels odb12v0_levels.tab.gz --level2species odb12v0_level2species.tab.gz --fasta odb12v0_aa_fasta.gz --output Mollusca.fa --report-species --report-species --speciesmap odb12v0_species.tab.gz > extract_mollusca.log 2>&1 &
```
## Run fasta simplification
Without this steps, Braker3 script is still executable, but it will complain annoying messages. 
```
perl {the path of this directory}/simplifyFastaHeaders.pl {combined protein evidence file name} prot {simplified header protein evidence file name} {simplified headers file file name}
``` 

## Prepare RNA evidence
### Illumina RNA-seq
```
#Download from SRA DB
nohup ~/tools/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump --split-files {SRA run ID} -O {outdir} --threads 128 > 05_raw_RNA_PacBio/download.log 2>&1 &

#Purge adapters from illumina RNA-seq, and those being too short or having a quality score below 20 from RNA will be removed.
nohup {the path of this directory}/run_cutadapt.sh -d {raw RNA-seq dir} -o {cleaned RNA-seq dir} -t {cpus} -l 5 -q 20 -Q 20 -f sample_trim.txt > run_cutadapt.log 2>&1 &

#Run RNA mapping with STAR and keep only reads flagged as 2.
nohup python {the path of this directory}/run_rna_mapping.py --mode short-read --genomepath {nuclear genome} --genomedir {genome index dir} --wd {dir for cleaned.fq/fastq} --out_dir {outdir} --threads {cpus} > run_rna_mapping.log 2>&1 &
```
### PacBio Iso-seq
In this step, cutadapt -g/-a/--poly-a was applied first to extract cDNA with 5'/3' primer and poly-A tails, and then those with 5'/3' primer inside the reads are discarded. Next, fastq file was transformed to unaligned bam with picard FastqToSam, then isoseq cluster (n*log(n) with qv guided info) was applied to generate longer reads. The short read data can be optionally polished isoseq data using the lordec-correct procedure to improve accuracy.

Map the reads to genome with minmap2 as team braker2 recommend

```
nohup python {the path of this directory}/isoseq_preprocess.py -i "A.fastq:B_1.fastq,B_2.fastq" -i -g "ADAPTER1{30}" -a ADAPTER3 

nohup python {the path of this directory}/run_rna_mapping.py --mode iso-seq --genomepath {nuclear genome} --wd {dir for cleaned.fq/fastq} --out_dir {outdir} --threads {cpus}  >run_rna_mapping_isoseq.log 2>&1 &
```
 
## Run Braker3
User can use TSEBRA to combine the prediction from isoseq and short-read.
```
#Run with detached docker container
---> Result with Isoseq evidence can combine with teh result from short-read seq  
docker run -d --rm -u root -v /home/abieskawa/output:/output --workdir /output/{file basenam} teambraker/braker3:tag(it should be isoseq if RNA data is isoseq) {the path of this directory}/run_braker3.sh -t 120 -g {TRF masked genome} -p {protein evidence} -s {parameter set name} -w {outdir} -b {bam file dir} -l {BUSCO odb}

#Run within docker
nohup bash -c "time ({the path of this directory}/run_braker3.sh -t 120 -g {TRF masked genome} -p {protein evidence} -s {parameter set name} -w {outdir} -b {bam file dir} -l {BUSCO odb}  > run_braker3.log 2>&1)" 2> run_braker3_time.log > /dev/null &
```

## Check Braker3 result
In this script, I include calculate the total number of gene, mRNA, exon, protein, single-exon, and median length of gene, gff_QC check, BUSCO, OMArk check. Since the gff file version may be different in newer version of braker3, please check my source code and results before and after applying it. 
In this stage, our target is to check basic statistics of nuclear protein-coding genes.
```
#-p: please input nuclear protein (braker has braker.aa output)
#Note: BUSCO does not have non-nuclear database, so do not combine mitoz annotation result and check BUSCO value. 
nohup python {the path of this directory}/eva_annotation.py structural -g 06_braker3/braker.gff3 -f {genome fasta} -b {file basename} -o 06_braker3/ -p {braker protein sequence}  --busco_lineage {ex.actinopterygii_odb10} -t 128 --busco_out_path BUSCO/ --busco_docker_image ezlabgva/busco:v5.8.2_cv1 --busco_workdir /home/abieskawa/output/{file basename}/ --omark_dir omark_omamer_output > run_eva_braker.log 2>&1 &
```

## Run gtf2gff
Remember to make everything 777 or any similar allow BUSCO docker to access the data.
```
sudo chmod 777 {wd}
perl {the path of this directory}/gtf2gff.pl <braker.gtf --out=braker.gff --printExon --printUTR --printIntron --gff3
```

## Preprocess the output before running functional annotation
```
#mitos
python {the path of this directory}/fix_merge_sort_gff.py -I infernal 06_ncRNA/{file basename}_ncRNA.tblout -I trnascan-se 06_ncRNA/{file basename}_tRNAscan_add_intron.gff -I braker3 {braker3 gff ex.06_braker3/braker.gff} -I mitos 06_mitos/result.gff --fmt2 --ignore-trna --basename {ex.T_jap} --source cmscan --outdir 07_merge_gff/

#mitohifi/mitoz
python {the path of this directory}/fix_merge_sort_gff.py -I infernal 06_ncRNA/{file basename}_ncRNA.tblout -I
trnascan-se 06_ncRNA/{file basename}_tRNAscan_add_intron.gff -I braker3 {braker3 gff ex.06_braker3/braker.gff} -I gb 06_mitoz_adjusted/{file basename}_mito_adjusted.gff --fmt2 --ignore-trna --basename {ex.O_nil} --source cmscan --outdir 07_merge_gff/

#Check format of merged file
nohup python {the path of this directory}/eva_annotation.py structural -g 07_merge_gff/merged_fix.gff -f {genome fasta} -b {file basename} -o 07_merge_gff/ --split_scaffolds {the scaffold desired to be ccalculated separately, ex. mitochondria} > run_eva_combined.log 2>&1 &

python {the path of this directory}/molas_fasta-intron-lowercase.py -g 00_annotation_start_genome/MTW_genome_clean_nuclear_mito_adjusted.fasta -r 07_merge_gff/merged_fix.gff -p M_tai -a combined -v --pep_mitos 06_mitos_adjusted/result.faa --pep_braker 06_braker3/braker.aa --nameprefix M_tai -d 08_FASTA

#mitohifi/mitoz
python {the path of this directory}/molas_fasta-intron-lowercase.py -g 00_annotation_start_genome/{genome fasta} -r 07_merge_gff/merged_fix.gff -p {prefix. ex.M_tai} -a combined -v --nmitopep --table 5 --pep_braker 06_braker3/braker.aa --nameprefix {prefix. ex.M_tai} -d 08_FASTA
```

## Transform gff to gtf3
```
agat levels --expose
agat_convert_sp_gff2gtf.pl --gff merged_fix.gff -o merged_fix.gtf --gtf_version 3
-----
#Insert the missing LEVEL 2 entries manually, and make sure name of entries are all lowercase 
-----
# Keep in mind that the transformation will sjip tRNA pseudogene and biological region
agat_convert_sp_gff2gtf.pl --gff merged_fix.gff -o merged_fix.gtf --gtf_version 3
```

## Prepare NR protein database
### Downolad NR DB
```
nohup update_blastdb.pl --source "ncbi" --decompress nr --num_threads 0 > download_nr.log 2>&1 &
``` 
### Extract fasta with blastdbcmd
Run the command line inside the directory nr
```
blastdbcmd -db nr -entry all -out nr000-122.fasta
```
### Prepare the database with diamond makedb
Please refer to the manual described here: https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#makedb-options
```
nohup diamond makedb --in nr/nr.000-125.fasta --threads 128 -d diamond_nr/nr.000-125.dmnd --taxonmap nr
/prot.accession2taxid.FULL.gz --taxonnodes nr/nodes.dmp --taxonnames nr/names.dmp > diamond_nr/diamond_makedb.log 2>&1 &
```
### Run Diamond blastp NR
```
nohup bash -c "time (diamond blastp --header simple --max-target-seqs 1 --outfmt 6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids -q 08_FASTA/{protein fasta} -d ~/output/diamond_nr/nr.000-125.dmnd -o 11_diamond_blastp_mito_nuclear/{output prefix}_diamond.tsv --threads 128 > run_diamond_blastp.log 2>&1)" 2> run_diamond_blastp_time.log > /dev/null &

python {the path of this directory}/eva_annotation.py diamond-nr -p {input protein sequence} -d 11_diamond_blastp{_mito}_nuclear/{output prefix}_diamond.tsv -o 11_diamond_blastp{_mito}_nuclear/{output prefix}
```

## Run DeepTMHMM
Make sure the server has NVIDIA GPU. It can use cpu though, but it will take longer time to finish the job.
```
nohup python ../run_deeptmhmm.py 08_FASTA/{protein sequence fasta} 09_deeptmhmm_mito_nuclear/ {file basename} > run_deeptmhmm.log 2>&1 &
```

## Download KAAS result
Note that it can take days(~2 days) to download, please stay calm and be patient.
*plesae check the download.log, sometimes there are error, and it will give up downloading that entry!*
*Besises, the KAAS file structure can be modified during the update of the DB, this script will be updated over time.* 
```
nohup python3 {the path of this directory}/downloadkaas.py {"download link from kaas"} > download.log 2>&1 &

python {the path of this directory}/eva_annotation.py kaas 10_kaas{_mito}_nuclear/query.ko 10_kaas{_mito}_nuclear
```

## Run InterProScan
We used v5.73-104.0, so the script eva_annotation.py and interproscan_extract might require modification if the version is not this one.
Sometimes, it might require to adjust the memory limitation, be careful to any error messages.
```
nohup bash -c 'time ( _JAVA_OPTIONS="-Xmx1536g" interproscan.sh -i {input peptide seq. ex.08_FASTA/marine_tilapia_MTW02_combined_pep.fa} -goterms -f tsv --output-file-base 12_interproscan_mito_nulcear/interproscan_result -cpu 128 > run_interproscan.log 2>&1 )' 2>run_interproscan_time.log > /dev/null &

python {the path of this directory}/eva_annotation.py interpro -p {input peptide seq.} -i 12_interproscan_mito_nulcear/interproscan_result.tsv -o 12_interproscan_mito_nulcear/{file basename}

#Extract the info that MOLAS ask for
python {the path of this directory}/interproscan_extract.py -in interproscan_result.tsv -p {file basename}
``` 
## Venn diagram of Annotated Gene
``` 
[within R conda env]{the path of this directory}/vennplot.r -i GO 10_interproscan_mito_nulcear/{file prefix}_mito_nuclear_annotated_gene_go.tsv -i IPR 10_interproscan_mito_nulcear/{file prefix}_mito_nuclear_annotated_gene_ipr.tsv -i diamond 11_diamond_blastp_mito_nulcear/{file prefix}_mito_nuclear.informative_gene_list.tsv -i KAAS 12_kaas_mito_nuclear/kaas_annotated_genes.tsv -a {file prefix}_MOLAS_input/FASTA/{file prefix}_combined_pep.fa -o functional_venn.png
``` 
## Run RNA analysis
This process can generate table for gene/transcript level for further analysis
#Although STAR can handle gff, it requires to be modified the normal gff format further
```
nohup python {the path of this directory}/run_rna_analysis.py --reads-dir {short-read cleaned seq dir} --genome-fasta {nuclear(+mito) genome fasta} --annotation {merged.gtf} --out-dir RNAseq_ana --threads 128 --transcript-fpkm --transcript-counts --gene-fpkm --gene-tpm --gene-counts --merge-expression > RNAseq_ana.log 2>&1 &

```

## Upload to MOLAS genome browser
``` 
python {the path of this directory}/molas_genome_browser.py -g {genome fasta} -r {gff} -p {file prefix} -n {diamond tsv}
``` 

## Detect centromere, gap and telomere


## Draw Circos plot
```
python {the path of this directory}/gff2circos.py --input ../07_merge_gff/merg
ed_fix.gff --out gene_w100000.hist.txt --window 100000 --biotype protein_coding --genome_fasta ../04_SplitMitoNuclear/{file basename}_genome_1_nuclear.fasta --seqids circos_target_seq.txt

python scripts/repeatout2circos.py {RepeatMasker Output, ex.beltfish_genome_1.fasta.out} 17_circos/repeat_dna_w100000.hist.txt --targets {a one-column txt} --window 100000 --pattern DNA --fasta {genome fasta}

docker run --rm -v $(pwd):/data alexcoppe/circos -conf /data/circos.conf -outputdir /data
```

## Citation and the tools used in this pipeline
### FCS scripts
fcs.py/run_fcsadaptor.sh [[link](https://github.com/ncbi/fcs)]

### Repeat Annotation/masking
TETools (v1.89.2) [[link](https://github.com/Dfam-consortium/TETools)]
TRF (v4.09) [[link](https://github.com/Benson-Genomics-Lab/TRF)]
bedtools (v2.31.1) [[link](https://github.com/arq5x/bedtools2)]
splitMfasta.pl [[link](https://github.com/Gaius-Augustus/Augustus/blob/487b12b40ec3b4940b6b07b72bbb443f011f1865/scripts/splitMfasta.pl)]
parseTrfOutput.py [[link](https://github.com/gatech-genemark/BRAKER2-exp/blob/34e9d1dfd7228128968063f76b37d29c73a39efc/bin/trf-scripts/parseTrfOutput.py)]

### protein preprocessing
simplifyFastaHeaders.pl[[link](https://github.com/Gaius-Augustus/Augustus/blob/487b12b40ec3b4940b6b07b72bbb443f011f1865/scripts/simplifyFastaHeaders.pl#L7)]

### Braker (gtf to gff)
gtf2gff.pl[[link](https://github.com/Gaius-Augustus/Augustus/blob/487b12b40ec3b4940b6b07b72bbb443f011f1865/scripts/gtf2gff.pl#L347)]
genome_anno.py[[link](https://github.com/Gaius-Augustus/Tiberius/blob/bfa9b37eaeca0794dd2b508c32e3f59bd28ec479/bin/genome_anno.py#L5)]

### Get aa seq from Tsebra
Fix_Augustus_gtf.pl[[link](https://github.com/Gaius-Augustus/BRAKER/issues/457#issuecomment-1050475171)]

### RNAseq analysis
stringtie (v2.2.3), prepDE.py[[link](https://github.com/gpertea/stringtie)] Note: prepDE.py3 in that repository was used here.
STAR (v2.7.11b)[[link](https://github.com/alexdobin/STAR)]