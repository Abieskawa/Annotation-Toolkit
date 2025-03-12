# Annotation-Toolkit
This repository is the home place and notes for any tool applied to annotated the genome, including the command line, evaluation and download tool. Although NCBI can do annotation for the users, they may take some time to wait. Also, their annotation tools does not make it open source until now, so a toolkit for annotation can allow the users to adjust according to their requirement. This toolkit includes structural annotation for protein-coding gene, ncRNA, mitochondria gene. Codes for functional annotation for protein-coding gene is recorded here, too. We recommend read and run everything with your genome first to know every details, and any adjustment might require to your genome, if you are familar with the workflow, then switch to the nextflow and docker.

This script has been tested on beltfish, marine tialpia, asian hard clam.

# Check the received genome quality
Check the assembly information, BUSCO, and FCS-gx/adapter (and clean if it is necessary)
### Check basics number of genome assembly
```
python assembly_quality_check_v5.py {fasta file} 
```
### Check with BUSCO genome mode
The miniprot is recommended to apply here. Please start a docker of BUSCO first and run the command below. It takes less an hour (threads = 128).
```
nohup busco -i {genome fasta} -m genome -l {ex.mollusca_odb10} -c {cpu numbers} -o asian_hard_clam_clean_genome_miniprot --miniprot{metaeuk} --out_path {output directory} > run_busco_genome_miniprotlog 2>&1 &
```
### Check with FCS adapter/gx
It requires some efforts to configure, users can follow the instruction recorded in the link.  
gx manual and example: https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart
gx report info: https://github.com/ncbi/fcs/wiki/FCS-GX-output
```
python3 {the path of this directory}/fcs.py screen genome --fasta {genome fasta} --out-dir {output directory} --gx-db "$GXDB_LOC/gxdb" --tax-id {tax id}

#replace all in fifth column to EXCLUDE (depends on user).
#genome fasta basename, ex. MTW_genome.fa -> MTW_genome
#tax id: hypothesized species or the related one

awk -F'\t' -v OFS='\t' '{ $5 = "EXCLUDE"; print }' {output directory}/{genome fasta basename}.{tax id}.fcs_gx_report.txt > {output directory}/{genome fasta basename}.{tax id}.fcs_gx_report_modified.txt

#Clean genome
cat MTW_genome.fasta | python3 {the path of this directory}/fcs.py clean genome --action-report ./{output directory}/{genome fasta basename}.{tax id}.fcs_gx_report_modified.txt --output {output clean directory}/{genome fasta basename}_clean.fasta --contam-fasta-out {output clean directory}/contam.fasta
```

## (Optional) Check which scaffold is mitochondria genome
```
mitoz findmitoscaf --fastafile {genome, gz supported here} --outprefix {prefix} --workdir {dir} --thread_number {threads} --clade {Chordata,Arthropoda,Echinodermata,Annelida-segmented-worms,Bryozoa,Mollusca,Nematoda,Nemertea-ribbon-worms,Porifera-sponges}
```

## (Optional) If one wants to extract mitochondria and nucleus chromosomes (Longest)
```
python {the path of this directory}/extract_nucleus_query_genome.py -i {the genome after cleaning with fcs} -on {nuclear fasta name} -os {fasta name of mitochondria or chloreplast} -n {the longest n seqs} -s {mitochondria,chloroplast,or the name of the target genome}
``` 
# Repeat soft masking
Run repeat masking with RepeatModler, RepeatMasker and second-time TRF (followed the instruction from Braker3 2024 paper)
```
#Repeat masking with RepeatMasker and RepeatModeler
nohup bash -c "time (../Annotation-Toolkit/repeatmask.sh -i 00_keep_nuclear_mito_genome/MTW_genome_keep_nuclear_mito.fasta -o /output/asian_hard_clam -t 128 -l /output/Lib_fish/famdb -s 'Bivalvia' -f 'asian_hard_clam' --log repeatmask.log 2>&1)" 2> run_repeatmask_time.log > /dev/null &

#Run repeat masking again with TRF
nohup bash -c "time (../Annotation-Toolkit/trf_mask.sh -i 03_SoftMask/MTW_genome_keep_nuclear_mito.fasta.masked -o 04_TRF_mask -t 128 > run_trf_mask.log 2>&1)" 2> run_trf_time.log > /dev/null &
```
# Structural Annotation 
## Run mitochondria Annotation with mitos
Version: runmitos.py 2.1.9
```
runmitos.py -i 00_annotation_start_genome/MTW_genome_clean_mitochondria.fasta -c 5 -R ../MITOS2-refdata/ -r ../MITOS2-refdata/refseq89m -o 06_mitos/ --debug
```

## Run ncRNA Annotation
```
nohup bash -c "time cmscan --cpu 128 -Z 2081.464022 --cut_ga --rfam --
nohmmonly --tblout 06_ncRNA/asian_hard_clam_ncRNA.tblout --fmt 2 --clanin ~/output/Rfam_db/Rfam.clanin --oskip ~/output/Rfam_db/Rfam.cm 00_annotation_start_genome/MTW_genome_clean_nuclear.fasta > 06_ncRNA/asian_hard_clam_ncRNA.cmscan 2> run_cmscan.log" > run_cmscan_time.log 2>&1 &

nohup bash -c "time tRNAscan-SE -E -I --gff 06_ncRNA/asian_hard_clam_tRNAscan.gff --stats 06_ncRNA/asian_hard_clam_tRNAscan.stats -d --thread 128 00_annotation_start_genome/MTW_genome_clean_nuclear.fasta > run_asian_hard_clam_tRNA.log 2>&1" > run_asian_hard_clam_tRNAscan_time.log 2>&1 &
```

## Prepare protein evidence:
The script was modified with some informatiuon provided in https://github.com/gatech-genemark/ProtHint/issues/39
Users should download odb12v0_level2species.tab.gz, odb12v0_levels.tab.gz, odb12v0_aa_fasta.gz, odb12v0_species.tab.gz can be downloaded from OrthoDB v12 website first.
```
nohup python ../Annotation-Toolkit/extract_clade_proteins.py --clade Mollusca --levels odb12v0_levels.tab.gz --level2species odb12v0_level2species.tab.gz --fasta odb12v0_aa_fasta.gz --output Mollusca.fa --report-species --report-species --speciesmap odb12v0_species.tab.gz > extract_mollusca.log 2>&1 &
```
## Run fasta simplification
Without this steps, Braker3 script is still executable, but it will complain annoying messages. 
```
perl ../Annotation-Toolkit/simplifyFastaHeaders.pl braker_odbv12v0_Metazoa_UniProtKB_6447_ncbi_6580.fa prot braker_odbv12v0_Metazoa_UniProtKB_6447_ncbi_6580_s.fa simplified_asian_hard_clam.headers
``` 

## Prepare RNA evidence
```
#Purge adapters from sequences, and those being too short or having a quality score below 20 from RNA.
nohup ../Annotation-Toolkit/run_cutadapt.sh -d 05_raw_RNA_MTW/ -o 05_cleaned_RNA_MTW/ -t 128 -l 5 -q 20 -Q 20 -f sample_trim.txt > run_cutadapt.log 2>&1 &

#Run RNA mapping with STAR and keep only reads flagged as 2.
nohup python ../Annotation-Toolkit/run_rna_mapping_star_2pass.py --genomepath 00_keep_nuclear_mito_genome/MTW_genome_keep_nuclear_mito.fasta --genomedir star_genome_index/ --wd 05_raw_RNA_MTW/ --out_dir 05_RNA_mapping_star_2pass_MTW/ --threads 128 > run_rna_mapping.log 2>&1 &
```
 
## Run Braker3
```
nohup bash -c "time (../Annotation-Toolkit/run_braker3.sh -t 120 -g 04_TRF_mask/MTW_geno
me_clean_nuclear.combined.masked -p ../braker3_protein_evidence/odbv12v0_Mollusca_UniProtKB_6447_ncbi_6580_s.fa -s meretrix_taiwanica_v3 -w 
06_braker3_v3/ -b 05_RNA_mapping_star_2pass_MTW/ -l mollusca_odb10  > run_braker3_v3.log 2>&1)" 2> run_braker3_v3_time.log > /dev/null &
```

## Run gtf2gff
```
sudo chmod 777 {wd}
perl gtf2gff.pl <braker.gtf --out=braker.gff --printExon --printUTR --gff3
```

## Run evaluation of Braker annotation 
In this script, I include calculate the total number of gene, mRNA, exon, protein, single-exon, and median length of gene, gff_QC check, BUSCO, OMArk check. Since the regulation of gff file version is different version, please check my source code before applying my source code. 
```
nohup python ../Annotation-Toolkit/eva_annotation.py psa -g 06_braker3_v3/braker.gff -f 00_annotation_start_genome/MTW_genome_clean_nuclear.fasta -b asian_hard_clam_v3 -o 06_braker3_v3 -p 06_braker3_v3/braker.aa --busco_lineage mollusca_odb10 -t 128 --busco_out_path BUSCO/ --busco_docker_image ezlabgva/busco:v5.8.2_cv1 --busco_workdir /home/abieskawa/output/asian_hard_clam/ --omark_dir omark_omamer_output > run_eva_braker3_v3.log 2>&1 &
```
## Prepare protein before running Functional annotation
Rename the output with specified prefix (prefix_{entry name}) and rename protein sequence and remove its stop codon.   
```
python ../Annotation-Toolkit/rename_gtf.py --gtf 06_braker3_v3/braker.gtf --prefix M_tai --out 07_renamed_braker3_v3/braker_renamed.gtf

perl ../Annotation-Toolkit/gtf2gff.pl <braker_renamed.gtf --out=braker_renamed.gff --printExon --printUTR --gff3

python ../../Annotation-Toolkit/rename_protein_cds.py ../06_braker3_v3/braker.aa braker_renamed.aa M_tai

python ../../Annotation-Toolkit/remove_protein_star.py -I braker_renamed.aa -O braker_renamed_wout_asterisk.aa -L remove_star.log
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
nohup bash -c "time (diamond blastp --header simple --max-target-seqs 1 --outfmt 6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids -q ../07_braker3_v7_renamed/braker_renamed_wout_asterisk.aa -d ~/output/diamond_nr/nr.000-125.dmnd -o marine_tilapia_diamond.tsv --threads 128 > run_diamond_blastp.log 2>&1)" 2> time_output.log > /dev/null &
```

## Run DeepTMHMM
Make sure the server has NVIDIA GPU. It can use cpu though, but it will take longer time to finish the job.
```
nohup python ../run_deeptmhmm.py 07_braker3_v7_renamed/braker_renamed_wout_asterisk.aa 08_deeptmhmm/ marine_tilapia > run_deeptmhmm.log 2>&1 &
```

## Download KAAS result
Note that it can take days(~2) to download, please stay calm and be patient.
```
nohup python3 ../../Annotation-Toolkit/downloadkaas.py "https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1739414180&key=9Ze6Wzzp" > download.log 2>&1 &
```

## Run InterProScan
We used v5.73-104.0, so the script eva_annotation.py and interproscan_extract might require modification if the version is not this one.
Sometimes, it might require to adjust the memory limitation, be careful to any error messages.
```
nohup bash -c 'time ( _JAVA_OPTIONS="-Xmx1536g" interproscan.sh -i 07_renamed_braker3_v3/braker_renamed_wout_asterisk.aa -goterms -f tsv --output-file-base 11_interproscan/interproscan_result -cpu 128 > run_interproscan.log 2>&1 )' 2>run_interproscan_time.log > /dev/null &

python ../Annotation-Toolkit/eva_annotation.py interpro -p 07_renamed_braker3_v3/braker_renamed_wout_asterisk.aa -i 11_interproscan/interproscan_result.tsv -o 11_interproscan/asian_hard_clam

#Extract the info MOLAS ask for
python ../../Annotation-Toolkit/interproscan_extract.py -in interproscan_result.tsv -p marine_tilapia
``` 
