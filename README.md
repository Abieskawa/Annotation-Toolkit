# Annotation-Toolkit
This repository is the home place for any tool I used to annotated the genome, including the command line, evaluation and download tool

## Run fasta simplification
Without this steps, Braker3 script is still executable, but it will complain annoying messages. 
```
perl ../Annotation-Toolkit/simplifyFastaHeaders.pl braker_odbv12v0_Metazoa_UniProtKB_6447_ncbi_6580.fa prot braker_odbv12v0_Metazoa_UniProtKB_6447_ncbi_6580_s.fa simplified_asian_hard_clam.headers
``` 
## Run Braker3


## Run gtf2gff
```
perl ~/tools/braker_tools/gtf2gff.pl <braker.gtf --out=braker.gff --printExon --printUTR --gff3
```

## Run BUSCO proteome
```
busco -i 06_braker3/braker.aa -m proteins -l mollusca_odb10 -c 128 -o asian_hard_clam_v1 --out_path BUSCO/ > run_busco_v1.log
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
```
diamond makedb --in nr/nr000-122.fasta --threads 128 -d diamond_nr
```

## Run DeepTMHMM
Make sure the server has NVIDIA GPU. 
```
nohup python ../run_deeptmhmm.py 07_braker3_v7_renamed/braker_renamed_wout_asterisk.aa 08_deeptmhmm/ marine_tilapia > run_deeptmhmm.log 2>&1 &
```

## Download KAAS result
Note that it can take days(~2) to download, please stay calm and be patient.
```
nohup python3 ../../Annotation-Toolkit/downloadkaas.py "https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1739414180&key=9Ze6Wzzp" > download.log 2>&1 &
```

## Run InterProScan
```
nohup bash -c "time (interproscan.sh -i ../07_braker3_v7_renamed/braker_renamed_removestar.aa -goterms -f tsv --output-file-base interproscan_result -cpu 128 > ../run_interproscan.log 2>&1)" 2>../run_interproscan_time.log > /dev/null & 
``` 
