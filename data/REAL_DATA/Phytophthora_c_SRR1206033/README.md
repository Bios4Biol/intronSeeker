Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/GX/GBGX01/GBGX01.1.fsa_nt.gz
gzip -d GBGX01.1.fsa_nt.gz
```

### Single reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/003/SRR1206033/SRR1206033.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GBGX01.1.fsa_nt -1 SRR1206033.fastq.gz  --prefix GBGX01 -o GBGX01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GBGX01/hisat2_GBGX01.sort.bam -r GBGX01.1.fsa_nt --prefix GBGX01 --output splitReadSearch_GBGX01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano SRR1206033.cfg
```

```diff
[Defaults]
fasta:GBGX01.1.fsa_nt
r1:SRR1206033.fastq.gz
flagstat:GBGX01/hisat2_GBGX01.sort.flagstat.txt
candidat:splitReadSearch_GBGX01/srs_GBGX01_candidates.txt
split:splitReadSearch_GBGX01/srs_GBGX01_split_alignments.txt
prefix:GBGX01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR1206033.cfg;

```

HTML report is available in public directory and here http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Phytophtora_c_GBGX01.html
