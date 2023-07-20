Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GA/PK/GAPK01/GAPK01.1.fsa_nt.gz
gzip -d GAPK01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/004/SRR1023744/SRR1023744_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/004/SRR1023744/SRR1023744_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GAPK01.1.fsa_nt -1 SRR1023744_1.fastq.gz -2 SRR1023744_2.fastq.gz --prefix GAPK01 -o GAPK01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GAPK01/hisat2_GAPK01.sort.bam -r GAPK01.1.fsa_nt --prefix GAPK01 --output splitReadSearch_GAPK01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano SRR1023744.cfg
```

```diff
[Defaults]
fasta:GAPK01.1.fsa_nt
r1:SRR1023744_1.fastq.gz
r2:SRR1023744_2.fastq.gz
flagstat:GAPK01/hisat2_GAPK01.sort.flagstat.txt
candidat:splitReadSearch_GAPK01/srs_GAPK01_candidates.txt
split:splitReadSearch_GAPK01/srs_GAPK01_split_alignments.txt
prefix:GAPK01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR1023744.cfg;

```

HTML report is available in this directory.
