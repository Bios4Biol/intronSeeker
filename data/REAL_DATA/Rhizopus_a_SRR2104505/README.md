Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GD/UK/GDUK01/GDUK01.1.fsa_nt.gz
gzip -d GDUK01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/005/SRR2104505/SRR2104505_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR210/005/SRR2104505/SRR2104505_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GDUK01.1.fsa_nt -1 SRR2104505_1.fastq.gz -2 SRR2104505_2.fastq.gz --prefix GDUK01 -o GDUK01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GDUK01/hisat2_GDUK01.sort.bam -r GDUK01.1.fsa_nt --prefix GDUK01 --output splitReadSearch_GDUK01
```

### Step 3: run simulation2HTML

Configuration file:

```dif
nano SRR2104505.cfg
```

```diff
[Defaults]
fasta:GDUK01.1.fsa_nt
r1:SRR2104505_1.fastq.gz
r2:SRR2104505_2.fastq.gz
flagstat:GDUK01/hisat2_GDUK01.sort.flagstat.txt
candidat:splitReadSearch_GDUK01/srs_GDUK01_candidates.txt
split:splitReadSearch_GDUK01/srs_GDUK01_split_alignments.txt
prefix:GDUK01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR2104505.cfg;

```

HTML report is available in this directory.
