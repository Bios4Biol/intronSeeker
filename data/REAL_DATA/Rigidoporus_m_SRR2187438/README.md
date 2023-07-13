Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GD/MN/GDMN01/GDMN01.1.fsa_nt.gz
gzip -d GDMN01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/008/SRR2187438/SRR2187438_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/008/SRR2187438/SRR2187438_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GDMN01.1.fsa_nt -1 SRR2187438_1.fastq -2 SRR2187438_2.fastq --prefix GDMN01 -o GDMN01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GDMN01/hisat2_GDMN01.sort.bam -r GDMN01.1.fsa_nt --prefix GDMN01 --output splitReadSearch_GDMN01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GDMN01.1.fsa_nt
r1:SRR2187438_1.fastq
r2:SRR2187438_2.fastq
flagstat:hisat2_SRR2187438.sort.flagstat.txt
candidat:srs_SRR2187438_candidates.txt
split:srs_SRR2187438_split_alignments.txt
prefix:GDMN01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR2187438.cfg;

```

HTML report is available in this directory.
