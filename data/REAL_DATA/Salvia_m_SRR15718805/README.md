Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/JN/GJJN01/GJJN01.1.fsa_nt.gz
gzip -d GJJN01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR15718805/SRR15718805_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR157/005/SRR15718805/SRR15718805_1.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJJN01.1.fsa_nt -1 SRR15718805_1.fastq.gz -2 SRR15718805_2.fastq.gz --prefix GJJN01 -o GJJN01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJJN01/hisat2_GJJN01.sort.bam -r GJJN01.1.fsa_nt --prefix GJJN01 --output splitReadSearch_GJJN01
```

### Step 3: run simulation2HTML

Configuration file:
```diff
nano SRR15718805.cfg
```

```diff
[Defaults]
fasta:GJJN01.1.fsa_nt
r1:SRR15718805_1.fastq.gz
r2:SRR15718805_2.fastq.gz
flagstat:GJJN01/hisat2_GJJN01.sort.flagstat.txt
candidat:splitReadSearch_GJJN01/srs_GJJN01_candidates.txt
split:splitReadSearch_GJJN01/srs_GJJN01_split_alignments.txt
prefix:GJJN01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR15718805.cfg;

```

HTML report is available in public directory and here http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Salvia_m_GJJN01.html
