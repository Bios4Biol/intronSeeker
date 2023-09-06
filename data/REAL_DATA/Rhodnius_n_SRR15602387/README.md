Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/JI/GJJI01/GJJI01.1.fsa_nt.gz
gzip -d GJJI01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR156/087/SRR15602387/SRR15602387_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR156/087/SRR15602387/SRR15602387_1.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJJI01.1.fsa_nt -1 SRR15602387_1.fastq.gz -2 SRR15602387_2.fastq.gz --prefix GJJI01 -o GJJI01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJJI01/hisat2_GJJI01.sort.bam -r GJJI01.1.fsa_nt --prefix GJJI01 --output splitReadSearch_GJJI01
```

### Step 3: run simulation2HTML

Configuration file:
```diff
nano SRR15602387.cfg
```


```diff
[Defaults]
fasta:GJJI01.1.fsa_nt
r1:SRR15602387_1.fastq.gz
r2:SRR15602387_2.fastq.gz
flagstat:GJJI01/hisat2_GJJI01.sort.flagstat.txt
candidat:splitReadSearch_GJJI01/srs_GJJI01_candidates.txt
split:splitReadSearch_GJJI01/srs_GJJI01_split_alignments.txt
prefix:GJJI01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR15602387.cfg;

```

HTML report is available in public directory and here https://emilien.lasguignes.pages.mia.inra.fr/intronSeeker/report_Rhodnius_n__SRR15602387.html

