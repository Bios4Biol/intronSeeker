Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GH/FS/GHFS01/GHFS01.1.fsa_nt.gz
gzip -d GHFS01.1.fsa_nt.gz
```

### Paired reads:


```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR614/004/SRR6148374/SRR6148374_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR614/004/SRR6148374/SRR6148374_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GHFS01.1.fsa_nt -1 SRR6148374_1.fastq.gz -2 SRR6148374_2.fastq.gz --prefix GHFS01 -o GHFS01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GHFS01/hisat2_GHFS01.sort.bam -r GHFS01.1.fsa_nt --prefix GHFS01 --output splitReadSearch_GJJD01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano  SRR6148374.cfg
```


```diff
[Defaults]
fasta:GHFS01.1.fsa_nt
r1:SRR6148374_1.fastq.gz
r2:SRR6148374_2.fastq.gz
flagstat:GHFS01/hisat2_GHFS01.sort.flagstat.txt
candidat:splitReadSearch_GJJD01/srs_GHFS01_candidates.txt
split:splitReadSearch_GJJD01/srs_GHFS01_split_alignments.txt
prefix:GHFS01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  SRR6148374.cfg;

```

HTML report is available in public directory and here  http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Bombus_t_GHFS01.html
