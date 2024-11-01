Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GG/UN/GGUN01/GGUN01.1.fsa_nt.gz
gzip -d GGUN01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR760/006/SRR7601946/SRR7601946_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR760/006/SRR7601946/SRR7601946_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GGUN01.1.fsa_nt -1 SRR7601946_1.fastq.gz -2 SRR7601946_2.fastq.gz --prefix GGUN01 -o GGUN01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GGUN01/hisat2_GGUN01.sort.bam -r GGUN01.1.fsa_nt --prefix GGUN01 --output splitReadSearch_GGUN01
```

### Step 3: run simulation2HTML

Configuration file:
```diff
nano SRR7601946.cfg
```

```diff
[Defaults]
fasta:GGUN01.1.fsa_nt
r1:SRR7601946_1.fastq.gz
r2:SRR7601946_2.fastq.gz
flagstat:GGUN01/hisat2_GGUN01.sort.flagstat.txt
candidat:splitReadSearch_GGUN01/srs_GGUN01_candidates.txt
split:splitReadSearch_GGUN01/srs_GGUN01_split_alignments.txt
prefix:GGUN01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR7601946.cfg;

```

HTML report is available in public directory and here http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Goniomonas_a_GGUN01.html
