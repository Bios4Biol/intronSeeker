Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/GG/XH/GGXH01/GGXH01.1.fsa_nt.gz
gzip -d GGXH01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR781/005/SRR7819335/SRR7819335_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR781/005/SRR7819335/SRR7819335_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GGXH01.1.fsa_nt -1 SRR7819335_1.fastq.gz -2 SRR7819335_2.fastq.gz --prefix GGXH01 -o GGXH01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GGXH01/hisat2_GGXH01.sort.bam -r GGXH01.1.fsa_nt --prefix GGXH01 --output splitReadSearch_GGXH01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano SRR7819335.cfg
```

```diff
[Defaults]
fasta:GGXH01.1.fsa_nt
r1:SRR7819335_1.fastq.gz
r2:SRR7819335_2.fastq.gz
flagstat:GGXH01/hisat2_GGXH01.sort.flagstat.txt
candidat:splitReadSearch_GGXH01/srs_GGXH01_candidates.txt
split:splitReadSearch_GGXH01/srs_GGXH01_split_alignments.txt
prefix:GGXH01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR7819335.cfg;

```

HTML report is available in public directory and here http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Piromyces_sp_GGXH01.html

