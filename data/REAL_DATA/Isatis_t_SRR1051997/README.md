Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GA/RR/GARR01/GARR01.1.fsa_nt.gz
gzip -d GARR01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/007/SRR1051997/SRR1051997_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/007/SRR1051997/SRR1051997_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GARR01.1.fsa_nt -1 SRR1051997_1.fastq.gz -2 SRR1051997_2.fastq.gz --prefix GARR01 -o GARR01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GARR01/hisat2_GARR01.sort.bam -r GARR01.1.fsa_nt --prefix GARR01 --output splitReadSearch_GARR01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano  SRR1051997.cfg
```

```diff
[Defaults]
fasta:GARR01.1.fsa_nt
r1:SRR1051997_1.fastq.gz
r2:SRR1051997_2.fastq.gz
flagstat:GARR01/hisat2_GARR01.sort.flagstat.txt
candidat:splitReadSearch_GARR01/srs_GARR01_candidates.txt
split:splitReadSearch_GARR01/srs_GARR01_split_alignments.txt
prefix:GARR01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR1051997.cfg;

```

HTML report is available in public directory and here https://bios4biol.pages.mia.inra.fr/intronseeker/report_Isatis_t_GARR01.html
