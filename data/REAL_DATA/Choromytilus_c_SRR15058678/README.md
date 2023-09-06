Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/JD/GJJD01/GJJD01.1.fsa_nt.gz
gzip -d GJJD01.1.fsa_nt.gz
```

### Single reads:


```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/078/SRR15058678/SRR15058678.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJJD01.1.fsa_nt -1 SRR15058678_1.fastq.gz  --prefix GJJD01 -o GJJD01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJJD01/hisat2_GJJD01.sort.bam -r GJJD01.1.fsa_nt --prefix GJJD01 --output splitReadSearch_GJJD01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano  SRR15058678.cfg
```


```diff
[Defaults]
fasta:GJJD01.1.fsa_nt
r1:SRR15058678_1.fastq.gz
flagstat:GJJD01/hisat2_GJJD01.sort.flagstat.txt
candidat:splitReadSearch_GJJD01/srs_GJJD01_candidates.txt
split:splitReadSearch_GJJD01/srs_GJJD01_split_alignments.txt
prefix:GJJD01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR15058678.cfg;

```

HTML report is available in public directory and here https://emilien.lasguignes.pages.mia.inra.fr/intronSeeker/report_Choromytilus_c_GJJD01.html
