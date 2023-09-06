Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GA/QX/GAQX01/GAQX01.1.fsa_nt.gz
gzip -d GAQX01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/SRR857257/SRR857257_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/SRR857257/SRR857257_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GAQX01.1.fsa_nt -1 SRR857257_1.fastq.gz -2 SRR857257_2.fastq.gz --prefix GAQX01 -o GAQX01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GAQX01/hisat2_GAQX01.sort.bam -r GAQX01.1.fsa_nt --prefix GAQX01 --output splitReadSearch_GAQX01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano SRR857257.cfg
```


```diff
[Defaults]
fasta:GAQX01.1.fsa_nt
r1:SRR857257_1.fastq.gz
r2:SRR857257_2.fastq.gz
flagstat:GAQX01/hisat2_GAQX01.sort.flagstat.txt
candidat:splitReadSearch_GAQX01/srs_GAQX01_candidates.txt
split:splitReadSearch_GAQX01/srs_GAQX01_split_alignments.txt
prefix:GAQX01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR857257.cfg;

```

HTML report is available in public directory and here https://emilien.lasguignes.pages.mia.inra.fr/intronSeeker/report_Graminella_n_GAQX01.html
