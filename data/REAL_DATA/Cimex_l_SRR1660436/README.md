Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/YH/GBYH01/GBYH01.1.fsa_nt.gz
gzip -d GBYH01.1.fsa_nt.gz
```

### Single reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/006/SRR1660436/SRR1660436.fastq.gz
```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GBYH01.1.fsa_nt -1 SRR1660436.fastq.gz  --prefix GBYH01 -o GBYH01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GBYH01/hisat2_GBYH01.sort.bam -r GBYH01.1.fsa_nt --prefix GBYH01 --output splitReadSearch_GBYH01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano  SRR1660436.cfg
```


```diff
[Defaults]
fasta:GBYH01.1.fsa_nt
r:SRR1660436.fastq.gz
flagstat:GBYH01/hisat2_GBYH01.sort.flagstat.txt
candidat:splitReadSearch_GBYH01/srs_GBYH01_candidates.txt
split:splitReadSearch_GBYH01/srs_GBYH01_split_alignments.txt
prefix:GBYH01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR1660436.cfg;

```

HTML report is available in this directory.
