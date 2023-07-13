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
intronSeeker hisat2Alignment -r GAQX01.1.fsa_nt -1 SRR857257_1.fastq -2 SRR857257_2.fastq --prefix GAQX01 -o GAQX01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GAQX01/hisat2_GAQX01.sort.bam -r GAQX01.1.fsa_nt --prefix GAQX01 --output splitReadSearch_GAQX01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GAQX01.1.fsa_nt
r1:SRR857257_1.fastq
r2:SRR857257_2.fastq
flagstat:hisat2_SRR857257.sort.flagstat.txt
candidat:srs_SRR857257_candidates.txt
split:srs_SRR857257_split_alignments.txt
prefix:GAQX01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR857257.cfg;

```

HTML report is available in this directory.
