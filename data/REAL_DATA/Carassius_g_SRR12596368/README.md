Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/KR/GJKR01/GJKR01.1.fsa_nt.gz
gzip -d GJKR01.1.fsa_nt.gz
```

### Paired reads:


```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/068/SRR12596368/SRR12596368_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/068/SRR12596368/SRR12596368_1.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJKR01.1.fsa_nt -1 SRR12596368_1.fastq.gz -2 SRR12596368_2.fastq.gz --prefix GJKR01 -o GJKR01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJKR01/hisat2_GJKR01.sort.bam -r GJKR01.1.fsa_nt --prefix GJKR01 --output splitReadSearch_GJKR01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GJKR01.1.fsa_nt
r1:SRR12596368_1.fastq
r2:SRR12596368_2.fastq
flagstat:hisat2_SRR12596368.sort.flagstat.txt
candidat:srs_SRR12596368_candidates.txt
split:srs_SRR12596368_split_alignments.txt
prefix:GJKR01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR12596368.cfg;

```

HTML report is available in this directory.
