Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/NJ/GJNJ01/GJNJ01.1.fsa_nt.gz
gzip -d GJNJ01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR149/072/SRR14933372/SRR14933372_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR149/072/SRR14933372/SRR14933372_1.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJNJ01.1.fsa_nt -1 SRR14933372_1.fastq.gz -2 SRR14933372_2.fastq.gz --prefix GJNJ01 -o GJNJ01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJNJ01/hisat2_GJNJ01.sort.bam -r GJNJ01.1.fsa_nt --prefix GJNJ01 --output splitReadSearch_GJNJ01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano  SRR14933372.cfg
```


```diff
[Defaults]
fasta:GJNJ01.1.fsa_nt
r1:SRR14933372_1.fastq.gz
r2:SRR14933372_2.fastq.gz
flagstat:GJNJ01/hisat2_GJNJ01.sort.flagstat.txt
candidat:splitReadSearch_GJNJ01/srs_GJNJ01_candidates.txt
split:splitReadSearch_GJNJ01/srs_GJNJ01_split_alignments.txt
prefix:GJNJ01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR14933372.cfg;

```

HTML report is available in this directory.
