Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/GJ/DA/GJDA01/GJDA01.1.fsa_nt.gz
gzip -d GJDA01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/072/SRR14143372/SRR14143372_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/072/SRR14143372/SRR14143372_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJDA01.1.fsa_nt -1 SRR14143372_1.fastq -2 SRR14143372_2.fastq --prefix GJDA01 -o GJDA01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJDA01/hisat2_GJDA01.sort.bam -r GJDA01.1.fsa_nt --prefix GJDA01 --output splitReadSearch_GJDA01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GJDA01.1.fsa_nt
r1:SRR14143372_1.fastq
r2:SRR14143372_2.fastq
flagstat:hisat2_SRR14143372.sort.flagstat.txt
candidat:srs_SRR14143372_candidates.txt
split:srs_SRR14143372_split_alignments.txt
prefix:GJDA01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR14143372.cfg;

```

HTML report is available in this directory.
