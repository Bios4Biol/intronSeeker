Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/JN/GJJN01/GJJN01.1.fsa_nt.gz
gzip -d GJJN01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR1618559
```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJJN01.1.fsa_nt -1 SRR1618559.fastq.gz --prefix GJJN01 -o GJJN01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJJN01/hisat2_GJJN01.sort.bam -r GJJN01.1.fsa_nt --prefix GJJN01 --output splitReadSearch_GJJN01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GJJN01.1.fsa_nt
r1:SRR1618559_1.fastq
r2:SRR1618559_2.fastq
flagstat:hisat2_SRR1618559.sort.flagstat.txt
candidat:srs_SRR1618559_candidates.txt
split:srs_SRR1618559_split_alignments.txt
prefix:GJJN01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/simulation2HTML.py -F --config_file  SRR1618559.cfg;

```

HTML report is available in this directory.
