Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GH/FS/GHFS01/GHFS01.1.fsa_nt.gz
gzip -d GHFS01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR6148374

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GHFS01.1.fsa_nt -1 SRR6148374.fastq.gz --prefix GHFS01 -o GHFS01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GHFS01/hisat2_GHFS01.sort.bam -r GHFS01.1.fsa_nt --prefix GHFS01 --output splitReadSearch_GJJD01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GHFS01.1.fsa_nt
r1:SRR6148374_1.fastq
r2:SRR6148374_2.fastq
flagstat:hisat2_SRR6148374.sort.flagstat.txt
candidat:srs_SRR6148374_candidates.txt
split:srs_SRR6148374_split_alignments.txt
prefix:GHFS01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /path/to/simulation2HTML.py -F --config_file  SRR6148374.cfg;

```

HTML report is available in this directory.
