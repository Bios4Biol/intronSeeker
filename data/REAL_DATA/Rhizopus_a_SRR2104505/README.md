Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GD/UK/GDUK01/GDUK01.1.fsa_nt.gz
gzip -d GDUK01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR2104505
https://www.ebi.ac.uk/ena/browser/view/SRR2104505

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GDUK01.1.fsa_nt -1 SRR2104505_1.fastq -2 SRR2104505_2.fastq --prefix GDUK01 -o GDUK01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GDUK01/hisat2_GDUK01.sort.bam -r GDUK01.1.fsa_nt --prefix GDUK01 --output splitReadSearch_GDUK01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GDUK01.1.fsa_nt
r1:SRR2104505_1.fastq
r2:SRR2104505_2.fastq
flagstat:hisat2_SRR2104505.sort.flagstat.txt
candidat:srs_SRR2104505_candidates.txt
split:srs_SRR2104505_split_alignments.txt
prefix:GDUK01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR2104505.cfg;

```

HTML report is available in this directory.
