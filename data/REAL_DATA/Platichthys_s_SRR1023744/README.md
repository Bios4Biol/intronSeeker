Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GA/PK/GAPK01/GAPK01.1.fsa_nt.gz
gzip -d GAPK01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR1023744
https://www.ebi.ac.uk/ena/browser/view/SRR1023744

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GAPK01.1.fsa_nt -1 SRR1023744_1.fastq -2 SRR1023744_2.fastq --prefix GAPK01 -o GAPK01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GAPK01/hisat2_GAPK01.sort.bam -r GAPK01.1.fsa_nt --prefix GAPK01 --output splitReadSearch_GAPK01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GAPK01.1.fsa_nt
r1:SRR1023744_1.fastq
r2:SRR1023744_2.fastq
flagstat:hisat2_SRR1023744.sort.flagstat.txt
candidat:srs_SRR1023744_candidates.txt
split:srs_SRR1023744_split_alignments.txt
prefix:GAPK01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR1023744.cfg;

```

HTML report is available in this directory.
