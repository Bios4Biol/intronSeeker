Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GA/RR/GARR01/GARR01.1.fsa_nt.gz
gzip -d GARR01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR1051997
https://www.ebi.ac.uk/ena/browser/view/SRR1051997

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GARR01.1.fsa_nt -1 SRR1051997_1.fastq -2 SRR1051997_2.fastq --prefix GARR01 -o GARR01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GARR01/hisat2_GARR01.sort.bam -r GARR01.1.fsa_nt --prefix GARR01 --output splitReadSearch_GARR01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GARR01.1.fsa_nt
r1:SRR1051997_1.fastq
r2:SRR1051997_2.fastq
flagstat:hisat2_SRR1051997.sort.flagstat.txt
candidat:srs_SRR1051997_candidates.txt
split:srs_SRR1051997_split_alignments.txt
prefix:GARR01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR1051997.cfg;

```

HTML report is available in this directory.
