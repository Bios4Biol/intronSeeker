Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/YH/GBYH01/GBYH01.1.fsa_nt.gz
gzip -d GBYH01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR1660436
https://www.ebi.ac.uk/ena/browser/view/SRR1660436

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GBYH01.1.fsa_nt -1 SRR1660436.fastq  --prefix GBYH01 -o GBYH01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GBYH01/hisat2_GBYH01.sort.bam -r GBYH01.1.fsa_nt --prefix GBYH01 --output splitReadSearch_GBYH01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GBYH01.1.fsa_nt
r:SRR1660436.fastq
flagstat:hisat2_SRR1660436.sort.flagstat.txt
candidat:srs_SRR1660436_candidates.txt
split:srs_SRR1660436_split_alignments.txt
prefix:GBYH01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR1660436.cfg;

```

HTML report is available in this directory.