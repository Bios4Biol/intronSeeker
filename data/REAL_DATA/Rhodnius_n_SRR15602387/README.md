Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/JI/GJJI01/GJJI01.1.fsa_nt.gz
gzip -d GJJI01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR15602387
https://www.ebi.ac.uk/ena/browser/view/SRR15602387

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJJI01.1.fsa_nt -1 SRR15602387_1.fastq -2 SRR15602387_2.fastq --prefix GJJI01 -o GJJI01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJJI01/hisat2_GJJI01.sort.bam -r GJJI01.1.fsa_nt --prefix GJJI01 --output splitReadSearch_GJJI01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GJJI01.1.fsa_nt
r1:SRR15602387_1.fastq
r2:SRR15602387_2.fastq
flagstat:hisat2_SRR15602387.sort.flagstat.txt
candidat:srs_SRR15602387_candidates.txt
split:srs_SRR15602387_split_alignments.txt
prefix:GJJI01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR15602387.cfg;

```

HTML report is available in this directory.
