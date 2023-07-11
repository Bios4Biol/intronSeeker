Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/IW/GJIW01/GJIW01.1.fsa_nt.gz
gzip -d GJIW01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR15372351
https://www.ebi.ac.uk/ena/browser/view/SRR15372351

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJIW01.1.fsa_nt -1 SRR15372351_1.fastq -2 SRR15372351_2.fastq --prefix GJIW01 -o GJIW01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJIW01/hisat2_GJIW01.sort.bam -r GJIW01.1.fsa_nt --prefix GJIW01 --output splitReadSearch_GJIW01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GJIW01.1.fsa_nt
r1:SRR15372351_1.fastq
r2:SRR15372351_2.fastq
flagstat:hisat2_SRR15372351.sort.flagstat.txt
candidat:srs_SRR15372351_candidates.txt
split:srs_SRR15372351_split_alignments.txt
prefix:GJIW01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR15372351.cfg;

```

HTML report is available in this directory.
