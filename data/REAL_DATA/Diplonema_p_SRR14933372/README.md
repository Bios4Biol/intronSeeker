Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/NJ/GJNJ01/GJNJ01.1.fsa_nt.gz
gzip -d GJNJ01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR14933372
https://www.ebi.ac.uk/ena/browser/view/SRR14933372

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJNJ01.1.fsa_nt -1 SRR14933372_1.fastq -2 SRR14933372_2.fastq --prefix GJNJ01 -o GJNJ01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJNJ01/hisat2_GJNJ01.sort.bam -r GJNJ01.1.fsa_nt --prefix GJNJ01 --output splitReadSearch_GJNJ01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GJNJ01.1.fsa_nt
r1:SRR14933372_1.fastq
r2:SRR14933372_2.fastq
flagstat:hisat2_SRR14933372.sort.flagstat.txt
candidat:srs_SRR14933372_candidates.txt
split:srs_SRR14933372_split_alignments.txt
prefix:GJNJ01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR14933372.cfg;

```

HTML report is available in this directory.