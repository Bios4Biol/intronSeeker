Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/GX/GBGX01/GBGX01.1.fsa_nt.gz
gzip -d GBGX01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR1206033
https://www.ebi.ac.uk/ena/browser/view/SRR1206033

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GBGX01.1.fsa_nt -1 SRR1206033.fastq  --prefix GBGX01 -o GBGX01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GBGX01/hisat2_GBGX01.sort.bam -r GBGX01.1.fsa_nt --prefix GBGX01 --output splitReadSearch_GBGX01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GBGX01.1.fsa_nt
r:SRR1206033_1.fastq
flagstat:hisat2_SRR1206033.sort.flagstat.txt
candidat:srs_SRR1206033_candidates.txt
split:srs_SRR1206033_split_alignments.txt
prefix:GBGX01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR1206033.cfg;

```

HTML report is available in this directory.
