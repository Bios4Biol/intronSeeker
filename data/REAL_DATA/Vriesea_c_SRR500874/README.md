Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GH/CB/GHCB01/GHCB01.1.fsa_nt.gz
gzip -d GHCB01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR500874
https://www.ebi.ac.uk/ena/browser/view/SRR500874

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GHCB01.1.fsa_nt -1 SRR500874_1.fastq --prefix GHCB01 -o GHCB01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GHCB01/hisat2_GHCB01.sort.bam -r GHCB01.1.fsa_nt --prefix GHCB01 --output splitReadSearch_GHCB01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GHCB01.1.fsa_nt
r:SRR500874_1.fastq
flagstat:hisat2_SRR500874.sort.flagstat.txt
candidat:srs_SRR500874_candidates.txt
split:srs_SRR500874_split_alignments.txt
prefix:GHCB01
threads: 6                
output:HTML/
force: -F
```


```diff
python3 /PATH/TO/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR500874.cfg;

```

HTML report is available in this directory.
