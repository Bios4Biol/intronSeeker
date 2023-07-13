Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/GG/XH/GGXH01/GGXH01.1.fsa_nt.gz
gzip -d GGXH01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR781/005/SRR7819335/SRR7819335_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR781/005/SRR7819335/SRR7819335_2.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GGXH01.1.fsa_nt -1 SRR7819335_1.fastq -2 SRR7819335_2.fastq --prefix GGXH01 -o GGXH01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GGXH01/hisat2_GGXH01.sort.bam -r GGXH01.1.fsa_nt --prefix GGXH01 --output splitReadSearch_GGXH01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:GGXH01.1.fsa_nt
r1:SRR7819335_1.fastq
r2:SRR7819335_2.fastq
flagstat:hisat2_SRR7819335.sort.flagstat.txt
candidat:srs_SRR7819335_candidates.txt
split:srs_SRR7819335_split_alignments.txt
prefix:GGXH01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR7819335.cfg;

```

HTML report is available in this directory.
