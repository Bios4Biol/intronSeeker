Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/IW/GJIW01/GJIW01.1.fsa_nt.gz
gzip -d GJIW01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/051/SRR15372351/SRR15372351_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/051/SRR15372351/SRR15372351_1.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GJIW01.1.fsa_nt -1 SRR15372351_1.fastq.gz -2 SRR15372351_2.fastq.gz --prefix GJIW01 -o GJIW01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GJIW01/hisat2_GJIW01.sort.bam -r GJIW01.1.fsa_nt --prefix GJIW01 --output splitReadSearch_GJIW01
```

### Step 3: run simulation2HTML

Configuration file:

```diff
nano SRR15372351.cfg
```

```diff
[Defaults]
fasta:GJIW01.1.fsa_nt
r1:SRR15372351_1.fastq
r2:SRR15372351_2.fastq
flagstat:GJIW01/hisat2_SRR15372351.sort.flagstat.txt
candidat:splitReadSearch_GJIW01/srs_SRR15372351_candidates.txt
split:splitReadSearch_GJIW01/srs_SRR15372351_split_alignments.txt
prefix:GJIW01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR15372351.cfg;

```

HTML report is available in this directory.
