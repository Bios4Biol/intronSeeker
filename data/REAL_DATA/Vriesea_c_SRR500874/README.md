Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GH/CB/GHCB01/GHCB01.1.fsa_nt.gz
gzip -d GHCB01.1.fsa_nt.gz
```

### Single reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR500/SRR500874/SRR500874.fastq.gz

```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
intronSeeker hisat2Alignment -r GHCB01.1.fsa_nt -1 SRR500874.fastq.gz --prefix GHCB01 -o GHCB01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GHCB01/hisat2_GHCB01.sort.bam -r GHCB01.1.fsa_nt --prefix GHCB01 --output splitReadSearch_GHCB01
```

### Step 3: run simulation2HTML

Configuration file:
```diff
nano SRR500874.cfg
```

```diff
[Defaults]
fasta:GHCB01.1.fsa_nt
r1:SRR500874.fastq.gz
flagstat:GHCB01/hisat2_GHCB01.sort.flagstat.txt
candidat:splitReadSearch_GHCB01/srs_GHCB01_candidates.txt
split:splitReadSearch_GHCB01/srs_GHCB01_split_alignments.txt
prefix:GHCB01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport -F --config_file  SRR500874.cfg;

```

HTML report is available in public directory and here https://bios4biol.pages.mia.inra.fr/intronseeker/report_Vriesea_c_GHCB01.html

