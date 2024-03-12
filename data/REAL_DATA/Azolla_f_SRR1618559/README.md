Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/GB/TV/GBTV01/GBTV01.1.fsa_nt.gz
gzip -d GBTV01.1.fsa_nt.gz
```

### Paired reads:

```diff
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/009/SRR1618559/SRR1618559_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/009/SRR1618559/SRR1618559_1.fastq.gz
```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
 intronSeeker hisat2Alignment -r GBTV01.1.fsa_nt  -1 SRR1618559_1.fastq.gz -2 SRR1618559_2.fastq.gz --prefix GBTV01  -o GBTV01 -t 12
```

### Step2: splitReadSearch

```diff
intronSeeker splitReadSearch -a GBTV01/hisat2_GBTV01.sort.bam -r GBTV01.1.fsa_nt --prefix GBTV01 --output splitReadSearch_GBTV01
```

### Step 3: run simulation2HTML

Configuration file SRR1618559.cfg:

```diff
nano  SRR1618559.cfg
```

```diff
[Defaults]
fasta:GBTV01.1.fsa_nt
r1:SRR1618559_1.fastq.gz
r2:SRR1618559_2.fastq.gz
flagstat:GBTV01/hisat2_GBTV01.sort.flagstat.txt
candidat:splitReadSearch_GBTV01/srs_GBTV01_candidates.txt
split:splitReadSearch_GBTV01/srs_GBTV01_split_alignments.txt
prefix:GBTV01
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  SRR1618559.cfg;

```

HTML report is available in public directory and here : http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Azolla_f_GBTV01.html
