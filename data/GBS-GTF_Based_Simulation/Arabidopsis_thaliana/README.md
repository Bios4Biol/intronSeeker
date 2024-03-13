Data source:
============

### Contigs FASTA: 

```diff
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gzip -d Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```

### GTF (General Transfer Format):


```diff
wget -nc https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gtf.gz
gzip -d Arabidopsis_thaliana.TAIR10.57.gtf.gz

```

intronSeeker command lines
============================

### Step 1 : Generate contigs fasta

```diff
cd data/GBS-GTF_Based_Simulation/
intronSeeker GTFbasedSimulation -a Arabidopsis_thaliana.TAIR10.57.gtf -r Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -p Athal -o Athal
```

### Step2: Generate reads with intronSeeker simulateReads from reference fasta

```diff
intronSeeker simulateReads -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -c ../../../config/grinder_GBS.cfg -p Athal -o Athal
```

### Step 3: hisat2Alignment

```diff
intronSeeker hisat2Alignment -r Athal/gbs_Athal_transcripts-modified.fa -1 Athal/sr_Athal_R1.fastq.gz -2 Athal/sr_Athal_R2.fastq.gz -o Athal -p Athal

```

### Step 4: splitReadSearch

```diff
intronSeeker splitReadSearch -a Athal/hisat2_Athal.sort.bam -r Athal/gbs_Athal_transcripts-modified.fa -o Athal -p Athal

```

### Step 5: trimFastaFromTXT

```diff
intronSeeker trimFastaFromTXT -r Athal/gbs_Athal_transcripts-modified.fa -c Athal/srs_Athal_HISAT2_candidates.txt -o Athal/HISAT2_trim/ -p Athal
```

### Step 6: Simulation report


Configuration file:

```diff
nano  Athal.cfg
```


```diff
[Defaults]
mfasta:Athal/gbs_Athal_transcripts-modified.fa
fasta:Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
gtf:Athal/gbs_Athal_transcripts-modified.gtf
r1:Athal/sr_Athal_R1.fastq.gz
r2:Athal/sr_Athal_R2.fastq.gz
flagstat:Athal/hisat2_Athal.sort.flagstat.txt
candidat:Athal/srs_Athal_candidates.txt
split:Athal/srs_Athal_split_alignments.txt
ranks:Athal/sr_Athal_ranks.txt
prefix:Athal
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  Athal.cfg;
```

HTML report is available in intronSeeker public directory : http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Athal.html
