Data source:
============

### Contigs FASTA: 

```diff
wget https://ftp.ensembl.org/pub/release-110/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
gzip -d Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
```

### GTF (General Transfer Format):


```diff
wget -nc https://ftp.ensembl.org/pub/release-110/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf.gz
gzip -d Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf

```

intronSeeker command lines
============================

### Step 1 : Generate contigs fasta

```diff
cd data/GBS-GTF_Based_Simulation/
intronSeeker GTFbasedSimulation -a Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf -r Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -p Ggal -o Ggal
```

### Step2: Generate reads with intronSeeker simulateReads from reference fasta

```diff
intronSeeker simulateReads -f Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -c ../../../config/grinder_GBS.cfg -p Ggal -o Ggal
```

### Step 3: hisat2Alignment

```diff
intronSeeker hisat2Alignment -r Ggal/gbs_Ggal_transcripts-modified.fa -1 Ggal/sr_Ggal_R1.fastq.gz -2 Ggal/sr_Ggal_R2.fastq.gz -o Ggal -p Ggal

```

### Step 4: splitReadSearch

```diff
intronSeeker splitReadSearch -a Ggal/hisat2_Ggal.sort.bam -r Ggal/gbs_Ggal_transcripts-modified.fa -o Ggal -p Ggal

```

### Step 5: trimFastaFromTXT

```diff
intronSeeker.py trimFastaFromTXT -r Ggal/gbs_Ggal_transcripts-modified.fa -c Ggal/srs_Ggal_HISAT2_candidates.txt -o Ggal/HISAT2_trim/ -p Ggal
```

### Step 6: Simulation report


Configuration file:

```diff
nano  Ggal.cfg
```


```diff
[Defaults]
mfasta:Ggal/gbs_Ggal_transcripts-modified.fa
fasta:Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa
gtf:Ggal/gbs_Ggal_transcripts-modified.gtf
r1:Ggal/sr_Ggal_R1.fastq.gz
r2:Ggal/sr_Ggal_R2.fastq.gz
flagstat:Ggal/hisat2_Ggal.sort.flagstat.txt
candidat:Ggal/srs_Ggal_candidates.txt
split:Ggal/srs_Ggal_split_alignments.txt
ranks:Ggal/sr_Ggal_ranks.txt
prefix:Ggal
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  Ggal.cfg;
```

HTML report is available in intronSeeker public directory : http://htmlpreview.github.io/?https://github.com/Bios4Biol/intronSeeker/blob/master/public/report_Ggal.html
