Data source:
============

### Contigs FASTA: 

```diff
wget https://ftp.ensembl.org/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
gzip -d Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
```

### GTF (General Transfer Format):


```diff
wget -nc https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz
gzip -d Drosophila_melanogaster.BDGP6.46.110.gtf.gz

```

intronSeeker command lines
============================

### Step 1 : Generate contigs fasta

```diff
cd data/GBS-GTF_Based_Simulation/
intronSeeker GTFbasedSimulation -a Drosophila_melanogaster.BDGP6.46.110.gtf -r Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa -p Dmel -o Dmel
```

### Step2: Generate reads with intronSeeker simulateReads from reference fasta

```diff
intronSeeker simulateReads -f Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa -c ../../../config/grinder_GBS.cfg -p Dmel -o Dmel
```

### Step 3: hisat2Alignment

```diff
intronSeeker hisat2Alignment -r Dmel/gbs_Dmel_transcripts-modified.fa -1 Dmel/sr_Dmel_R1.fastq.gz -2 Dmel/sr_Dmel_R2.fastq.gz -o Dmel -p Dmel

```

### Step 4: splitReadSearch

```diff
intronSeeker splitReadSearch -a Dmel/hisat2_Dmel.sort.bam -r "+mref+" -o Dmel -p Dmel

```

### Step 5: trimFastaFromTXT

```diff
intronSeeker trimFastaFromTXT -r Dmel/gbs_Dmel_transcripts-modified.fa -c Dmel/srs_Dmel_HISAT2_candidates.txt -o Dmel/HISAT2_trim/ -p Dmel
```

### Step 6: Simulation report


Configuration file:

```diff
nano  Dmel.cfg
```


```diff
[Defaults]
mfasta:Dmel/gbs_Dmel_transcripts-modified.fa
fasta:Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa
gtf:Drosophila_melanogaster.BDGP6.46.110.gtf
r1:Dmel/sr_Dmel_R1.fastq.gz
r2:Dmel/sr_Dmel_R2.fastq.gz
flagstat:hisat2_Dmel.sort.flagstat.txt
candidat:Dmel/srs_Dmel_HISAT2_candidates.txt
split:Dmel/srs_Dmel_split_alignments.txt
rank:Dmel/sr_Dmel-ranks.txt
prefix:Dmel
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  Dmel.cfg;
```

HTML report is available in this directory.
