Data source:
============

### Contigs FASTA: 

```diff
wget https://ftp.ensembl.org/pub/release-110/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
gzip -d Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
```

### GTF (General Transfer Format):


```diff
wget -nc https://ftp.ensembl.org/pub/release-110/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.110.gtf.gz
gzip -d Caenorhabditis_elegans.WBcel235.110.gtf.gz

```

intronSeeker command lines
============================

### Step 1 : Generate contigs fasta

```diff
cd data/GBS-GTF_Based_Simulation/
intronSeeker GTFbasedSimulation -a Caenorhabditis_elegans.WBcel235.110.gtf -r Caenorhabditis_elegans.WBcel235.dna.toplevel.fa -p Cele -o Cele
```

### Step2: Generate reads with intronSeeker simulateReads from reference fasta

```diff
intronSeeker simulateReads -f  Caenorhabditis_elegans.WBcel235.dna.toplevel.fa -c ../../../config/grinder_GBS.cfg -p Cele -o Cele
```

### Step 3: hisat2Alignment

```diff
intronSeeker hisat2Alignment -r Cele/gbs_Cele_transcripts-modified.fa -1 Cele/sr_Cele_R1.fastq.gz -2 Cele/sr_Cele_R2.fastq.gz -o Cele -p Cele

```

### Step 4: splitReadSearch

```diff
intronSeeker splitReadSearch -a Cele/hisat2_Cele.sort.bam -r Cele/gbs_Cele_transcripts-modified.fa -o Cele -p Cele

```

### Step 5: trimFastaFromTXT

```diff
intronSeeker trimFastaFromTXT -r Cele/gbs_Cele_transcripts-modified.fa -c Cele/srs_Cele_HISAT2_candidates.txt -o Cele/HISAT2_trim/ -p Cele
```

### Step 6: Simulation report


Configuration file:

```diff
nano  Cele.cfg
```


```diff
[Defaults]
mfasta:Cele/gbs_Cele_transcripts-modified.fa
fasta:Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
gtf:Caenorhabditis_elegans.WBcel235.110.gtf
r1:Cele/sr_Cele_R1.fastq.gz
r2:Cele/sr_Cele_R2.fastq.gz
flagstat:hisat2_Cele.sort.flagstat.txt
candidat:Cele/srs_Cele_candidates.txt
split:Cele/srs_Cele_split_alignments.txt
rank:Cele/sr_Cele-ranks.txt
prefix:Cele
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  Cele.cfg;
```

HTML report is available in this directory.
