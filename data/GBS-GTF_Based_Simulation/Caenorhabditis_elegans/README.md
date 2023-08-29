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
intronSeeker.py GTFbasedSimulation -a Caenorhabditis_elegans.WBcel235.110.gtf -r Caenorhabditis_elegans.WBcel235.dna.toplevel.fa -p "Cele" -o Cele
```

### Step2: Generate reads with intronSeeker simulateReads from reference fasta

```diff
intronSeeker.py simulateReads -f "+ref+" -c "+grinder+" -p "+pipelineName+" -o "+pipelineName
```

### Step 3: hisat2Alignment

```diff
intronSeeker hisat2Alignment -r "+mref+" -1 Cele/sr_Cele_R1.fastq.gz -2 Cele/sr_Cele_R2.fastq.gz -o Cele -p Cele

```

### Step 4: splitReadSearch

```diff
intronSeeker splitReadSearch -a Cele/hisat2_Cele.sort.bam -r "+mref+" -o Cele -p Cele

```

### Step 5: trimFastaFromTXT

```diff
intronSeeker.py trimFastaFromTXT -r "+mref+" -c Cele/srs_Cele_HISAT2_candidates.txt -o Cele/HISAT2_trim/ -p Cele


```

### Step 6: Simulation report


Configuration file:

```diff
nano  Cele.cfg
```


```diff
[Defaults]
mfasta:gbs_Cele_transcripts-modified.fa
fasta:
gtf:Caenorhabditis_elegans.WBcel235.110.gtf
r1:Cele.fastq
r2:Cele.fastq
flagstat:hisat2_Cele.sort.flagstat.txt
candidat:srs_Cele_candidates.txt
split:srs_Cele_split_alignments.txt
prefix:Cele
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  Cele.cfg;

```

HTML report is available in this directory.
