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
intronSeeker.py GTFbasedSimulation -a Caenorhabditis_elegans.WBcel235.110.gtf -r Caenorhabditis_elegans.WBcel235.dna.toplevel.fa -p "Athal" -o Athal
```

### Step2: Generate reads with intronSeeker simulateReads from reference fasta

```diff
intronSeeker.py simulateReads -f "+ref+" -c "+grinder+" -p "+pipelineName+" -o "+pipelineName
```

### Step 3: hisat2Alignment

```diff
intronSeeker hisat2Alignment -r "+mref+" -1 Athal/sr_Athal_R1.fastq.gz -2 Athal/sr_Athal_R2.fastq.gz -o Athal -p Athal

```

### Step 4: splitReadSearch

```diff
intronSeeker splitReadSearch -a Athal/hisat2_Athal.sort.bam -r "+mref+" -o Athal -p Athal

```

### Step 5: trimFastaFromTXT

```diff
intronSeeker.py trimFastaFromTXT -r "+mref+" -c Athal/srs_Athal_HISAT2_candidates.txt -o Athal/HISAT2_trim/ -p Athal


```

### Step 6: Simulation report


Configuration file:

```diff
nano  Athal.cfg
```


```diff
[Defaults]
mfasta:gbs_Athal_transcripts-modified.fa
fasta:
gtf:Caenorhabditis_elegans.WBcel235.110.gtf
r1:Athal.fastq
r2:Athal.fastq
flagstat:hisat2_Athal.sort.flagstat.txt
candidat:srs_Athal_candidates.txt
split:srs_Athal_split_alignments.txt
prefix:Athal
threads: 6                
output:HTML/
force: -F
```


```diff
intronSeeker buildReport  -F --config_file  Athal.cfg;

```

HTML report is available in this directory.