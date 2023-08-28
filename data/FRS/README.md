FRS -Full Random Simulation
---------------------------

##### 1 - Step 1 : Full Random Simulation

```diff
source activate ISeeker_environment
(ISeeker_environment) bash-4.4$ intronSeeker fullRandomSimulation -n 1000 -r  -o sampleFRS/ -p sampleFRS
```

Results:

```diff
(ISeeker_environment) $ ls sampleFRS/
frs_sampleFRS_contigs.fa  frs_sampleFRS_contigs-modified.fa  frs_sampleFRS_contigs-modified.gtf
```

##### 2 - Step 2 : Simulate Reads

```diff
(ISeeker_environment) bash-4.4$ intronSeeker simulateReads -o ./ -f sampleFRS/frs_sampleFRS_contigs.fa -c grinder_frs.cfg
###  Start to simulate reads (grinder)   ###


```

##### 3 - Step 3 : STAR (or HISAT2) alignment and Split Read Search

```diff
(ISeeker_environment) $  intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ ls 
frs_sample1_contigs.fa  frs_sample1_contigs-modified.fa  frs_sample1_contigs-modified.gtf  sr.log  sr_R1.fastq.gz  sr_R2.fastq.gz  sr_ranks.txt  STAR_alignment
```



##### IntronSeeker simulation Report 

For each sample:

```diff
intronSeeker buildReport -F ...
```
