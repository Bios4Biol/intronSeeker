FRS - Full Random Simulation
----------------------------

##### 1 - Step 1 : Full Random Simulation

```diff
source activate ISeeker_environment
(ISeeker_environment) bash-4.4$ intronSeeker fullRandomSimulation -n 1000 -o sampleFRS/ -p sampleFRS
```

Results:

```diff
(ISeeker_environment) bash-4.4$ ls sampleFRS/
frs_sampleFRS_contigs.fa  frs_sampleFRS_contigs-modified.fa  frs_sampleFRS_contigs-modified.gtf
```

##### 2 - Step 2 : Simulate Reads

```diff
(ISeeker_environment) bash-4.4$ intronSeeker simulateReads -o ./ -f sampleFRS/frs_sampleFRS_contigs.fa -c grinder_frs.cfg
###  Start to simulate reads (grinder)   ###
###  End of reads simulation (grinder)   ###
(ISeeker_environment) bash-4.4$ ls
grinder_frs.cfg
sampleFRS
sr_ranks.txt
sr.log
sr_R2.fastq.gz
sr_R1.fastq.gz
```

##### 3 - Step 3 : STAR or HISAT2 alignment

STAR alignment:

```diff
(ISeeker_environment) bash-4.4$ intronSeeker starAlignment -r sampleFRS/frs_sampleFRS_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
###  Start to map with STAR   ###
##  Fasta indexing   ##
##  STAR alignement   ##
[E::idx_find_and_load] Could not retrieve index file for 'STAR_alignment/star.star.Aligned.sortedByCoord.out.bam'
##  BAM indexing   ##
##  BAM flagstat computation   ##
###  End of STAR mapping   ###
```
Results files for STAR alignment:
```diff
(ISeeker_environment) bash-4.4$ ls -ltrah STAR_alignment/
star_genomeRef
star.star.Log.out
star.sort.bam
star.sort.bai
star.sort.flagstat.txt
star.log
```

HiSAT2 alignment:
```diff
(ISeeker_environment) bash-4.4$ intronSeeker hisat2Alignment -r sampleFRS/frs_sampleFRS_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
###  Start to map with HiSat2   ###
##  Fasta indexing   ##
##  HiSat2 Alignement   ##
##  BAM indexing   ##
##  BAM flagstat computation   ##
###  End of HiSat2 mapping   ###
```

Results files for STAR alignment:
```diff
(ISeeker_environment) bash-4.4$ ls HISAT2_alignment/
hisat2_genomeRef
hisat2_build.log
hisat2_aln.log
hisat2.sort.bam
hisat2.sort.bai
hisat2.sort.flagstat.txt
hisat2.log
```

##### 4 - Step 4 : Split Read Search

```diff
(ISeeker_environment) bash-4.4$ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r sampleFRS/frs_sampleFRS_contigs-modified.fa -o SRS_HISAT2/
###  Start to search split read   ###
Begin candidates and splits search :  2023-08-28 11:12:00.358650
zip list :  2023-08-28 11:17:10.050756
##  Processing ## 
 Preview of candidates list before filter:  (                                       reference  start   end  depth split_borders  DP_before       DP_in  DP_after  filter
SEQUENCE1000.modif|13|303     SEQUENCE1000.modif     13   303    226         AT_CA      326.0  109.924399       3.0  SS;RDP
SEQUENCE1000.modif|13|370     SEQUENCE1000.modif     13   370     35         AT_AT      326.0   89.913408       3.0  SS;RDP
SEQUENCE1000.modif|13|649     SEQUENCE1000.modif     13   649    184         AT_AG      326.0   51.846154       3.0  SS;RDP
SEQUENCE1000.modif|14|681     SEQUENCE1000.modif     14   681    309         TC_AG      327.7   55.708084     187.4  SS;RDP
SEQUENCE1000.modif|14|756     SEQUENCE1000.modif     14   756    248         TC_GA      327.7   60.711978      66.7  SS;RDP
...                                          ...    ...   ...    ...           ...        ...         ...       ...     ...
SEQUENCE1000.modif|1367|2006  SEQUENCE1000.modif   1367  2006     39             _        0.0    0.000000       0.0      SS
SEQUENCE1000.modif|1405|1762  SEQUENCE1000.modif   1405  1762    145             _        0.0    0.000000       0.0      SS
SEQUENCE1000.modif|1419|2103  SEQUENCE1000.modif   1419  2103     86             _        0.0    0.000000       0.0      SS
SEQUENCE1000.modif|1420|1939  SEQUENCE1000.modif   1420  1939    134             _        0.0    0.000000       0.0      SS
SEQUENCE1000.modif|1422|1842  SEQUENCE1000.modif   1422  1842     50             _        0.0    0.000000       0.0      SS

[860 rows x 9 columns],) 

 Preview of splits list before filter:  (        start_split  end_split  split_length split_borders      read           reference strand
0               171        673           503         CT_AC  481634/1     SEQUENCE1.modif      +
1               171        673           503         CT_AC  145967/1     SEQUENCE1.modif      +
2               171        673           503         CT_AC  815467/2     SEQUENCE1.modif      +
3               171        673           503         CT_AC  164936/1     SEQUENCE1.modif      -
4               171        673           503         CT_AC  322273/2     SEQUENCE1.modif      +
...             ...        ...           ...           ...       ...                 ...    ...
184619           91        668           578         CT_AC  807846/1  SEQUENCE1000.modif      +
184620           91        668           578         CT_AC      74/2  SEQUENCE1000.modif      +
184621           91        668           578         CT_AC  208535/2  SEQUENCE1000.modif      +
184622           91        668           578         CT_AC  307374/2  SEQUENCE1000.modif      +
184623           91        668           578         CT_AC  985697/1  SEQUENCE1000.modif      +

[184624 rows x 7 columns],)
End of candidates and split search :  2023-08-28 11:17:10.079779
###  Filter   ###
Begin candidates and splits filter :  2023-08-28 11:17:10.102088
###  Focus on flagged PASS candidates and modify filter field (from PASS to OI) if two "retained introns" are overlapping  ###
End of candidates and split filter :  2023-08-28 11:17:10.133371
###  Write split alignment file   ###
###  Write candidates file   ###
###  End of searching split read   ###
```			     

SRS results files:

```diff
(ISeeker_environment) bash-4.4$ ls SRS_HISAT2/
srs_split_alignments.txt
srs_candidates.txt
```	
			     
##### 5 - Step 5 : IntronSeeker simulation Report 

```diff
more FRS.config
[Defaults]
fasta:sampleFRS/frs_sampleFRS_contigs.fa
mfasta:sampleFRS/frs_sampleFRS_contigs-modified.fa
r1:sr_R1.fastq.gz
r2:sr_R2.fastq.gz
flagstat:STAR_alignment/star.sort.flagstat.txt
candidat:SRS_HISAT2/srs_candidates.txt
split:SRS_HISAT2/srs_split_alignments.txt
ranks:sr_ranks.txt
prefix:FRS
threads:6             
output:HTML/
force: -F

(ISeeker_environment) bash-4.4$ intronSeeker buildReport  -F --config_file FRS.config
```

HTML report result file:
```diff
(ISeeker_environment) bash-4.4$ ls HTML/
report_FRS.html
```
