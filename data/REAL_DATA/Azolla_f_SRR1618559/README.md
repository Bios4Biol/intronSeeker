Data source:
============

### Contigs FASTA: 

```diff
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/GJ/JN/GJJN01/GJJN01.1.fsa_nt.gz
```

### Paired reads:

Use SRA Toolkit (https://github.com/ncbi/sra-tools/wiki) to download runs locally.

```diff
fastq-dump -X 5  --split-files  SRR1618559
```

intronSeeker command lines
============================

### Step 1 : hisat2Alignment

```diff
#Load modules
module load devel/Miniconda/Miniconda3
source activate ISeeker_environment

# working dir 
cd /PATH/TO/YOUR/WORKINGDIR/;

# Command line
intronSeeker hisat2Alignment -r GBTV01.1.fsa_nt -1 SRR1618559_1.fastq -2 SRR1618559_2.fastq  -p  SRR1618559  -o /PATH/TO/YOUR/OUTPUT/PATH/;
```

### Step2: splitReadSearch

```diff

#Load modules
module load devel/Miniconda/Miniconda3
source activate ISeeker_environment

# working dir 
cd /work/user/smaman/23032022_intronSeeker/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/;

# Command line
intronSeeker splitReadSearch -t 6 -a hisat2_SRR1618559.sort.bam -r GBTV01.1.fsa_nt  -o /work/user/smaman/23032022_intronSeeker/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/ -p SRR1618559 -d 1;

```

### Step 3: run simulation2HTML

Configuration file:

```diff
[Defaults]
fasta:/work/smaman/ENCOURS/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/GBTV01.1.fsa_nt
r1:/work/smaman/ENCOURS/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/SRR1618559_1.fastq
r2:/work/smaman/ENCOURS/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/SRR1618559_2.fastq
flagstat:/work/smaman/ENCOURS/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/hisat2_SRR1618559.sort.flagstat.txt
candidat:/work/smaman/ENCOURS/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/srs_SRR1618559_candidates.txt
split:/work/smaman/ENCOURS/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/srs_SRR1618559_split_alignments.txt
prefix:Azolla_f_SRR1618559
threads: 6                
output:/work/smaman/ENCOURS/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/HTML/
force: -F
```


```diff
#Load modules
module load devel/Miniconda/Miniconda3
source activate ISeeker_environment

cd /work/user/smaman/23032022_intronSeeker/intronSeeker/data/REAL_DATA/Azolla_f_SRR1618559/;
python3 /work/smaman/intronSeeker/scripts/simulation2HTML.py -F --config_file  SRR1618559.cfg;

```

HTML report is available in this directory.
