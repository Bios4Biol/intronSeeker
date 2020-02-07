intronSeeker's data simulation 


How it works ?
--------------

#### Full-random simulation 

To simulate Full-random RNA-seq data, several commands are implemeted :

1. Generate a FASTA file which contains random sequences. Simulation parameters are number of sequence(s),
maximum and minimum lengths of generated sequences.

2. Insert random introns in a FASTA file (i.e. the generated random contigs). 
This produces a FASTA file with contigs containing introns and a TEXT 
file containing inserted introns' description and location information.

3. Use [Grinder](https://sourceforge.net/projects/biogrinder/) program to 
generate an intron free simulated contigs read library file.

#### Data Simulation from a reference genome


How to use ? 
------------

This documentation, not exhaustive, explain intronSeeker's data simulation features. 
For an exhaustive documentation (functions description, options...), 
please see HOW TO USE file in doc/ directory. 
This is a fast and basic usage documentation with program examples. All  
input or output files presented here are availables in data directory. 

Before running any intronSeeker command, activate ISeeker environment with :

```diff
+ source activate ISeeker_environment;
```

#### Full-random data simulation and simulate reads

This program offer a full-random data simulaton functionality. This generate a set of 
random contigs (totally random sequences), insert random introns inside contigs
and generate a library of reads from these contigs. 

Reads simulation is done in 2 steps. First, run intronSeeker fullRandomSimulation (for FRS data) 
or intronSeeker GTFbasedSimulation (for GBS data) in order to write a FASTA file with contigs (*-modified.fa file). 
Then, run intronSeeker simulateReads to write sequences FASTQ files (R1 and R2) from your reference FASTA file.

##### FRS data

###### Step 1 : intronSeeker fullRandomSimulation

FASTA file with contigs result of intronSeeker fullRandomSimulation.

```diff
+ intronSeeker fullRandomSimulation -o DirName -p FileNameBeginning;
```


This command insert random introns in each contig. 
`-r` option : introns will be inserted in only 
half of the contigs (randomly choosen). 
For more informations (default values, parameters), please
see the HOW TO USE file in doc directory.

Outputs : 
A FASTA file with contigs.
A modified FASTA file with modified contigs with inserted introns.
A GTF file, description of introns (contig name, begin/end coordinates, strand).
Results are generated in a separate directory.

###### Step 2  :  intronSeeker simulateReads 

Reads results of intronSeeker simulateReads script, 
runned from reference fasta file (not from the modified FASTA file). 


```diff
+ intronSeeker simulateReads -o DirName -p FileNameBeginning -f frs*contigs.fa -c grinderFile; 
```

The library is generated from reference contigs (i.e. without introns). Grinder 
needs severals parameters (detailed in Grinder documentation and in Grinder section
in HOW TO USE file in doc/ repertory). An example of basic grinder.cfg file (i.e. the file 
which contains all grinder parameters) is available in config/ directory.

Outputs:
A TEXT file with library inormation (abundance percentage of each contig 
in term of reads).
One FASTQ file for single-end library or two R1 and R2 FASTQ files for paired-end library.


##### GBS data

###### Step 1 : intronSeeker GTFbasedSimulation

FASTA file with contigs result of intronSeeker GTFbasedSimulation.

 ```diff
 + intronSeeker GTFbasedSimulation -a gtfFile -r RefFastaFile -o DirName -p FileNameBeginning;
 ``` 
 
###### Step 2 Generate the reads library

Reads result of intronSeeker simulateReads script, 
runned from reference fasta file (not from the modified FASTA file)

 ```diff
 + intronSeeker simulateReads -f RefFastaFile -c $grinderFile -o DirName -p FileNameBeginning;
 ``` 


