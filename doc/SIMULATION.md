intronSeeker's data simulation 


How it works ?
--------------

#### Full-random simulation 

To simulate Full-random RNA-seq data, several commands are implemeted :

1. Generate a fasta file which contains random sequences. The simulation parameters are number of sequence,
maximum and minimum lengths of the generated sequences.

2. Insert random introns in a multi-fasta file (i.e. the generated random contigs). 
This produces a multi-fasta file with the contigs containing introns and a text 
file containing inserted introns' description and location information.

3. Use [Grinder](https://sourceforge.net/projects/biogrinder/) program to 
generate an intron free simulated contigs read library file.

#### Data Simulation from a reference genome


How to use ? 
------------

This documentation, not exhaustive, explain intronSeeker's data simulation features. 
For an exhaustive documentation about the different functions 
(description, options...), read the HOW TO USE file in doc directory. 
This is a fast and basic usage documentation with program examples. All  
input or output files presented here are availables in data directory. 

Before running any intronSeeker command, activate ISeeker environment with :

```diff
+ source activate ISeeker_environment;
```

#### Full-random data simulation and simulate reads

This program offer a full-random data simulaton functionality. It is possible to 
generate a set of random contigs (their sequences are totally random), insert random introns inside contigs
and generate a library of reads from these contigs. 

Reads simulation is done in 2 steps. First, run intronSeeker fullRandomSimulation (for FRS data) 
or intronSeeker GTFbasedSimulation (for GBS data) in order to write a fasta file with contigs (*-modified.fa file). 
Then, run intronSeeker simulateReads to write sequences fastq files (R1 and R2) from your reference fasta file.

##### FRS data

###### Step 1 : intronSeeker fullRandomSimulation

Fasta file with contigs result of intronSeeker fullRandomSimulation

```diff
+ intronSeeker fullRandomSimulation -o DirName -p FileNameBeginning;
```


Results generated in a separate directory.
This command insert random introns in each contig. 
`-r` option : introns will be inserted in only 
half of the contigs (randomly choosen). 
For the details of default values and available parameters, 
see the HOW TO USE file in doc directory.
Outputs : 
A FASTA file which contains contigs.
A modified FASTA file which contains contigs with inserted introns.
A GTF file which contains all the information about introns 
(contig name, begin/end coordinates, strand).


###### Step 2  :  intronSeeker simulateReads 

Reads results of intronSeeker simulateReads script, 
runned from reference fasta file (no modified) 


```diff
+ intronSeeker simulateReads -o DirName -p FileNameBeginning -f frs*contigs.fa -c grinderFile; 
```

The library is generated from reference contigs (i.e. without introns). Grinder 
needs a lot of parameters (detailed in Grinder documentation and in Grinder section
in HOW TO USE file in doc repertory). An example of basic grinder.cfg file (i.e. the file 
which contains all grinder parameters) is available in config directory.

Outputs:
A txt file which contains information about the library (abundance percentage of each contig 
in term of reads).
Other files are archive which contains reads : there is only 
one file for single-end library or two for paired-end library.


##### GBS data

###### Step 1 : intronSeeker GTFbasedSimulation

Fasta file with contigs result of intronSeeker GTFbasedSimulation

 ```diff
 + intronSeeker GTFbasedSimulation -a gtfFile -r RefFastaFile -o DirName -p FileNameBeginning;
 ``` 
 
###### Step 2 Generate the reads library

Reads results of intronSeeker simulateReads script, 
runned from reference fasta file (no modified)

 ```diff
 + intronSeeker simulateReads -f RefFastaFile -c $grinderFile -o DirName -p FileNameBeginning;
 ``` 


