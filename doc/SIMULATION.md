intronSeeker's data simulation module's documentation
=====================================================


How it works ?
--------------

#### Full-random simulation 

To simulate Full-random RNA-seq data, a module is implemeted. Hereafter, we'll review the 
different step performed to generate this type of data :

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

Here, we'll present a non-exhaustive documentation on how to use the intronSeeker's 
data simulation features. For an exhaustive documentation about the different functions 
(description, options...), read the HOW TO USE file in doc directory. 
We'll just present a fast and basic usage with examples of the program. All the 
input files or output files presented here are available in data directory. 

Before running any intronSeeker command, activate the conda environment with :

```diff
+ source activate Stalker_environment
```

#### Full-random data simulation

This program offer a full-random data simulaton functionality. It is possible to 
generate a set of random contigs (their sequences are totally random), insert random introns inside contigs
and generate a library of reads from these contigs. 

##### Contigs generation

To generate the set of contigs **without** introns, run the command : 
```diff
+ intronStalker contig -n 10 --min 150 --max 1000 --output OutputRandomContig.fa
```
Or just run :

```diff
+ intronStalker contig
```
The two commands above give the same results : the arguments given in the first 
command are the defaults values. The `-n` argument rules the total number of 
sequences (i.e. contigs), `--min` is the minimal length of the contigs and `--max`
the maximal length. Finally, `--output` (or `-o`) correspond to output filename.
Like all of this program's features, the result will be generate in a separate 
directory. Here, the results will be store in `OutputRandomContig_contig`.

##### Introns insertion 

Now, to insert random introns in all newly generated contigs, you have to run :

 ```diff
 + intronStalker intron -i OutputRandomContig_contig/OutputRandomContig.fa
 ``` 

This command insert random introns in each contig. If you don't want intron in 
each contig, you can use the `--rand` option : introns will be inserted in only 
half of the contigs (randomly choosen). For the details of default values and 
available parameters, see the HOW TO USE file in doc directory.
This command creates three files but two are interesting : the fasta file which 
contains the contigs with inserted introns and the `Coord.txt` file which contains 
all the information about introns (contig name, begin/end coordinates, strand).

##### Library generation : Grinder running

Now, we get the contigs and the introns, we have to generate the reads library. 
The library is generated from reference contigs (i.e. without introns). Grinder 
needs a lot of parameter (detailed in Grinder documentation and in Grinder section
in HOW TO USE file in doc repertory). An example of basic grinder.cfg file (i.e. the file 
which conatains all grinder parameters) is available in config directory. The command is : 

```diff
+ intronStalker simulatedReads -r ref.fa -p grinder.cfg -o outputDirName
```

This command generates two or three files : the first one is the `ranks.txt` file
which contains information about the library (abundance percentage of each contig 
in term of reads). Other files are archive which contains reads : there is only 
one file for single-end library or two for paired-end library.