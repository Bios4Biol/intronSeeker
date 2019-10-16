#!/usr/bin/env python3

import textwrap
import configparser
import os

config = configparser.RawConfigParser()
config.read(os.path.abspath(os.path.dirname(__file__))+"/../config/intronSeeker.properties")

# Metadata readable in intronSeeker/config/intronSeeker.properties
__author__ = config['Metadata']['author']
__copyright__ = config['Metadata']['copyright']
__license__ = config['Metadata']['license']
__version__ = config['Metadata']['version']
__email__ = config['Metadata']['email']
__status__ = config['Metadata']['status']

########################################
########### Version Printing ###########
########################################
def program_version() :
    print("intronSeeker v{version}".format(version=__version__),end="\n\n")

########################################
##### Global program help printing #####
########################################

def program_help() :
    text = '\
\n\
Program: intronSeeker\n\
Version: v{version}\n\n\
\
This tool identify  potentially retained introns  in de novo RNA-seq  assembly\
 in order to quantify and remove them. The intron detection is based on the\
 read splicing signal  which appears when an alignement is performed beetween\
 the reads library and the assembly.\n\
This program makes available  all the tools  to produce (STAR and Hisat2 aligners),\
 detect (reads splicing  identification)  and analyze  (ORF prediction,  protein\
 alignement,  data integration)  the signal  in  order to  give  the most  reliable\
 and consistent results as possible in term of retained intron detection.\n\
\n\
You can find also two types of RNA-seq data simulations to validate  the\
 detection process and measure the false positive detection rate.\n\
The first  simulation module uses random sequence  simulation in order  to check\
 if splice aligners are  able  to find  inserted  introns  when only  contigs  with\
 introns and reads without intron are used as well as when contigs with and\
 without introns and reads without introns are used.\n\
The second simulation is based on an existing genome and corresponding genome\
 annotation. In this case the simulator produces reads with and without intron\
 or spliced exons as well as transcripts with and whithout introns or spliced\
 exon.  This modules enables to  verify the  fraction  of retained  introns  which\
 can be detected in real condition  and  set the appropriate detection\
 thresholds.\n\
 \n'.format(version=__version__)

    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
    )
    
    cw = textwrap.TextWrapper(
        width=60,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    # Program Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]))
    
    # Usage
    print('Usage : intronSeeker <command> [arguments] [--help] [--version]')
    print('(To know the detailed usage of each sub-commands use \'intronSeeker <command> --help\')',end='\n\n')
    
    # Detail of the commands
    print('Commands : ')
    
    print(' -- Align')
    # hisat2Alignement
    print('   hisat2Alignement',end='')
    print(cw.fill(
'Calls Hisat2 to align the reads library and the reference assembly and produce\
 a BAM file on which reads split detection can be performed.'
    ))
    # starAlignement
    print('   starAlignement',end='')
    print(cw.fill(
'Calls STAR to align the reads library and the reference assembly and produce\
 a BAM file on which reads splicing detection can be performed.'
    ))
    print()
    
    print(' -- Analyze')
    # splitReadSearch
    print('   splitReadSearch',end='')
    print(cw.fill(
'Performs the extraction of read splicing detection on a BAM file in order to\
 identify and caracterize introns candidates.'
    ))
    # trimFastaFromGFF
    print('   trimFastaFromGFF',end='')
    print(cw.fill(
'From a reference Fasta file and a associated GFF, produces a new Fasta file\
 containing each sequence of the original Fasta spliced of the segments\
 corresponding to the features of the GFF file.'
    ))
    # analyzeORF
    print('   analyzeORF\t',end='')
    print(cw.fill(
'From a reference Fasta file, predicts the Open Reading Frames for each\
 sequence of the Fasta and returns  the characteristics of the best predicted\
 ORF for each sequence in a tabulated file. The ORF prediction is performed\
 with TransDecoder.'
    ))
    # analyzeProtein
    print('   analyzeProtein',end='')
    print(cw.fill(
'From a reference Fasta file and a proteic bank, performs the proteic alignement\
 all-vs-all beetween the Fasta and the bank and returns the characteristics\
 of the best alignement for each sequence of the Fasta in a tabulated file.\
 The proteic alignement is performed with Diamond.'
    ))
    print()
    
    print(' -- Simulation')
    # fullRandomSimulation
    print('   fullRandomSimulation',end='')
    print(cw.fill(
'Simulates a multi Fasta file which contains the pseudo-contigs and a TXT file\
 which contains the pseudo retained introns charesteristics. All the sequences\
 (contigs and introns) as well as the introns insertion are fully random.'
    ))
    # fullRandomSimulation
    print('   GTFbasedSimulation',end='')
    print(cw.fill(
'From a genome and an associated GTF file, generates pseudo-contigs (corresponds\
 to transcripts) with potentially retained introns or spliced exons. Three\
 files are generated : one corresponds to reference pseudo-assembly, another\
 one is used for dimulate the reads library and the third one gathers all\
 the simulated features (retained intorns or spliced exons). gffread program\
 is called during the simulation process.'
    ))
    # simulateReads
    print('   simulateReads',end='')
    print(cw.fill(
'From a Fasta file, calls Grinder program to simulates a corresponding reads library.'
    ))
    print()
    
    print(' -- Test')
    # checkInstall
    print('   checkInstall\t',end='')
    print(cw.fill(
'Checks the correct installation of the dependencies as well as the dependencies\' versions.'
    ))
    
    
    print()
    print('Program:   intronSeeker version {version}')
    print('Version:   {version}'.format(version=__version__))
    print('License:   {license}'.format(license=__license__))
    print('Copyright: {copyright}'.format(copyright=__copyright__))
    print('Authors:   {author}'.format(author=__author__))
    print('Support:   {email}\n'.format(email=__email__))


########################################
####### Command help dispatching #######
########################################

def command_help(command: str) :
    if command == "starAlignement" :
        star_help()
    elif command == "hisat2Alignement" :
        hisat2_help()
    elif command == "splitReadSearch" :
        split_help()
    elif command == "trimFastaFromTXT" :
        trim_help()
    elif command == "analyzeORF" :
        orf_help()
    elif command == "analyzeProtein" :
        protein_help()
    
########################################
##### starAlignement help printing #####
########################################

def star_help() :
    text='\
    \n\
starAlignement :\n\n\
Calls STAR to align the reads library and the reference assembly and produce\
 a BAM file on which reads splicing detection can be performed.\
 A directory <output>_star is created where the STAR indexing, STAR log and output BAM and BAI file will be stored\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
        drop_whitespace=True
    )

    # starAlignement Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print(' Usage:')
    print(
    textwrap.fill('\
 intronSeeker starAlignement -r <reference.fa> -1 <read1.fq>\
 -o <outfile_basename> [-2 <read2.fq>] [-t <num_threads>]',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=55,
        initial_indent="\t\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Arguments :')
    print('   -r/--reference',end='')
    print(cw.fill(
'Name of the reference FASTA file on wich the reads will be mapped.'
    ))
    print('   -1/--r1\t',end='')
    print(cw.fill(
'Name of the FASTQ file which contains the single-end reads library. If paired-end,\
 filename of #1 reads mates.'
    ))
    print('   -2/--r2\t',end='')
    print(cw.fill(
'Only for a paired-end library, filename of #2 reads mates '
    ))
    print('   -o/--output\t',end='')
    print(cw.fill(
'Basename of the output directory where the STAR index and ouput BAM fil will be stored.'
    ))
    print('   -t/--threads\t',end='')
    print(cw.fill(
'Number of alignment threads to lauch.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()

########################################
#### hisat2Alignement help printing ####
########################################

def hisat2_help() :
    text='\
    \n\
hisat2Alignement :\n\n\
Calls HiSat2 to align the reads library and the reference assembly and produce\
 a BAM file on which reads splicing detection can be performed.\
 A directory <output>_hisat2 is created where the HiSat2 indexing, HiSat2 log and output BAM and BAI file will be stored\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
        drop_whitespace=True
    )

    # hisat2Alignement Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print(' Usage:')
    print(
    textwrap.fill('\
 intronSeeker hisat2Alignement -r <reference.fa> -1 <read1.fq>\
 -o <outfile_basename> [-2 <read2.fq>] [-t <num_threads>]',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=55,
        initial_indent="\t\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Arguments :')
    print('   -r/--reference',end='')
    print(cw.fill(
'Name of the reference FASTA file on wich the reads will be mapped.'
    ))
    print('   -1/--r1\t',end='')
    print(cw.fill(
'Name of the FASTQ file which contains the single-end reads library. If paired-end,\
 filename of #1 reads mates.'
    ))
    print('   -2/--r2\t',end='')
    print(cw.fill(
'Only for a paired-end library, filename of #2 reads mates '
    ))
    print('   -o/--output\t',end='')
    print(cw.fill(
'Basename of the output directory where the HiSat2 index and ouput BAM fil will be stored.'
    ))
    print('   -t/--threads\t',end='')
    print(cw.fill(
'Number of alignment threads to lauch.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()

########################################
#### splitReadSearch help printing #####
########################################

def split_help() :
    text='\
    \n\
splitReadSearch :\n\n\
Performs the split read signal extraction from an alignement (BAM File) beetween an assembly (FASTA file) and its reads library.\
 A split read is interpreted like a potential splicing event and, so, a candidate for the intron identification.\
 Several splicing events could be merged in one entity if their splicing limits are the same (+/-2b).\
 A directory named <output>_splitReadSearch is created where the TXT output file is stored. Each line of this file characterize an intron candidate.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
        drop_whitespace=True
    )

    # splitReadSearch Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print(' Usage:')
    print(
    textwrap.fill('\
 intronSeeker splitReadSearch -a <alignement.bam> -r <reference.fa> \
 -o <outfile_basename> ',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=55,
        initial_indent="\t\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Arguments :')
    print('   -a/--alignement',end='')
    print(cw.fill(
'Name of the alifnement BAM file from which the split read signal has to be extracted.'
    ))
    print('   -r/--reference',end='')
    print(cw.fill(
'Name of the reference FASTA file of the alignement and in which the introns have to be identified.'
    ))
    print('   -o/--output\t',end='')
    print(cw.fill(
'Basename of the output directory where the TXT ouput file will be stored.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()

########################################
#### trimFastaFromTXT help printing ####
########################################

def trim_help() :
    text='\
    \n\
trimFastaFromTXT :\n\n\
From a tabulated TXT which contains all the features to trim and the reference FASTA where the features have to be trimmed,\
 produces a new FASTA reference file where the features were trimmed. The TXT file must hold at least 4 necessary fields :\n\
-reference : ID of the sequence to splice\n\
-start : feature\'s starting coordinate on the reference\n\
-end : ending\'s coordinate on the reference\n\
-to_trim : 0 or 1 value, 1 if the feature must be trimmed, 0 otherwise\n\n\
The TXT file can hold other fields but they will not influence trimming process.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
        drop_whitespace=True
    )

    # trimFastaFromTXT Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print(' Usage:')
    print(
    textwrap.fill('\
 intronSeeker trimFastaFromTXT -r <reference.fa> -f <features.txt>\
 -o <outfile_basename> ',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=55,
        initial_indent="\t\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Arguments :')
    print('   -r/--reference',end='')
    print(cw.fill(
'Name of the reference FASTA file where the features must be spliced.'
    ))
    print('   -f/--features',end='')
    print(cw.fill(
'Name of the features TXT file which has to hold at least the \'reference\', \'start\', \'end\' and \'to_trim\' fields (above-detailed).'
    ))
    print('   -o/--output\t',end='')
    print(cw.fill(
'Basename of the output directory where the FASTA ouput file will be stored.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()


########################################
####### analyzeORF help printing #######
########################################

def orf_help() :
    text='\
    \n\
analyzeORF :\n\n\
Performs the ORF analysis on a reference FASTA file : predicts, with TransDecoder software, the ORF of each sequence of the FASTA\
 and produces a file where each line characterizes the best predicted ORF of each sequence. This module is involved two times in intorn identification :\
 first, directly on the original reference file to predict ORFs and a second time after the intron candidates identification :\
 the candidates are trimmed from the original reference FASTA and the ORFs prediction is performed on the resulting file. The goal of this approach\
 is to analyze the potential ORFs modifications in order to help us in intorn identification.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
        drop_whitespace=True
    )

    # analyzeORF Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print(' Usage:')
    print(
    textwrap.fill('\
 intronSeeker analyzeORF -r <reference.fa> -o <outfile_basename>\
 [-k] [--no-refine-starts]',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=55,
        initial_indent="\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Arguments :')
    print('   -r/--reference',end='')
    print(cw.fill(
'\tName of the reference FASTA file on which the ORF prediction must be performed.'
    ))
    print('   -o/--output\t',end='')
    print(cw.fill(
'\tBasename of the output directory where the ouput file and potential TransDEcoder\'s files will be stored.'
    ))
    print('   -k/--keep-intermediate',end='')
    print(cw.fill(
'Set up this option, if you want to keep the TransDecoder\'s intermediate files.'
    ))
    print('   --no-refine-starts',end='')
    print(cw.fill(
'\tAllow TransDecoder to not perform ORF\'s starts refining. Use this option only if TransDecoder fails at this step of ORF prediction.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'\tPrint this help message.'
    ))
    print()


########################################
##### analyzeProtein help printing #####
########################################

def protein_help() :
    text='\
    \n\
analyzeProtein :\n\n\
Performs a proteic alignement on a reference FASTA file with Diamond. Produces a file where each line characterize the best alignement for each sequence of the FASTA.\
 This module take part of the same approach of the analyzeORF module : the proteic alignements before and after candidate trimming are compared to support the intron identification.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
        drop_whitespace=True
    )

    # analyzeProtein Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print(' Usage:')
    print(
    textwrap.fill('\
 intronSeeker analyzeProtein -r <reference.fa> -o <outfile_basename>\
 [-k] [--no-refine-starts]',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=55,
        initial_indent="\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Arguments :')
    print('   -r/--reference',end='')
    print(cw.fill(
'\tName of the reference FASTA file on which theproteic alignement must be performed.'
    ))
    print('   -p/--db-proteins',end='')
    print(cw.fill(
'\t.'
    ))
    print('   -o/--output\t',end='')
    print(cw.fill(
'\tBasename of the output directory where the ouput file and potential TransDEcoder\'s files will be stored.'
    ))
    print('   -k/--keep-intermediate',end='')
    print(cw.fill(
'Set up this option, if you want to keep the TransDecoder\'s intermediate files.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'\tPrint this help message.'
    ))
    print()
