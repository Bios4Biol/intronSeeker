#!/usr/bin/env python3

import textwrap
import os
import sys

# Metadata
__author__    = "Sarah Maman, Philippe Bardou, Emilien Lasguignes, Faustine Oudin, Floréal Cabanettes, Christophe Klopp."
__copyright__ = "Copyright (C) 2019 INRAE"
__license__   = "GNU General Public License"
__version__   = "1.0"
__email__     = "sigenae(at)inrae.fr"
__status__    = "dev"


########################################
###########   Logs Printing  ###########
########################################
def print_to_stdout(*a):
 
    # Here a is the array holding the objects
    # passed as the argument of the function
    print(*a, file = sys.stdout)

########################################
########### Version Printing ###########
########################################
def program_version() :
    print("intronSeeker {version}".format(version=__version__),end="\n")
    print(__copyright__) 

########################################
##### Global program help printing #####
########################################

def program_help() :
    text = '\
\n\
Program: intronSeeker\n\
Version: {version}\n\n\
\
This tool identify  potentially retained introns  in de novo RNA-seq  assembly\
 in order to quantify and remove them. The intron detection is based on the\
 read splicing signal  which appears when  an alignment is performed beetween\
 the reads library and the assembly.\n\
This program makes available  all the tools  to produce (STAR and Hisat2 aligners),\
 detect (reads splicing  identification)  and analyze  (ORF prediction,  protein\
 alignment,  data  integration)  the signal  in  order to  give  the most  reliable\
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
    print('Usage: intronSeeker <command> [arguments] [--help] [--version]')
    print('(To know the detailed usage of each sub-commands use \'intronSeeker <command> --help\')',end='\n\n')
    
    # Detail of the commands
    print('Commands: ')
    
    print(' -- Align')
    # hisat2Alignment
    print('   hisat2Alignment',end='')
    print(cw.fill(
'Calls Hisat2 to align the reads library and the reference assembly  and  produce\
 a BAM file on  which  reads  split detection can be performed.'
    ))
    # starAlignment
    print('   starAlignment',end='')
    print(cw.fill(
'Calls STAR to align  the reads library and  the reference assembly and produce\
 a BAM file  on which  reads splicing detection can be performed.'
    ))
    print()
    
    print(' -- Analyze')
    # splitReadSearch
    print('   splitReadSearch',end='')
    print(cw.fill(
'Performs the extraction of read splicing detection on BAM file in order to\
 identify intron candidates.'
    ))
    # trimFastaFromGFF
    print('   trimFastaFromTXT',end='')
    print(cw.fill(
'From a reference Fasta file and a associated TXT produces  a new Fasta file\
 containing each sequence of the original Fasta  spliced  of  the  segments \
 corresponding  to  the features of the TXT file.'
    ))
    # findEvidence
    print('   findEvidence\t',end='')
    print(cw.fill(
'From the reference  and trimmed Fasta files, predicts the Open  Reading  Frames\
  (getorf)  and  find  proteic  hits (diamond) to find evidences.'
    ))
    print()
    
    print(' -- Simulation')
    # fullRandomSimulation
    print('   fullRandomSimulation',end='')
    print(cw.fill(
'Simulates two  Fasta  files which  contains  the  pseudo-contigs and a GTF file\
 which contains the pseudo retained  introns  charesteristics.  All the sequences\
 (contigs and introns)  as  well  as  the  introns  insertion are fully random.'
    ))
    # GTFbasedSimulation
    print('   GTFbasedSimulation',end='')
    print(cw.fill(
'From  a genome  and an  associated  GTF  file,  generates   pseudo-contigs   (corresponds \
  to   transcripts)    with potentially  retained  introns  or  spliced exons.  Three\
 files  are  generated:   one  corresponds  to   reference pseudo-assembly, another\
 one  is used  for  simulate  the reads library and the third one gathers all\
 the simulated features  (retained  introns  or  spliced exons). gffread program\
 is called during the simulation process.'
    ))
    # simulateReads
    print('   simulateReads',end='')
    print(cw.fill(
'From a Fasta file,  calls Grinder program  to simulate a corresponding reads library.'
    ))
    print()
    
    print(' -- Test')
    # checkInstall
    print('   checkInstall\t',end='')
    print(cw.fill(
'Checks the correct  installation  of the dependencies and versions.'
    ))
    
    
    print()
    print('Program:   intronSeeker')
    print('Version:   {version}'.format(version=__version__))
    print('License:   {license}'.format(license=__license__))
    print('Copyright: {copyright}'.format(copyright=__copyright__))
    print('Authors:   {author}'.format(author=__author__))
    print('Support:   {email}\n'.format(email=__email__))


########################################
####### Command help dispatching #######
########################################

def command_help(command: str) :
    if command == "starAlignment" :
        star_help()
    elif command == "hisat2Alignment" :
        hisat2_help()
    elif command == "splitReadSearch" :
        split_help()
    elif command == "trimFastaFromTXT" :
        trim_help()
    elif command == "findEvidence" :
        findEvidence_help()
    elif command == "fullRandomSimulation" :
        frs_help()
    elif command == "GTFbasedSimulation" :
        gbs_help()
    elif command == "simulateReads" :
        sr_help()
    elif command == "checkInstall" :
        checki_help()

########################################
##### starAlignement help printing #####
########################################

def star_help() :
    text='\
\nDescription:\n\
Calls STAR to align the reads library and the reference assembly and produce\
 a BAM file on which  reads  splicing  detection  can be performed.'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # starAlignement Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print('intronSeeker starAlignment -r <reference.fa> -1 <reads1.fq> -o STR')
    print('\t\t\t   [-2 <reads2.fq>] [-F, -k] [-t INT] [-p STR]\n')
    
    cw = textwrap.TextWrapper(
        width=67,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -r/--reference FILE',end='')
    print(cw.fill(
'Name of the reference FASTA file on wich the reads will be mapped.'
    ))
    print('   -1/--r1 FILE\t',end='')
    print(cw.fill(
'Name of the  FASTQ  file  which  contains  the  single-end   reads library. If paired-end,\
 filename of #1 reads mates.'
    ))
    print('   -2/--r2 FILE\t',end='')
    print(cw.fill(
'Only for a paired-end library, filename of #2 reads mates.'
    ))
    print('   -p/--prefix STR',end='')
    print(cw.fill(
'Prefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'Output dir name.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\tForce to overwrite output files.'
    ))
    print('   -k/--keep',end='')
    print(cw.fill(
'\tKeep all intermediate files.'
    ))
    print('   -t/--threads INT',end='')
    print(cw.fill(
'Number of alignment threads to launch [1].'
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
\nDescription:\n\
Calls HiSat2  to align the reads library and the reference assembly and produce\
 a BAM file on which reads splicing detection can be performed.'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # hisat2Alignement Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print('intronSeeker hisat2Alignment -r <reference.fa> -1 <reads1.fq> -o STR')
    print('\t\t\t     [-2 <reads2.fq>] [-t INT] [-F, -p STR]')
    
    cw = textwrap.TextWrapper(
        width=67,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -r/--reference FILE',end='')
    print(cw.fill(
'Name of the reference FASTA file on wich the reads will be mapped.'
    ))
    print('   -1/--r1 FILE\t',end='')
    print(cw.fill(
'Name of the  FASTQ  file  which  contains  the  single-end   reads library. If paired-end,\
 filename of #1 reads mates.'
    ))
    print('   -2/--r2 FILE\t',end='')
    print(cw.fill(
'Only for a paired-end library, filename of #2 reads mates.'
    ))
    print('   -p/--prefix STR',end='')
    print(cw.fill(
'Prefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'Output dir name.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\tForce to overwrite output files.'
    ))
    print('   -t/--threads INT',end='')
    print(cw.fill(
'Number of alignment threads to launch [1].'
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
\nDescription:\n\
Performs  the split read  signal  extraction  from  an   alignment (BAM File)  beetween an assembly (FASTA file) and its reads library.\
 A split read is interpreted  like a potential splicing event and, so, a candidate for the intron identification.\
 Several splicing events could be merged in one entity if their splicing limits are  the same (+/-2b).\
  A directory named <output>_splitReadSearch is created where  the TXT output file is stored.  Each line of this file characterize an intron candidate.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # splitReadSearch Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print('intronSeeker splitReadSearch -a <alignment.bam> -r <reference.fa>')
    print('\t\t\t     -o STR [-F, -p STR, -d INT, -l FLOAT, -t INT]')
	
    cw = textwrap.TextWrapper(
        width=67,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -a/--alignment FILE',end='')
    print(cw.fill(
'Name of the  alignment BAM file  from which the split  read signal has to be extracted.'
    ))
    print('   -r/--reference FILE',end='')
    print(cw.fill(
'Name of the reference FASTA file  of the  alignment  and  in which the introns have to be identified.'
    ))
    print('   -d/--min-depth INT',end='')
    print(cw.fill(
'Candidates will be flagged if not supported by more than [1] read.'
    ))
    print('   -l/--max-length INT',end='')
    print(cw.fill(
'Candidates will be flagged if length greater  than [80] percent of contig length.'
    )) 
    print('   -f/--foot-size INT',end='')
    print(cw.fill(
'Splited reads of candidates will be removed if their minimal foot size is less than [5].'
    ))    
    print('   -p/--prefix STR',end='')
    print(cw.fill(
'Prefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'Output dir name.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\tForce to overwrite output files.'
    ))
    print('   -t/--threads INT',end='')
    print(cw.fill(
'Number of threads [1].'
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
\nDescription:\n\
From a tabulated TXT which contains all the features to trim and the reference FASTA where the features have to be trimmed,\
 produces  a new FASTA  reference file where  the features were trimmed.  The TXT file  must hold at least  4  mandatory fields\
 (supplementary fields will not take into account by the trimming process):\n\
  - reference:\tID of the sequence to splice\n\
  - start:\tFeature\'s starting coordinate on the reference\n\
  - end:\tEnding\'s coordinate on the reference\n\
  - to_trim:\t0 or 1 value, 1 if the feature must be trimmed, 0 otherwise\n'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # trimFastaFromTXT Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print(
    textwrap.fill('\
intronSeeker trimFastaFromTXT -r <ref.fa> -c <candidates.txt> -o STR [-m, -F, -p STR]',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=72,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -r/--reference FILE',end='')
    print(cw.fill(
'Name of the reference FASTA file where features must be spliced.'
    ))
    print('   -c/--candidates FILE',end='')
    print(cw.fill(
'Name of the features TXT file  which has  to  hold at  least the \'reference\', \'start\', \'end\', \'selected\' fields (detailed above).'
    ))
    print('   -p/--prefix STR',end='')
    print(cw.fill(
'Prefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'Output dir name.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\tForce to overwrite output files.'
    ))
    print('   -m/--multi',end='')
    print(cw.fill(
'\tAllow  multiple trimming for contigs with several candidates : only one ouput contig is saved with all trimmed candidates.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()


##########################################
####### findEvidence help printing #######
##########################################

def findEvidence_help() :
    text='\
\nDescription:\n\
Analyze and compare ORF (getorf - emboss) and proteic hits (Diamond) between the reference FASTA file and the trimmed reference FASTA file (output of trimFastaFromTXT command).\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # findEvidence Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print('intronSeeker findEvidence -r <ref.fa> -t <trim.fa> -d <db.fa> -c <cand.txt> -o STR')
    print('                          [-p STR] [-F,-k]')
    
    cw = textwrap.TextWrapper(
        width=59,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -r/--reference FILE',end='')
    print(cw.fill(
'Name of the reference FASTA file.'
    ))
    print('   -t/--trimref FILE',end='')
    print(cw.fill(
'Name of the trimmed reference FASTA file (output of trimFastaFromTXT command).'
    ))
    print('   -d/--database FILE',end='')
    print(cw.fill(
'Proteic FASTA file used by Diamond to perform a proteic alignment against "reference" and "trimref" FASTA files.'
    ))
    print('   -c/--candidates FILE',end='')
    print(cw.fill(
'Name of the candidates TXT file.'
    ))
    print('   -p/--prefix STR',end='')
    print(cw.fill(
'Prefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'Output dir name.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\tForce to overwrite output files.'
    ))
    print('   -k/--keep',end='')
    print(cw.fill(
'\tKeep all intermediate files.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()


################################################
##### Full Random Simulation help printing #####
################################################

def frs_help() :
    text='\
\nDescription:\n\
Simulates two Fasta and one GTF files.\n\
Three files are generated: one corresponds  to reference  pseudo-assembly\
 used to simulate the reads library another one contains the modified transcripts (retained introns)\
 used as reference of the  alignment step  and  the third  one  which contains the pseudo  retained\
 intron charesteristics (GTF).\n\
All the sequences are fully random.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print(
    textwrap.fill('\
intronSeeker fullRandomSimulation [-n INT] [-m,-M INT] [-rF] [-l,-H INT] [-p STR] -o STR',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=63,
        initial_indent="\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -n/--nb_contigs INT',end='')
    print(cw.fill(
'\tNumber of contig/sequence randomly generated. [2000].'
    ))
    print('   -m/--min-contig-len INT',end='')
    print(cw.fill(
'Minimal length of random contigs [250].'
    ))
    print('   -M/--max-contig-len INT',end='')
    print(cw.fill(
'Maximal length of random contigs [1500].'
    ))
    print('   -r/--random-part INT\t',end='')
    print(cw.fill(
'Insert  intron  in all  or a part  (percentage)  of  the simulated contigs [100].'
    ))
    print('   --mix/--mix-state\t',end='')
    print(cw.fill(
'Boolean which rules if the modified contigs is mixed: if    the fasta file contains the transcript in the two states  (with AND without retained introns).'
    ))
    print('   -l/--lower-intron-len INT',end='')
    print(cw.fill(
'Minimal length of random intron [250].'
    ))
    print('   -H/--higher-intron-len INT',end='')
    print(cw.fill(
'Maximal length of random intron [750].'
    ))
    print('   -p/--prefix STR',end='')
    print(cw.fill(
'\tPrefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'\tOutput dir name.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\t\tForce to overwrite output files.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'\tPrint this help message.'
    ))
    print()


############################################
##### GTFbasedSimulation help printing #####
############################################

def gbs_help() :
    text='\
\nDescription:\n\
From a genome and associated ensembl GTF file, generates pseudo-contigs (corresponds  to transcripts) with potentially\
 retained introns or spliced exons.\nThree files are generated: one corresponds to reference pseudo-assembly\
 used to simulate the reads library another one contains the modified transcripts (retained introns and or spliced exons)\
 used as reference of the alignment step  and  the third  one gathers  all the  simulated features \
 (retained introns or  spliced  exons).\ngffread  program  is called during  the simulation process.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print('intronSeeker GTFbasedSimulation -a <annotation.gtf> -r <reference.fa> -o STR')
    print('\t\t\t\t[-N|-n INT] [-m] [-F] [-p STR]\n')
    
    cw = textwrap.TextWrapper(
        width=59,
        initial_indent="\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -a/--annotation FILE',end='')
    print(cw.fill(
'\tGFF filename which contains the genome annotation.'
    ))
    print('   -r/--reference FILE',end='')
    print(cw.fill(
'\tName of the reference FASTA file.'
    )) 
    print('   -n/--nb-transcripts INT',end='')
    print(cw.fill(
'Total number of transcript to simulate.'
    ))
    print('   -N/--all-transcripts\t',end='')
    print(cw.fill(
'All genes from GFF have to be transcripted [DEFAULT].'
    ))
    print('   -m/--mix-state\t',end='')
    print(cw.fill(
'Boolean which rules if the modified contigs is mixed:    if the fasta file contains the transcripts in the two states  (with  AND  without  retained introns  and/or spliced exons).'
    ))
    print('   -u/--uniq-transcript\t',end='')
    print(cw.fill(
'Generate transcripts only from uniq transcript genes.'
    ))
    print('   -p/--prefix STR',end='') 
    print(cw.fill(
'\tPrefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'\tBasename for generated files.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\t\tForce to overwrite output files.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'\tPrint this help message.'
    ))
    print()


#######################################
##### simulateReads help printing #####
#######################################

def sr_help() :
    text='\
\nDescription:\n\
From a Fasta file, calls Grinder program to simulate a corresponding reads library.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print(
    textwrap.fill('\
intronSeeker simulateReads -f <file.fa> -c <grinder.cfg> [-p STR] [-F] -o STR',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=59,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -f/--fasta FILE',end='')
    print(cw.fill(
'Name of the FASTA file.'
    ))
    print('   -c/--cfg FILE',end='')
    print(cw.fill(
'Grinder config file (see config/grinder.cfg for example).'
    ))
    print('   -p/--prefix STR',end='')
    print(cw.fill(
'Prefix for output files name.'
    ))
    print('   -o/--output STR',end='')
    print(cw.fill(
'Output dir name.'
    ))
    print('   -F/--force',end='')
    print(cw.fill(
'\tForce to overwrite output files.'
    ))
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()


########################
##### checkInstall #####
########################

def checki_help() :
    text='\
\nDescription:\n\
Checks the correct installation of the dependencies as well as the dependencies\' versions.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="",
        drop_whitespace=True
    )

    # Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]),end='\n\n')
    
    # Usage
    print('Usage:')
    print(
    textwrap.fill('\
intronSeeker checkInstall',
    width=90
    ))
    print()
    
    cw = textwrap.TextWrapper(
        width=59,
        initial_indent="\t",
        subsequent_indent="\t\t\t",
        break_long_words=False
    )
    print('Options:')
    print('   -h/--help\t',end='')
    print(cw.fill(
'Print this help message.'
    ))
    print()
