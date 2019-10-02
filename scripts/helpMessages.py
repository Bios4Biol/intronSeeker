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

def program_help() :
    text = '\
\n\
intronSeeker v{version}\n\
-------------------\n\n\
This tool identify potentially retained introns in de novo RNA-seq assembly\
 in order to quantify and remove them. The intron detection is based on the\
 read splicing signal which appears when an alignement is performed beetween\
 the reads library and the assembly.\n\
This program makes available all the tools to produce (STAR and Hisat2 aligners),\
 detect (reads splicing identification) and analyze (ORF prediction, protein\
 alignement, data integration) the signal in order to give the most reliable\
 and consistent results as possible in term of retained intron detection.\n\
\n\
You can find also two types of RNA-seq data simulations to validate the\
 detection process and measure the false positive detection rate.\n\
The first simulation module uses random sequence simulation in order to check\
 if splice aligners are able to find inserted introns when only contigs with\
 introns and reads without intron are used as well as when contigs with and\
 without introns and reads without introns are used.\n\
The second simulation is based on an existing genome and corresponding genome\
 annotation. In this case the simulator produces reads with an without intron\
 or spliced exons as well as transcripts with and whithout introns or spliced\
 exon. This modules enables to verify the fraction of retained introns which\
 can be detected in real condition and and set the appropriate detection\
 thresholds.\n\
 \n\n'.format(version=__version__)

    tw = textwrap.TextWrapper(
        width=80,
        initial_indent="        ",
        subsequent_indent="    "
    )
    
    cw = textwrap.TextWrapper(
        width=45,
        initial_indent="\t\t",
        subsequent_indent="\t\t\t\t\t"
    
    )
    # Program Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]))
    
    # Usage
    print('\tUsage : intronSeeker <command> [options]\n')
    
    # Detail of the commands
    print('\tCommands : ')
    
    print('\t --Align')
    # hisat2Alignement
    print('\t   hisat2Alignement',end='')
    print(cw.fill(
'Call Hisat2 to align the reads library and the reference assembly and produce\
 a BAM file on which reads split detection can be performed.'
    ))
    # starAlignement
    print('\t   star2Alignement',end='')
    print(cw.fill(
'Call STAR to align the reads library and the reference assembly and produce\
 a BAM file on which reads splicing detection can be performed.'
    ))
    print()
    
    print('\t --Analyze')
    # splitReadSearch
    print('\t   splitReadSearch',end='')
    print(cw.fill(
'Perform the extraction of read splicing detection on a BAM file in order to\
 identify and caracterize introns candidates.'
    ))
    # trimFastaFromGFF
    print('\t   trimFastaFromGFF',end='')
    print(cw.fill(
'From a reference Fasta file and a associated GFF, produce a new Fasta file\
 containing each sequence of the original Fasta spliced of the segments\
 corresponding to the features of the GFF file.'
    ))
    # analyzeORF
    print('\t   analyzeORF\t',end='')
    print(cw.fill(
'From a reference Fasta file, predicts the Open Reading Frames for each\
 sequence of the Fasta and returns  the characteristics of the best predicted\
 ORF for each sequence in a tabulated file. The ORF prediction is performed\
 with TransDecoder.'
    ))
    # analyzeProtein
    print('\t   analyzeProtein',end='')
    print(cw.fill(
'From a reference Fasta file and a proteic bank, perform the proteic alignement\
 all-vs-all beetween the Fasta and the bank and returns the characteristics\
 of the best alignement for each sequence of the Fasta in a tabulated file.\
 The proteic alignement is performed with Diamond.'
    ))
    print()
    
    print('\t --Simulation')
    # fullRandomSimulation
    print('\t   fullRandomSimulation',end='')
    print(cw.fill(
'TO FILL'
    ))
    # fullRandomSimulation
    print('\t   GTFbasedSimulation',end='')
    print(cw.fill(
'TO FILL'
    ))
    # simulateReads
    print('\t   simulateReads',end='')
    print(cw.fill(
'TO FILL'
    ))
    print()
    
    print('\t --Test')
    # checkInstall
    print('\t   checkInstall\t',end='')
    print(cw.fill(
'TO FILL'
    ))
    print()
    

def command_help(command: str ,verbosity: int) :
    pass
     # ~ if command == "starAlignement" :
         # ~ if verbosity == 0 :
            # ~ print(command)
    # ~ elif command == "hisat2Alignement" :
        # ~ if verbosity == 0 :
            
