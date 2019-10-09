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

def program_version() :
    print("intronSeeker v{version}".format(version=__version__),end="\n\n")

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
 \n'.format(version=__version__)

    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
    )
    
    cw = textwrap.TextWrapper(
        width=60,
        initial_indent="\t\t",
        subsequent_indent="\t\t\t\t",
        break_long_words=False
    )
    # Program Description
    print("\n".join([tw.fill(line) for line in text.splitlines()]))
    
    # Usage
    print('Usage : intronSeeker <command> [arguments] [--help] [--version]')
    print('(To know the detailed usage of each sub-commands use \'intronSeeker <command> --help\')',end='\n\n')
    
    # Detail of the commands
    print('Commands : ')
    
    print(' --Align')
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
    
    print(' --Analyze')
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
    
    print(' --Simulation')
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
    
    print(' --Test')
    # checkInstall
    print('   checkInstall\t',end='')
    print(cw.fill(
'Checks the correct installation of the dependencies as well as the dependencies\' versions.'
    ))
    
    
    print()
    print('Program: intronSeeker version {version}'.format(version=__version__))
    print('License : {license}'.format(license=__license__))
    print('Copyright : {copyright}'.format(copyright=__copyright__))
    print('Author(s): {author}'.format(author=__author__))
    print('Support : {email}'.format(email=__email__))


def command_help(command: str) :
    if command == "starAlignement" :
        star_help()
    

def star_help() :
    text='\
\n\
intronSeeker starAlignement :\n\n\
Calls STAR to align the reads library and the reference assembly and produce\
 a BAM file on which reads splicing detection can be performed.\
'
    tw = textwrap.TextWrapper(
        width=90,
        initial_indent="    ",
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
    
    
