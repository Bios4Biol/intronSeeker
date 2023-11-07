#!/usr/bin/env python3

#<IntronSeeker searches introns by splice-realigning reads on contigs.>
#Copyright (C) <2019-2024> INRAE
#<Sarah Maman, Philippe Bardou, Emilien Lasguignes, Faustine Oudin, Floréal Cabanettes, Christophe Klopp>
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Modules : 
import argparse 
import sys 
import os 
import argparse 
from argparse import ArgumentParser
import configparser # To parse parameters file
from intronSearch import findEvidence,trimFastaFromTXT,splitReadSearch
from readsMapping import star,hisat2
from dataSimulation import full_random_simulation,gtf_based_simulation,grinder
from checkInstall import checkInstall
from helpMessages import program_help,command_help,program_version, print_to_stdout
from buildReport import simulationReport


def parse_arguments() :
    
    ### Creation of the Argument Parser 
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-v','--version',action='store_true',required=False,dest='version')
    parser.add_argument('-h','--help',action='store_true',required=False,dest='help')
    
    subparser = parser.add_subparsers()
        
    # subparser to call STAR aligner
    parser_star = subparser.add_parser('starAlignment',add_help=False)
    parser_star.add_argument('-r','--reference', type=argparse.FileType('r'), required=True, dest='reference')
    parser_star.add_argument('-1','--r1', type=argparse.FileType('r'), required=True, dest='r1')
    parser_star.add_argument('-2','--r2', type=argparse.FileType('r'), required=False, dest='r2')
    parser_star.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_star.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_star.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_star.add_argument('-k','--keep', action='store_false',required=False, default=False, dest='rm')
    parser_star.add_argument('-t','--threads', type=int, default=1, dest='threads')
    parser_star.add_argument('-h','--help',action='store_const', const = parser_star.prog.split()[-1], dest='c_help')
    parser_star.set_defaults(func=star)

    # subparser to call Hisat2 aligner
    parser_hisat2 = subparser.add_parser('hisat2Alignment',add_help=False)
    parser_hisat2.add_argument('-r','--reference', type=argparse.FileType('r'), required=True, dest='reference')
    parser_hisat2.add_argument('-1','--r1', type=argparse.FileType('r'), required=True,  dest='r1')
    parser_hisat2.add_argument('-2','--r2', type=argparse.FileType('r'), required=False, dest='r2')
    parser_hisat2.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_hisat2.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_hisat2.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_hisat2.add_argument('-t','--threads', type=int, default=1, dest='threads')
    parser_hisat2.add_argument('-h','--help',action='store_const', const = parser_hisat2.prog.split()[-1], dest='c_help')
    parser_hisat2.set_defaults(func=hisat2)
    
    # subparser for the split read search
    parser_split = subparser.add_parser('splitReadSearch',add_help=False)
    parser_split.add_argument('-a', '--alignment', type=argparse.FileType('r'), required=True, dest='bamfile')
    parser_split.add_argument('-r', '--reference', type=argparse.FileType('r'), required=True, dest='fastafile')
    parser_split.add_argument('-d', '--min-depth', type=int, default=1, dest='mindepth')
    parser_split.add_argument('-l', '--max-length', type=int, default=80, dest='maxlen')
    parser_split.add_argument('-f', '--foot-size', type=int, default=10, dest='minfootsize')
    parser_split.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_split.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_split.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_split.add_argument('-t','--threads', type=int, default=1, dest='threads')
    parser_split.add_argument('-h','--help',action='store_const', const = parser_split.prog.split()[-1],dest='c_help')
    parser_split.set_defaults(func=splitReadSearch)

    # subparser for writing fasta of spliced sequence
    parser_trim = subparser.add_parser('trimFastaFromTXT',add_help=False)
    parser_trim.add_argument('-r', '--reference', type=argparse.FileType('r'), required=True, dest='reference')
    parser_trim.add_argument('-c', '--candidates', type=argparse.FileType('r'),  required=True, dest='cand_file')
    parser_trim.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_trim.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_trim.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_trim.add_argument('-m', '--multi', action='store_true', default=False, dest='multi')
    parser_trim.add_argument('-h','--help',action='store_const', const = parser_trim.prog.split()[-1],dest='c_help')
    parser_trim.set_defaults(func=trimFastaFromTXT)

    # subparser for HTML simulation report
    parser_html = subparser.add_parser('buildReport',add_help=False)
    parser_html.add_argument('--config_file')
    args_html, left_html_argv = parser_html.parse_known_args()
    
    if args_html.config_file:
        with open(args_html.config_file, 'r') as f:
            config_html = configparser.ConfigParser()
            config_html.read([args_html.config_file])

    parser_html.add_argument('--config-file', type=argparse.FileType('r'), required=False, help="Provide a config file")
    parser_html.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='fasta', help="Path to the reference FASTA file.")
    parser_html.add_argument('-m','--modified-fasta', type=argparse.FileType('r'), required=False, dest='mfasta', help="Path to the modified FASTA file.")
    parser_html.add_argument('-g','--gtf', type=argparse.FileType('r'), required=False, dest='gtf', help="GTF filename which contains the genome annotation.")
    parser_html.add_argument('-1','--R1', type=argparse.FileType('r'), required=True, dest='r1', help="Name of the  FASTQ  file  which  contains  the  single-end   reads library. If paired-end, filename of #1 reads mates")
    parser_html.add_argument('-2','--R2', type=argparse.FileType('r'), required=False, dest='r2', help="Only for a paired-end library, filename of #2 reads mates.")
    parser_html.add_argument('--flagstat', type=argparse.FileType('r'), required=False, dest='flagstat', help="Path to flagstat file.")
    parser_html.add_argument('-r','--ranks', type=argparse.FileType('r'), required=False, dest='ranks', help="Path to ranks file.")
    parser_html.add_argument('-c','--candidat', type=argparse.FileType('r'), required=False, dest='candidat', help="Path to candidat file.")
    parser_html.add_argument('-s','--split', type=argparse.FileType('r'), required=False, dest='split', help="Path to split file.")
    parser_html.add_argument('-o','--output', type=str, required=True, dest='output', help="Output dir name.")
    parser_html.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix', help="Prefix for output files name.")
    parser_html.add_argument('-t','--threads', type=int, default=1, required=False, dest='threads', help="Number of threads [1]")
    parser_html.add_argument('-F', '--force', action='store_true', default=False, dest='force', help="Force to overwrite output files.")
         
    
    try:
        config_html
    except NameError:
        pass
    else:
        for k, v in config_html.items("Defaults"):
            config_html_args={str(k): str(v)}
            # Use values from configuration file by default
            parser_html.set_defaults(**config_html_args)
            # Reset `required` attribute when provided from config file
            for action in parser_html._actions:
                if action.dest in config_html_args:
                    action.required = False   
    parser_html.set_defaults(func=simulationReport) 		
    parser_html.add_argument('-h','--help',action='store_const', const = parser_html.prog.split()[-1],dest='c_help')

    # subparser for analyze candidates (findEvidence)
    parser_findEvidence = subparser.add_parser('findEvidence',add_help=False)
    parser_findEvidence.add_argument('-r', '--reference', type=argparse.FileType('r'), required=True, dest='reference')
    parser_findEvidence.add_argument('-t', '--trim-reference', type=argparse.FileType('r'), required=True, dest='trim_ref')
    parser_findEvidence.add_argument('-d', '--database', type=argparse.FileType('r'),  required=True, dest='db_prot')
    parser_findEvidence.add_argument('-c', '--candidates', type=argparse.FileType('r'),  required=True, dest='cand_file')
    parser_findEvidence.add_argument('-k','--keep', action='store_false',required=False, default=False, dest='rm')
    parser_findEvidence.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_findEvidence.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_findEvidence.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_findEvidence.add_argument('-h','--help',action='store_const', const = parser_findEvidence.prog.split()[-1],dest='c_help')
    parser_findEvidence.set_defaults(func=findEvidence)

    # subparser for Full Random Simulation
    parser_frs = subparser.add_parser('fullRandomSimulation',add_help=False)
    parser_frs.add_argument('-n','--nb_contigs', type=int, default=2000, dest='nb')
    parser_frs.add_argument('-m','--min-contig-len', type=int, default=250, dest='mini')
    parser_frs.add_argument('-M','--max-contig-len', type=int, default=1500, dest='maxi')
    parser_frs.add_argument('-r','--random-part', type=int, default=100, dest='part')
    parser_frs.add_argument('--mix','--mix-state', action='store_true', default=False, dest='mix')
    parser_frs.add_argument('-l', '--lower-intron-len',  type=int, default=250, dest='lower')
    parser_frs.add_argument('-H', '--higher-intron-len', type=int, default=750, dest='upper')
    parser_frs.add_argument('-o', '--outputDir', type=str, required=True, dest='output')
    parser_frs.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_frs.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_frs.add_argument('-h','--help',action='store_const', const = parser_frs.prog.split()[-1],dest='c_help')
    parser_frs.set_defaults(func=full_random_simulation)

    # subparser for annotated genome-based data simulation (annoToData) 
    parser_gbs = subparser.add_parser('GTFbasedSimulation',add_help=False)
    parser_gbs.add_argument('-a','--annotation', type=argparse.FileType('r'), required=True, dest='annotation')
    parser_gbs.add_argument('-r', '--reference', type=argparse.FileType('r'), required=True, dest='fasta')
    group_nb = parser_gbs.add_mutually_exclusive_group() 
    group_nb.add_argument('-n','--nb-transcripts', type=int, required=False, default=0, dest='nb')
    group_nb.add_argument('-N','--all-transcripts', action='store_const', const=0, default=False, dest='nb')
    parser_gbs.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_gbs.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_gbs.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_gbs.add_argument('-m', '--mix-state', action='store_true', default=False, dest='mix')
    parser_gbs.add_argument('-u', '--uniq-transcript', action='store_true', default=False, dest='uniq')
    parser_gbs.add_argument('-h','--help',action='store_const', const = parser_gbs.prog.split()[-1],dest='c_help')
    parser_gbs.set_defaults(func=gtf_based_simulation)

    # subparser for grinder
    parser_grinder = subparser.add_parser('simulateReads',add_help=False, help='Grinder help. Needs software Grinder-v0.5.4 or more recent version')
    parser_grinder.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='rf')
    parser_grinder.add_argument('-c','--cfg', type=argparse.FileType('r'), required=True, dest='pf')
    parser_grinder.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_grinder.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_grinder.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_grinder.add_argument('-h','--help',action='store_const', const = parser_grinder.prog.split()[-1],dest='c_help')
    parser_grinder.set_defaults(func=grinder)

    # subparser for checkInstall
    parser_install = subparser.add_parser('checkInstall',add_help=False, help='Check module to verify if the running environment and the dependencies are correctly installed')
    parser_install.add_argument('-h','--help',action='store_const', const = parser_install.prog.split()[-1],dest='c_help')
    parser_install.set_defaults(func=checkInstall) ;
    
    ### Parsing of command line
    
    # Printing of help message without '-h' option
    if len(sys.argv) == 1 : 
        program_help()
        exit()
    elif len(sys.argv) > 1 :
        h_command = sys.argv[1]
    
    # Trying to parse the arguments.
    # If errors are raised, printing of error message
    try :
        args = vars(parser.parse_args())
    except :
        command_help(h_command)
        exit(2)

    #if args.pop('checkInstall') :
    #    print_to_stdout("### checkInstallcheckInstallcheckInstallcheckInstallcheckInstall ###")


    # Printing of help or version messages with the 'help' or 'version' option.
    if args.pop('help') :
        program_help()
        exit()
    elif args.pop('version') :
        program_version()
        exit()
    else :
        h_command = args.pop('c_help')
        if h_command :
            command_help(h_command)
            exit()
        else :
            return args


if __name__ == '__main__':
    args = parse_arguments()
    
    # Run the program
    func = args.pop('func')
    func(**args)
    

