#!/usr/bin/env python3

# Modules : 
import argparse 
import sys 
import os 
from intronSearch import searchProtein,predictORF,truncate,split_research
from readsMapping import star,hisat2
from dataSimulation import full_random_simulation,gtf_based_simulation,grinder
from checkInstall import checkInstall
from helpMessages import program_help,command_help,program_version

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
    parser_split.add_argument('-a', '--alignment', type=argparse.FileType('r'), required=True, dest='bamfilename')
    parser_split.add_argument('-r', '--reference', type=argparse.FileType('r'), required=True, dest='fastafilename')
    parser_split.add_argument('-o','--output', type=str, required=True, dest='output')
    parser_split.add_argument('-F', '--force', action='store_true', default=False, dest='force')
    parser_split.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_split.add_argument('-h','--help',action='store_const', const = parser_split.prog.split()[-1],dest='c_help')
    parser_split.set_defaults(func=split_research)

    # subparser for writing fasta of spliced sequence
    parser_trim = subparser.add_parser('trimFastaFromTXT',add_help=False)
    parser_trim.add_argument('-r', '--reference', type=argparse.FileType('r'), dest='fasta_file')
    parser_trim.add_argument('-f', '--features', type=argparse.FileType('r'), dest='gff_feature')
    parser_trim.add_argument('-o','--output', type=str, default='sequences', dest='output')
    parser_trim.add_argument('-h','--help',action='store_const', const = parser_trim.prog.split()[-1],dest='c_help')
    parser_trim.set_defaults(func=truncate)

    # subparser for searching ORF on sequences
    parser_orf = subparser.add_parser('analyzeORF',add_help=False)
    parser_orf.add_argument('-r', '--reference', type=argparse.FileType('r'), dest='fasta_file')
    parser_orf.add_argument('-k','--keep', action='store_false',required=False, default=False, dest='rm')
    parser_orf.add_argument('-o','--output', type=str, dest='output', default='predicted_orfs')
    parser_orf.add_argument('--no-refine-starts', action='store_false', required=False, default=True, dest='refine')
    parser_orf.add_argument('-h','--help',action='store_const', const = parser_orf.prog.split()[-1],dest='c_help')
    parser_orf.set_defaults(func=predictORF)

    # subparser for aligning contigs and proteins
    parser_protein = subparser.add_parser('analyzeProtein',add_help=False)
    parser_protein.add_argument('-r', '--reference', type=argparse.FileType('r'), dest='fasta')
    parser_protein.add_argument('-p','--dbprotein', type=str, dest='dbprotein')
    parser_protein.add_argument('-o', '--output', type=str, required=False, default='LongGaps.gff', dest='output')
    parser_protein.add_argument('-k','--keep', action='store_true',required=False, default=False, dest='rm')
    parser_protein.add_argument('-t','--threads', type=int, default=1, dest='threads')
    parser_protein.add_argument('-h','--help', action='store_const', const = parser_protein.prog.split()[-1],dest='c_help')
    parser_protein.set_defaults(func=searchProtein)

    # subparser for Full Random Simulation
    parser_frs = subparser.add_parser('fullRandomSimulation',add_help=False)
    parser_frs.add_argument('-n','--nb_contigs', type=int, default=2000, dest='nb')
    parser_frs.add_argument('-m','--min-contig-len', type=int, default=250, dest='mini')
    parser_frs.add_argument('-M','--max-contig-len', type=int, default=1500, dest='maxi')
    parser_frs.add_argument('-r','--random-half', default=False, action='store_true',dest='half')
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
    parser_gbs.add_argument('-h','--help',action='store_const', const = parser_gbs.prog.split()[-1],dest='c_help')
    parser_gbs.set_defaults(func=gtf_based_simulation)

    # subparser for grinder
    parser_grinder = subparser.add_parser('simulateReads',add_help=False, help='Grinder help. Needs software Grinder-v0.5.4 or more recent version')
    parser_grinder.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='rf')
    parser_grinder.add_argument('-c','--cfg', type=argparse.FileType('r'), required=True, dest='pf')
    parser_grinder.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser_grinder.add_argument('-o','--outputDir', type=str, required=True, dest='output')
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

