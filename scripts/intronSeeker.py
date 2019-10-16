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
    parser = argparse.ArgumentParser(
        add_help=False
        )
    parser.add_argument('-v','--version',action='store_true',required=False,dest='version')
    parser.add_argument('-h','--help',action='store_true',required=False,dest='help')
    
    subparser = parser.add_subparsers()
    
    # subparser to call STAR aligner
    parser_star = subparser.add_parser('starAlignement',add_help=False)
    parser_star.add_argument('-r','--reference', type=str, dest='reference')
    parser_star.add_argument('-1','--r1', type=str, dest='r1')
    parser_star.add_argument('-2','--r2', type=str, default='',dest='r2')
    parser_star.add_argument('-o','--output', type=str, dest='prefix')
    parser_star.add_argument('-t','--threads', type=int, default=1,dest='threads')
    parser_star.add_argument('-h','--help',action='store_const', const = parser_star.prog.split()[-1],dest='c_help')
    parser_star.set_defaults(func=star)

    # subparser to call Hisat2 aligner
    parser_hisat2 = subparser.add_parser('hisat2Alignement',add_help=False)
    parser_hisat2.add_argument('-r','--reference', type=str, dest='reference')
    parser_hisat2.add_argument('-1','--r1', type=str, dest='r1')
    parser_hisat2.add_argument('-2','--r2', type=str, default='', dest='r2')
    parser_hisat2.add_argument('-o','--output', type=str, default='HiSat2', dest='prefix')
    parser_hisat2.add_argument('-t','--threads', type=int, default=1, dest='threads')
    parser_hisat2.add_argument('-h','--help',action='store_const', const = parser_hisat2.prog.split()[-1],dest='c_help')
    parser_hisat2.set_defaults(func=hisat2)
    
    # subparser for the split read research
    parser_split = subparser.add_parser('splitReadSearch',add_help=False)
    parser_split.add_argument('-a', '--alignement', dest='bamfilename', type=str, )
    parser_split.add_argument('-r', '--reference', type=str, dest='fastafilename')
    parser_split.add_argument('-o','--output', type=str, dest='basename')
    parser_split.add_argument('-h','--help',action='store_const', const = parser_split.prog.split()[-1],dest='c_help')
    parser_split.set_defaults(func=split_research)

    # subparser for writing fasta of spliced sequence
    parser_trim = subparser.add_parser('trimFastaFromTXT',add_help=False)
    parser_trim.add_argument('-r', '--reference', dest='fasta_file', type=str, )
    parser_trim.add_argument('-f', '--features', dest='gff_feature', type=str, )
    parser_trim.add_argument('-o','--output', dest='output', type=str, default='sequences')
    parser_trim.add_argument('-h','--help',action='store_const', const = parser_trim.prog.split()[-1],dest='c_help')
    parser_trim.set_defaults(func=truncate)

    # subparser for searching ORF on sequences
    parser_orf = subparser.add_parser('analyzeORF',add_help=False)
    parser_orf.add_argument('-r', '--reference', dest='fasta_file', type=str, )
    parser_orf.add_argument('-k','--keep-intermediate', dest='rm', action='store_false',required=False, default=True)
    parser_orf.add_argument('-o','--output', dest='output', type=str, default='predicted_orfs')
    parser_orf.add_argument('--no-refine-starts', dest='refine',action='store_false', required=False, default=True)
    parser_orf.add_argument('-h','--help',action='store_const', const = parser_orf.prog.split()[-1],dest='c_help')
    parser_orf.set_defaults(func=predictORF)

    # subparser for aligning contigs and proteins
    parser_protein = subparser.add_parser('analyzeProtein',add_help=False, help='Run Diamond to search long gaps (i.e. potential introns) in contigs-protein alignement.'
                                                     'Returns a GFF file with all the found gaps'
                                                     'Needs Diamond-v0.9.9')
    parser_protein.add_argument('-i', '--fasta',
                                help='Contigs sequences', dest='fasta', type=str, )
    parser_protein.add_argument('-p','--dbprotein',
                                help='Name of the Diamond database containing the indexed proteic sequences.',type=str, dest='dbprotein')
    parser_protein.add_argument('-o', '--output',
                                help='Output filename',required=False,default='LongGaps.gff', dest='output', type=str)
    parser_protein.add_argument('-k','--keep_intermediate',
                                help='Boolean which rules intermediate files erasure (default False)', dest='rm', action='store_true',required=False)
    parser_protein.add_argument('-t','--threads',
                                help='Integer which indicates the numberof CPU to use for alignement (default 1)', dest='threads', type=int, required=False,default = 1)
    parser_protein.add_argument('-h','--help',action='store_const', const = parser_protein.prog.split()[-1],dest='c_help')
    parser_protein.set_defaults(func=searchProtein)

    # subparser for Full Random Simulation
    parser_frs = subparser.add_parser('fullRandomSimulation',add_help=False, help='contig help')
    parser_frs.add_argument('-n','--nb_contigs',
                            help='Number of contig/sequence randomly generated. Default 10', type=int, default=10, dest='nb')
    parser_frs.add_argument('-m','--min-contig-length',
                            help='Minimal length of random contigs. Default 150', type=int, default=150, dest='mini')
    parser_frs.add_argument('-M','--max-contig-length', 
                            help='Maximal length of random contigs. Default 1000', type=int, default=1000, dest='maxi')
    parser_frs.add_argument('--random-half',
                            help='Insert intron in random half of the simulated contigs if the option is specified. Default False',default=False, action='store_true',dest='half')
    parser_frs.add_argument('-l', '--lower-intron-length',
                            help='Minimal length of random intron. Default 150', type=int, default=150, dest='lower')
    parser_frs.add_argument('-H', '--higher-intron-length',
                            help='Maximal length of random intron. Default 1000', type=int, default=1000, dest='upper')
    parser_frs.add_argument('-o', '--output',
                            help='Basename of the generated files. Default [FullRandomSimulation] ', type=str, default='FullRandomSimulation', dest='output')
    parser_frs.add_argument('-h','--help',action='store_const', const = parser_frs.prog.split()[-1],dest='c_help')
    parser_frs.set_defaults(func=full_random_simulation)

    # subparser for annotated genome-based data simulation (annoToData) 
    parser_gbs = subparser.add_parser('GTFbasedSimulation',add_help=False, help='Data simulation (genome assembly and reads library) from a genome annotation (GFF or GTF file)')
    parser_gbs.add_argument('-i','--annotation',
                            help='Filename of GFF file which contains the genome annotation', type=str, dest='annotation')
    parser_gbs.add_argument('-f','--fasta',
                            type=str, dest='fasta')
    group_nb = parser_gbs.add_mutually_exclusive_group() 
    group_nb.add_argument('-n','--nb_genes',
                          help='Total number of genes to transcript ', type=float, required=False, default=0, dest='nb')
    group_nb.add_argument('-a','--all',
                          help='Flag which says if all genes from GFF have to be transcripted', action='store_const', const=0, default=False, dest='nb')
    parser_gbs.add_argument('-o','--output',
                            help='prefix of ouput files', type=str, dest='output')
    parser_gbs.add_argument('--mix-library',
                            help='Boolean which rules if the generated library is mixed i.e. if the library contains the transcript in two state when a intron is retained or an exon is spliced', action='store_true', default=False, dest='mix')
    parser_gbs.add_argument('-h','--help',action='store_const', const = parser_gbs.prog.split()[-1],dest='c_help')
    parser_gbs.set_defaults(func=gtf_based_simulation)

    # subparser for grinder
    parser_grinder = subparser.add_parser('simulateReads',add_help=False, help='Grinder help. Needs software Grinder-v0.5.4 or more recent version')
    parser_grinder.add_argument('-i','--rf',
                                help='reference file', type=str, dest='rf')
    parser_grinder.add_argument('-p','--pf',
                             help='Profile file : all arguments for grinder', type=str, required=False, default=os.path.abspath(os.path.dirname(sys.argv[0]) + '/../config/profile_file.txt'),dest='pf')
    parser_grinder.add_argument('-o','--pref',
                             help='Prefix of the output files', type=str, default='Grinder', dest='pref')
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
    
    # Trying to parse the arguments.
    # If errors are raised, printing of error message
    try :
        args = vars(parser.parse_args())
    except :
        print("\n***",file=sys.stderr)
        print("To know how to call intronSeeker program, use 'intronSeeker --help'.",file=sys.stderr)
        print("***",file=sys.stderr,end="\n\n")
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
        if len(sys.argv) == 2 :
           h_command = sys.argv[1]
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

