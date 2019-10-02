#!/usr/bin/env python3

# Modules : 
import argparse 
import sys 
import os 
from intronResearch import searchProtein,predictORF,truncate,split_research
from readsMapping import star,hisat2
from dataSimulation import full_random_simulation,gtf_based_simulation,grinder
from checkInstall import checkInstall
from helpMessages import program_help,command_help




def parse_arguments() :
    
    if len(sys.argv) == 1 : 
        program_help()
        exit()
    elif len(sys.argv) ==  2 :
        pass
        
    parser = argparse.ArgumentParser(
        add_help=False
        )
    
    subparser = parser.add_subparsers()
    
    # subparser to call STAR aligner
    parser_star = subparser.add_parser('starAlignement',add_help=False, help='STAR help. Needs aligner STAR-2.6.0c or more recent version')
    parser_star.add_argument('-i','--reference',
                             help='reference file', type=str, required=True, dest='reference')
    parser_star.add_argument('-1','--r1',
                             help='fasta file of reads', type=str, required=True,dest='r1')
    parser_star.add_argument('-2','--r2',
                             help='second fasta file if reads if paired-end', type=str, default='',dest='r2')
    parser_star.add_argument('-o','--prefix',
                             help='output prefix', type=str, default='STAR',dest='prefix')
    parser_star.add_argument('-t','--threads',
                             help='number of threads used to perform alignement', type=int, default=1,dest='threads')
    parser_star.set_defaults(func=star)

    # subparser to call Hisat2 aligner
    parser_hisat2 = subparser.add_parser('hisat2Alignement',add_help=False, help='HiSat2 help. Needs aligner HiSat2-2.1.0 or more recent version')
    parser_hisat2.add_argument('-i','--reference',
                               help='reference file', type=str, required=True, dest='reference')
    parser_hisat2.add_argument('-1','--r1',
                               help='fasta file of reads', type=str, required=True, dest='r1')
    parser_hisat2.add_argument('-2','--r2',
                               help='second fasta file if reads if paired-end', type=str, required=True, default='', dest='r2')
    parser_hisat2.add_argument('-o','--prefix',
                               help='output prefix', type=str, default='HiSat2', dest='prefix')
    parser_hisat2.add_argument('-t','--threads',
                               help='number of threads used to perform alignement', type=int, default=1, dest='threads')
    parser_hisat2.set_defaults(func=hisat2)
    
    # subparser for the split read research
    parser_split = subparser.add_parser('splitReadSearch',add_help=False, help='In the alignment file, look for split read')
    parser_split.add_argument('-i', '--input',
                              help='input bam file from alignment', dest='bamfilename', type=str, required=True)
    parser_split.add_argument('-r', '--reference',
                              help='name of the reference file used for the alignment', type=str, required=True, dest='fastafilename')
    parser_split.add_argument('-o','--prefix',
                              help='prefix of ouput files',type=str, required=True, dest='basename')
    parser_split.set_defaults(func=split_research)

    # subparser for writing fasta of spliced sequence
    parser_trunc = subparser.add_parser('trimFastaFromGFF',add_help=False,help='Write fasta file where each sequence is spliced if it bears a feature (present in the gff file provided). If a sequence bears several features, a different spliced sequence is written for each feature')
    parser_trunc.add_argument('-i', '--fasta',
                              help='Fasta with sequences to splice', dest='fasta_file', type=str, required=True)
    parser_trunc.add_argument('-g', '--gff', 
                              help='File with features to delete (GFF format only)', dest='gff_feature', type=str, required=True)
    parser_trunc.add_argument('-o','--output',
                              help='Basename of the output file', dest='output', type=str, default='sequences')
    parser_trunc.set_defaults(func=truncate)

    # subparser for searching ORF on sequences
    parser_orf = subparser.add_parser('analyzeORF',add_help=False, help='Run TransDecoder to search Open Reading Frame on sequences provided in fasta format.'
                                                     'Returns the predicted ORFs in BED, nucleic FASTA, proteic FASTA and TSV format.'
                                                     'Needs software TransDecoder-v5.5.0')
    parser_orf.add_argument('-i', '--fasta', 
                            help='Sequences on which the prediction will be performed', dest='fasta_file', type=str, required=True)
    parser_orf.add_argument('-t','--truncate-search',
                            help='Boolean which says if, after the ORF prediction is performed on fasata files, a truncation (truncate module) is performed and another ORF prediction on truncated fasta',
                            action='store_true',dest='trunc')
    parser_orf.add_argument('-g', '--gff', 
                            help='File with features to delete (GFF format only). Only used if -t option is present',dest='gff', type=str, required=False, default=None)
    parser_orf.add_argument('-k','--keep_intermdiate',
                            help='Boolean which rules intermediate files erasure (default False)', dest='rm', action='store_false',required=False, default=True)
    parser_orf.add_argument('-o','--output',
                            help='Basename of the output file', dest='output', type=str, default='predicted_orfs')
    parser_orf.add_argument('--no_refine_starts', 
                            help='Boolean which rules the refine start step of TransDecoder. Use it only if TransDecoder fails at this point.',dest='refine',action='store_false', required=False, default=True)
    parser_orf.set_defaults(func=predictORF)

    # subparser for aligning contigs and proteins to search long gaps (i.e. potential introns)
    parser_protein = subparser.add_parser('analyzeProtein',add_help=False, help='Run Diamond to search long gaps (i.e. potential introns) in contigs-protein alignement.'
                                                     'Returns a GFF file with all the found gaps'
                                                     'Needs Diamond-v0.9.9')
    parser_protein.add_argument('-i', '--fasta',
                                help='Contigs sequences', dest='fasta', type=str, required=True)
    parser_protein.add_argument('-p','--dbprotein',
                                help='Name of the Diamond database containing the indexed proteic sequences.',type=str, required=True, dest='dbprotein')
    parser_protein.add_argument('-o', '--output',
                                help='Output filename',required=False,default='LongGaps.gff', dest='output', type=str)
    parser_protein.add_argument('-k','--keep_intermediate',
                                help='Boolean which rules intermediate files erasure (default False)', dest='rm', action='store_true',required=False)
    parser_protein.add_argument('-t','--threads',
                                help='Integer which indicates the numberof CPU to use for alignement (default 1)', dest='threads', type=int, required=False,default = 1)
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
    parser_frs.add_argument('-h', '--higher-intron-length',
                            help='Maximal length of random intron. Default 1000', type=int, default=1000, dest='upper')
    parser_frs.add_argument('-o', '--output',
                            help='Basename of the generated files. Default [FullRandomSimulation] ', type=str, default='FullRandomSimulation', dest='output')
    parser_frs.set_defaults(func=full_random_simulation)

    # subparser for annotated genome-based data simulation (annoToData) 
    parser_gbs = subparser.add_parser('GTFbasedSimulation',add_help=False, help='Data simulation (genome assembly and reads library) from a genome annotation (GFF or GTF file)')
    parser_gbs.add_argument('-i','--annotation',
                            help='Filename of GFF file which contains the genome annotation', type=str, required=True, dest='annotation')
    parser_gbs.add_argument('-f','--fasta',
                            type=str, required=True, dest='fasta')
    group_nb = parser_gbs.add_mutually_exclusive_group(required=True) 
    group_nb.add_argument('-n','--nb_genes',
                          help='Total number of genes to transcript ', type=float, required=False, default=0, dest='nb')
    group_nb.add_argument('-a','--all',
                          help='Flag which says if all genes from GFF have to be transcripted', action='store_const', const=0, default=False, dest='nb')
    parser_gbs.add_argument('-o','--output',
                            help='prefix of ouput files', type=str, required=True, dest='output')
    parser_gbs.add_argument('--mix-library',
                            help='Boolean which rules if the generated library is mixed i.e. if the library contains the transcript in two state when a intron is retained or an exon is spliced', action='store_true', default=False, dest='mix')
    parser_gbs.set_defaults(func=gtf_based_simulation)

    # subparser for grinder
    parser_grinder = subparser.add_parser('simulateReads',add_help=False, help='Grinder help. Needs software Grinder-v0.5.4 or more recent version')
    parser_grinder.add_argument('-i','--rf',
                                help='reference file', type=str, required=True, dest='rf')
    parser_grinder.add_argument('-p','--pf',
                             help='Profile file : all arguments for grinder', type=str, required=False, default=os.path.abspath(os.path.dirname(sys.argv[0]) + '/../config/profile_file.txt'),dest='pf')
    parser_grinder.add_argument('-o','--pref',
                             help='Prefix of the output files', type=str, default='Grinder', dest='pref')
    parser_grinder.set_defaults(func=grinder)

    # subparser for checkInstall
    parser_install = subparser.add_parser('checkInstall',add_help=False, help='Check module to verify if the running environment and the dependencies are correctly installed')
    parser_install.set_defaults(func=checkInstall) ;
    
    
    # ~ return vars(parser.parse_args())




if __name__ == '__main__':
    
    args = parse_arguments()

    # ~ if len(args) < 1 :
        # ~ print('ici')
    # ~ elif len(args) == 1 :
        # ~ print(args)
    # ~ else:
        # ~ func = args.pop('func')
        # ~ func(**args)
