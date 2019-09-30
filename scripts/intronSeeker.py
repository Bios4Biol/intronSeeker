#!/usr/bin/env python3

# Modules : 
import argparse ; 
import sys ; 
import os ;
import intronResearch as ir ;
import readsMapping as rmg ;
import dataSimulation as ds ;
import checkInstall as ci ;


if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		description="Functions to find and analyse split in alignment file.")
	subparser = parser.add_subparsers(help="sub-command help")
	
	# subparser for the split read research
	parser_split = subparser.add_parser("split", help="In the alignment file, look for split read")
	parser_split.add_argument("-i", "--input", help="input bam file from alignment", dest="bamfilename", type=str,
							   metavar='', required=True)
	parser_split.add_argument("-r", "--reference", help="name of the reference file used for the alignment", type=str,
							   required=True, metavar='', dest="fastafilename")
	parser_split.add_argument("-o","--prefix",help="prefix of ouput files",type=str, required=True, metavar='',
							   dest="basename")
	parser_split.set_defaults(func=ir.split_research)

	# subparser for writing fasta of spliced sequence
	parser_trunc = subparser.add_parser("truncate", help="Write fasta file where each sequence is spliced if it bears a feature (present in the gff file provided). "
														  "If a sequence bears several features, a different spliced sequence is written for each feature")
	parser_trunc.add_argument("-i", "--fasta", help="Fasta with sequences to splice", dest="fasta_file", type=str,
							   metavar='', required=True)
	parser_trunc.add_argument("-g", "--gff", help="File with features to delete (GFF format only)", dest="gff_feature", type=str,
							   metavar='', required=True)
	parser_trunc.add_argument("-o","--output", help="Basename of the output file", dest="output",
							   type=str, default="sequences", metavar='')
	parser_trunc.set_defaults(func=ir.truncate)

	# subparser for searching ORF on sequences
	parser_orf = subparser.add_parser("orf", help="Run TransDecoder to search Open Reading Frame on sequences provided in fasta format."
													 "Returns the predicted ORFs in BED, nucleic FASTA, proteic FASTA and TSV format."
													 "Needs software TransDecoder-v5.5.0")
	parser_orf.add_argument("-i", "--fasta", help="Sequences on which the prediction will be performed", dest="fasta_file", type=str,
							   metavar='', required=True)
	parser_orf.add_argument("-t","--truncate-search",help="Boolean which says if, after the ORF prediction is performed on fasata files," 
							   "a truncation (truncate module) is performed and another ORF prediction on truncated fasta",action="store_true",dest="trunc")
	parser_orf.add_argument("-g", "--gff", help="File with features to delete (GFF format only). Only used if -t option is present",
							   dest="gff", type=str, metavar='', required=False,default=None)
	parser_orf.add_argument("-k","--keep_intermdiate", help="Boolean which rules intermediate files erasure (default False)", dest="rm",
							   action="store_false",required=False, default=True)
	parser_orf.add_argument("-o","--output", help="Basename of the output file", dest="output",
							   type=str, default="predicted_orfs", metavar='')
	parser_orf.add_argument("--no_refine_starts", help="Boolean which rules the refine start step of TransDecoder. Use it only if TransDecoder fails at this point.",
							   dest="refine",action="store_false", required=False, default=True)
	parser_orf.set_defaults(func=ir.predictORF)

	# subparser for aligning contigs and proteins to search long gaps (i.e. potential introns)
	parser_protein = subparser.add_parser("protein", help="Run Diamond to search long gaps (i.e. potential introns) in contigs-protein alignement."
													 "Returns a GFF file with all the found gaps"
													 "Needs Diamond-v0.9.9")
	parser_protein.add_argument("-i", "--fasta", help="Contigs sequences", dest="fasta", type=str,
							   metavar='', required=True)
	parser_protein.add_argument("-p","--dbprotein",help="Name of the Diamond database containing the indexed proteic sequences.",type=str, 
							   metavar='',required=True, dest="dbprotein")
	parser_protein.add_argument("-o", "--output", help="Output filename",required=False,default="LongGaps.gff",
							   dest="output", type=str, metavar='')
	parser_protein.add_argument("-k","--keep_intermediate", help="Boolean which rules intermediate files erasure (default False)", dest="rm",
							   action="store_true",required=False)
	parser_protein.add_argument("-t","--threads", help="Integer which indicates the numberof CPU to use for alignement (default 1)", dest="threads",
							   type=int, required=False,default = 1)
	parser_protein.set_defaults(func=ir.searchProtein)

	# subparser for the contig creation
	parser_contig = subparser.add_parser("contig", help="contig help")
	parser_contig.add_argument("-n", help="Number of contig/sequence randomly generated. Default 10", type=int,
							   default=10, metavar='',dest="n")
	parser_contig.add_argument("--min", help="lower limits of the sequence(s). Default 150", type=int, default=150,
							   metavar='',dest="minimum")
	parser_contig.add_argument("--max", help="upper limits of the sequence(s). Default 1000", type=int, default=1000,
							   metavar='',dest="maximum")
	parser_contig.add_argument("-o", "--output", help="output name/directory, default=OutputRandomContig.fa", type=str,
							   default="OutputRandomContig.fa", metavar='',dest="output")
	parser_contig.set_defaults(func=ds.contig)

	# subparser for the input insertion
	parser_intron = subparser.add_parser("intron", help="intron help")
	parser_intron.add_argument("-i", "--input", help="input File", type=str, required=True, metavar='',dest="input_file")
	parser_intron.add_argument("-l", "--lower", help="lower limits of intron length.", type=int, default=150,
							   metavar='',dest="lower")
	parser_intron.add_argument("-u", "--upper", help="upper limits of intron length.", type=int, default=1000,
							   metavar='',dest="upper")
	parser_intron.add_argument("--bi", help="5 intron limit. Default GT", type=str, default="GT", metavar='',dest="bi")
	parser_intron.add_argument("--bs", help="3 intron limit. Default AG", type=str, default="AG", metavar='',dest="bs")
	parser_intron.add_argument("-o", "--output", help="output name/directory, default=OutputIntron.fa", type=str,
							   default="OutputIntron.fa", metavar='',dest="output")
	parser_intron.add_argument("-c", "--coorf", help="output name/directory of the intron information file", type=str,
							   default="IntronCoord.txt", metavar='',dest="coorf")
	parser_intron.add_argument("-b", "--begin", help="file with the beginning of intron for all sequences and if they "
													 "are reverse", type=str,
							   default="", metavar='',dest="begin")
	parser_intron.add_argument("-r", "--rand", help="insert intron in randomly half "
													"of the sequences if the option is specified", action="store_true",dest="rand")
	parser_intron.set_defaults(func=ds.write_intron)

	# subparser for grinder
	parser_grinder = subparser.add_parser("grinder", help="Grinder help. Needs software Grinder-v0.5.4 or more recent version")
	parser_grinder.add_argument("-i","--rf", help="reference file", type=str, required=True, metavar='',dest="rf")
	parser_grinder.add_argument("-p","--pf",
							 help="Profile file : all arguments for grinder", type=str, required=False, 
							 default=os.path.abspath(os.path.dirname(sys.argv[0]) + "/../config/profile_file.txt"),
							 metavar='',dest="pf")
	parser_grinder.add_argument("-o","--pref",
							 help="Prefix of the output files", type=str,
							 default="Grinder", metavar='',dest="pref")
	parser_grinder.set_defaults(func=ds.grinder)

	# subparser for STAR
	parser_star = subparser.add_parser("star", help="STAR help. Needs aligner STAR-2.6.0c or more recent version")
	parser_star.add_argument("-i","--reference",
							 help="reference file", type=str, required=True, metavar='',dest="reference")
	parser_star.add_argument("-1","--r1",
							 help="fasta file of reads", type=str, required=True, metavar='',dest="r1")
	parser_star.add_argument("-2","--r2",
							 help="second fasta file if reads if paired-end", type=str, default="", metavar='',dest="r2")
	parser_star.add_argument("-o","--prefix",
							 help="output prefix", type=str, default="STAR", metavar='',dest="prefix")
	parser_star.add_argument("-t","--threads",
							 help="number of threads used to perform alignement", type=int, default=1, metavar='',dest="threads")
	parser_star.set_defaults(func=rmg.star)

	# subparser for HiSat2
	parser_hisat2 = subparser.add_parser("hisat2", help="HiSat2 help. Needs aligner HiSat2-2.1.0 or more recent version")
	parser_hisat2.add_argument("-i","--reference",
							 help="reference file", type=str, required=True, metavar='',dest="reference")
	parser_hisat2.add_argument("-1","--r1",
							 help="fasta file of reads", type=str, required=True, metavar='',dest="r1")
	parser_hisat2.add_argument("-2","--r2",
							 help="second fasta file if reads if paired-end", type=str, required=True, default="", metavar='',dest="r2")
	parser_hisat2.add_argument("-o","--prefix",
							 help="output prefix", type=str, default="HiSat2", metavar='',dest="prefix")
	parser_hisat2.add_argument("-t","--threads",
							 help="number of threads used to perform alignement", type=int, default=1, metavar='',dest="threads")
	parser_hisat2.set_defaults(func=rmg.hisat2) ;

	# subparser for checkInstall
	parser_install = subparser.add_parser("checkInstall", help="Check module to verify if the running environment and the dependencies are correctly installed")
	parser_install.set_defaults(func=ci.checkInstall) ;

	# subparser for Annotation-based data simulation (annoToData) 
	parser_annotodata = subparser.add_parser("annoToData", help="Data simulation (genome assembly and reads library) from a genome annotation (GFF or GTF file)")
	parser_annotodata.add_argument("-i","--annotation", help="Filename of GFF file which contains the genome annotation",type=str,metavar='',required=True,dest="annotation") ;
	parser_annotodata.add_argument("-f","--fasta",type=str,metavar='',required=True,dest="fasta")
	group_nb = parser_annotodata.add_mutually_exclusive_group(required=True) 
	group_nb.add_argument("-n","--nb_genes", help="Total number of genes to transcript ",type=float, required=False, default = 0,dest="nb",metavar='') ;
	group_nb.add_argument("-a", "--all", help="Flag which says if all genes from GFF have to be transcripted", action="store_const",const=0,default = False, dest="nb") ;
	parser_annotodata.add_argument("-o","--output", help="prefix of ouput files", type=str, required=True,dest="output",metavar='') ; 
	parser_annotodata.add_argument("--no-grinder",help="Boolean which rules if let the library at temporary state (i.e. before the reads library generation)",
								 action="store_false",default = True, dest="grinder") ;
	parser_annotodata.add_argument("--mix-library",help="Boolean which rules if the generated library is mixed i.e. if the library contains the transcript in two state when a intron is retained or an exon is spliced",
								 action="store_true",default = False, dest="mix") ;
	parser_annotodata.set_defaults(func=ds.annoToData)

	args = parser.parse_args()


	if len(vars(args)) < 1 :
		parser.print_help()
	else:
		# need to change it
		param = {a: b for a, b in vars(args).items() if a != "func"}
		args.func(**param)
