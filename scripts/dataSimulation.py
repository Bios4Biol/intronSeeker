#!/usr/bin/env python3

# Modules
try :
	import gzip ;
	import os ;
	import sys ;
	import subprocess ;
	import random ;
	import numpy as np ;
	import numpy.random as rd ;
	import configparser ;
	import re ;
	import pandas as pd ;
	from Bio.Seq import Seq ;
	from Bio import SeqIO ;
	from Bio.SeqRecord import SeqRecord ;
	from collections import defaultdict ;
except ImportError as error :
	print(error) ;
	exit(1) ;

from pprint import pprint

def _grinder(input_file : str, profile_file : str, output_file : str) :

	os.system("grinder -rf {input_file} -pf {profile_file} -bn {output_file}".format(
		input_file=input_file, profile_file=profile_file, output_file=output_file
		)) ;

def random_seq(length: int, letter: str) :
	"""
		Create a sequence based on letter and as long as length

		:param length: length of the sequence
		:type length: int
		:param letter: letter used in the sequence
		:type letter: str
		:return: sequence
		:rtype: str
	"""
	sequence = "" ;
	for l in range(0, length) :
		sequence += random.choice(letter) ; # Random choice of each letter 
	return sequence ;


def contig(n: int, minimum: int, maximum: int, output: str) :
	"""
		Write n sequence of [min, max] length in the output file

		:param n: number of sequence to create
		:type n: int
		:param minimum: minimal length of sequences
		:type minimum: int
		:param maximum: maximal length of sequences
		:type maximum: int
		:param output: name of the output fasta file
		:type output: str
		:return: nothing
	"""
	outdir = output.split(".")[0] + "_contig" ; # Output directory name (where this function will write 
												# all the results and tmp files).
	os.system("mkdir " + outdir) ;
	output_path = outdir + "/" + output ; # Path to output file

	if minimum > maximum :
		print('Maximal limit is inferior to minimal limit. Please correct it and try again') ;
	else:
		print("Writing of Fasta file...") ;
		with open(output_path, 'w') as out :
			for p in range(0, n):
				length = random.randint(minimum, maximum) ; # Random contig length 
				name = "SEQUENCE" + str(p+1) ; # Contig name creation
				seq = SeqRecord(Seq(random_seq(length, "ATCG")), id=name, description="") ; # Random sequence creation
				SeqIO.write(seq, out, "fasta") ; # Write Sequence in fasta format


def intron_begin_position(input_file: str, output: str):
	"""
	
	:param input_file:
	:param output:
	:return:
	"""
	with open(output, 'w') as out, open(input_file, 'rU') as inp:
		for record in SeqIO.parse(inp, "fasta"):
			name = record.id
			index = random.randint(1, len(record.seq))
			rev = random.choice(['nr', 're'])
			out.write("\t".join(map(str, [name, index, rev]))+"\n")


def intron(seq, name: str, coor: list, minimum: int, maximum: int, rand: bool, bi: str, bs: str, begin: dict):
	"""
	Create the intron, place it in the sequence
	and write its position in the sequence name and in a the coor list
	:param begin:
	:param seq: sequence to add an intron
	:param name: name of the sequence
	:param coor: list of all previous intron position
	:param minimum: minimal length of intron
	:param maximum: maximal length of intron
	:param rand: if TRUE, must insert intron only in randomly half of the sequence
	:param bi: 5' limit of the intron
	:param bs: 3' limit of the intron
	:return: seqn = the new sequence, name = the new name, coor = the update coor list
	:rtype : seqn = Seq, name = str, coor = list
	"""
	name_seq = name
	seqn = seq
	insert = 1  # 1 : insert intron/ 0 : no intron insertion
	if rand:
		insert = random.randint(0, 1)
	if insert == 1:
		seqn = str(seqn)
		rev = begin[name][1]  # nr : no reverse transcription / re : reverse transcription
		intron_length = random.randint(minimum, maximum) - 3
		intron_seq = "".join([bi, random_seq(intron_length, "ATCG"), bs])
		intron_length += len(bi) + len(bs)
		index = int(begin[name][0]) - 1
		name_seq = "|".join(map(str, [name, index + 1, index + intron_length, rev]))
		list_coor = [name_seq, index + 1, index + intron_length, rev]
		seqn = Seq("".join([seqn[:index], intron_seq, seqn[index:]]))
		if rev == 're':
			seqn = Seq.reverse_complement(seqn)
			name_seq = "|".join(map(str, [name, len(seqn) - index - intron_length + 1, len(seqn) - index, rev]))
			list_coor = [name_seq, len(seqn) - index - intron_length + 1, len(seqn) - index, rev]
		coor.append(list_coor)
	return seqn, name_seq, coor


def write_intron(input_file: str, output: str, lower: int, upper: int, rand: bool, coorf: str, bi: str, bs: str,
				 begin: str):
	"""
	Write the output fasta file and the intron-location txt file
	:param begin:
	:param input_file: name/path of the fasta file
	:param output: name of the fasta output file
	:param lower: minimal intron length
	:param upper: maximal intron length
	:param rand: if TRUE, must insert intron only in randomly half of the sequence
	:param coorf: name of the intron-position file
	:param bi: 5' intron limit
	:param bs: 3' intron limit
	:return: nothing 
	"""
	outdir = output.split(".")[0] + "_intron" ; # Output directory name (where this function will write 
												# all the results and tmp files).
	os.system("mkdir " + outdir) ;

	coor = []
	if lower > upper:
		print('Limit max is inferior to limit min. Please correct it and try again')
	else:
		if begin == "":
			begin = "Intron_begin_rev.txt"
			intron_begin_position(input_file, outdir + "/" + begin)
		begin_dict = {}
		with open(outdir + "/" + begin, "r") as file_b:
			for line in file_b:
				info = line[:-1].split("\t")
				begin_dict[info[0]] = info[1:]
		with open(outdir + "/" + output, 'w') as out, open(input_file, 'rU') as inp:
			for record in SeqIO.parse(inp, "fasta"):
				name = record.id
				sequence = record.seq
				new_seq, name, coor = intron(sequence, name, coor, lower, upper, rand, bi, bs, begin_dict)
				new_seq = SeqRecord(new_seq, id=name, description="")
				SeqIO.write(new_seq, out, "fasta")

		# save the intron position in coorf file
		with open(outdir + "/" + coorf, 'w') as coord:
			coord.write("\n".join("\t".join(map(str, value)) for value in coor))


def grinder(rf: str, pf: str, pref: str):
	"""
	Generate reads for a reference file depending on grinder parameters from the profile file.
	Call the software Grinder.
	:param rf: reference file
	:param pf:profile file containing all grinder parameters
	:param pref: prefix of the output files
	:return:
	"""
	outdir = pref + "_grinder" ;
	os.system("mkdir " + outdir ) ;
	output_path = outdir + "/" + pref ;

	_grinder(rf, pf, output_path) ;
	split_read(output_path)


def split_read(output_path : str):
	"""
	Split a paired-end fasta file in two fasta files.q
	:param input_file: name/path of the fasta file
	:param pref: prefix oh the two output file
	:return:
	"""
	outdir = os.path.dirname(output_path) ;
	pref = os.path.basename(output_path) ;
	input_file = outdir + "/" + [n for n in os.listdir(outdir) if n.startswith(pref+"-reads.")].pop() ;
	pile = 0 ;
	if input_file.endswith(".fa") :
		with open(input_file, 'rU') as filer, gzip.open(output_path + "_read_1.fa.gz", 'wt') as read1, \
				gzip.open(output_path + "_read_2.fa.gz", 'wt') as read2:
			for record in SeqIO.parse(filer, "fasta"):
				if pile % 2 == 0:
					SeqIO.write(record, read1, "fasta") ;
				else:
					SeqIO.write(record, read2, "fasta") ;
				pile += 1
	elif input_file.endswith(".fastq") :
		with open(input_file, 'rU') as filer, gzip.open(output_path + "_read_1.fastq.gz", 'wt') as read1, \
				gzip.open(output_path + "_read_2.fastq.gz", 'wt') as read2:
			for record in SeqIO.parse(filer, "fastq"):
				if pile % 2 == 0:
					SeqIO.write(record, read1, "fastq") ;
				else:
					SeqIO.write(record, read2, "fastq") ;
				pile += 1 ; 
	os.remove(input_file) ;


def intron_unique(pref: str):
	"""
	Create the gff file of the intron from star alignment
	:param pref: output prefix
	:return:
	"""
	data = pd.read_table(pref + ".txt", names=['ref', 'start', 'end', 'query', 'qstart', 'qend'])
	data = data.drop_duplicates(subset=('ref', 'start', 'end', 'qstart', 'qend'))
	data = pd.DataFrame(
		{'col0_ref': data.ref, 'col1': ".", "col2_class": "Intron", "col3_start": data.qstart, "col4_end": data.qend,
		 'col5_quality': "5000", 'col6': ".", 'col7': ".", 'col8': "."})
	data.to_csv(pref+".gff", header=False, sep="\t", index=False)



def countGenes(gff : str) :
	"""
	Count the protein-coding genes in the GFF file in order to randomly
	 pick among them for simulate the transcripts.
	
	:param gff: Name of the genome reference GFF file.
	:return: Number of genes. 
	"""
	
	grep = subprocess.Popen(["grep", "gene.*=protein_coding;",gff], stdout=subprocess.PIPE) ;
	wcline = subprocess.Popen(["wc", "-l"], stdin=grep.stdout, stdout=subprocess.PIPE) ;
	# Allow grep to receive a SIGPIPE if wcline exits before grep
	grep.stdout.close() ; 
	# Run the commands 
	sdo = wcline.communicate()[0] ;
	# Extract the result from standard output
	nb_genes = int(re.search(r'(\d+)',str(sdo)).group(1)) ;
	
	return nb_genes ;


def chooseGenes(gff : str, nb_genes = int) :
	"""
	Randomly pick 
	"""
	try :
		total_genes = countGenes(gff) ;
		if nb_genes < 0 :
			raise ValueError ;
		elif nb_genes >= total_genes or nb_genes == 0 :
			choosen = range(total_genes);
		else :
			choosen = rd.choice(range(total_genes),int(nb_genes),replace=False) ;
	except ValueError :
		print("*Value Error* : Number of genes to select can't be negative.") ;
		exit(1);

	return choosen ;


def generateTranscripts(gff_file : str, choosen_genes : list, output_path : str, mix : bool) :
	classes_transcripts = makeDensityLaw() ;
	reference = "" ;
	library = "" ;
	all_features = "" ;
	with open(gff_file,"r") as gff :
		ligne = gff.readline().rstrip() ;
		
		i = 0 ;
		num_transcript = 0 ;
		while ligne :
			if ligne.find("\tgene\t") != -1 and ligne.find("=protein_coding;") != -1 :
				if i in choosen_genes :
					introns = None ; exon = None ;
					gene = [] ;
					nb_mRNA = 0 ;
					ligne = gff.readline().rstrip() ;
					while ligne and (ligne.find("\tgene\t") == -1 and ligne.find("\tpseudogene\t") == -1 ) :
						if ligne.find("\ttranscript\t") != -1 or ligne.find("\ttRNA\t") != -1 or ligne.find("\trRNA\t") != -1:
							while ligne and (ligne.find("\tgene\t") == -1 and ligne.find("\tmRNA\t") == -1 ) :
								ligne = gff.readline().rstrip() ;
						if ligne.find("\tmRNA\t") != -1  :
							nb_mRNA += 1 ;
						gene.append(ligne) ;
						ligne = gff.readline().rstrip() ;
					transcript = "" ;
					num_transcript += 1 ;
					mRNA,nb_exons = parseGene(gene,nb_mRNA,num_transcript) ;
					classe = int(rd.choice(classes_transcripts, 1)) ;
					if classe == -1 :
						fragments,exon,real_class = constructTranscript(list(tuple(mRNA)),nb_exons, classe) ;
						library += "\n".join(fragments) + "\n" ;
						transcript += "\n".join(addClasse(mRNA,real_class)) ;
						if mix :
							mRNA_copy = [] ;
							for ele_fea in list(mRNA) :
								new = ele_fea.split("\t") ;
								attr = new[-1].split(";") ;
								attr[1] += ".1" ;
								new[-1] = ";".join(attr) ;
								mRNA_copy.append("\t".join(new)) ; 
							library += "\n".join(mRNA_copy) + "\n" ;
						transcript += "\n".join(addClasse(mRNA,real_class)) ;
					elif classe == 0 :
						library += "\n".join(addClasse(mRNA,classe)) + "\n"; 
						transcript = "\n".join(addClasse(mRNA,classe)) ;
					else :
						fragments,introns,real_class = constructTranscript(list(tuple(mRNA)), nb_exons, classe) ;
						transcript = "\n".join(fragments) ;
						library += "\n".join(addClasse(mRNA,real_class)) + "\n"; 
						if mix :
							mRNA_copy = [] ;
							for ele_fea in list(fragments) :
								new = ele_fea.split("\t") ;
								attr = new[-1].split(";") ;
								attr[1] += ".1" ;
								new[-1] = ";".join(attr) ;
								mRNA_copy.append("\t".join(new)) ; 
							library += "\n".join(mRNA_copy) + "\n" ;
					reference += transcript +"\n";
					if introns :
						all_features += "\n".join(introns) + "\n" ; 
					elif exon :
						all_features += "\n".join(exon) + "\n" ;
				
				else :
					ligne = gff.readline().rstrip() ;
				i += 1 ;
				
			else :
				ligne = gff.readline().rstrip() ;
				
	writeGFF(all_features,output_path + "_Features_of_interest") ;
	writeGFF(reference,output_path + "_reference_transcripts.tmp") ;
	writeGFF(library, output_path + "_library_transcripts.tmp") ;


def makeDensityLaw() :
	law = [] ;
	config_path = os.path.abspath(os.path.dirname(sys.argv[0]) + "/../config/intronStalker.properties")

	config = configparser.RawConfigParser() ;
	config.read(config_path) ;
	for classe in config["Density"] :
		effectif = int(config["Density"][classe]) ;
		law += [int(classe)]*effectif ;
	return tuple(law) ;


def parseGene(gene : list, nb_mRNA : int ,num_transcript : int) :
	# We randomly pick one mRNA
	choosen_mRNA = int(rd.choice(range(1,nb_mRNA+1),1)) ;
	j = 0 ;
	current_mRNA = [] ;
	nb_exons = 0 ;
	while j <= choosen_mRNA :
		ligne = gene.pop(0) ;
		if ligne.find("\tmRNA\t") != -1 or not gene:
			j += 1 ;
			old_mRNA = current_mRNA ;
			old_nb = nb_exons ;
			current_mRNA = [] ;
			nb_exons = 0 ;
		elif ligne.find("\texon\t") != -1  :
			nb_exons += 1 ;
			exon = changeAttributes(ligne, str(num_transcript), nb_exons) ;
			current_mRNA.append(exon) ;
	return old_mRNA, old_nb ;


def changeAttributes(exon : str, num_transcript : str, num_exon : int) :
	transcript_id = "T"+num_transcript ;
	exon_id = "exon"+str(num_exon) ;
	items = exon.split("\t") ;
	old_attributes = items[-1].split(";") ;
	new_attributes = "ID="+exon_id+";Parent="+transcript_id+";"+(";".join(old_attributes[2:])) ;
	return "\t".join(items[:-1]) + "\t" + new_attributes ;


def addClasse(mRNA : list, classe : int) :
	for e in range(len(mRNA)) :
		mRNA[e] += ";transcript_class="+str(classe) ;
	return mRNA ;


def constructTranscript(mRNA : list, nb_exons : int, classe : int) :
	features_interest=[] ;
	if classe == -1 and nb_exons > 2:
		spliced_exon = int(rd.choice(range(1,nb_exons-1),1)) ;
		mRNA,exon = spliceExon(mRNA,spliced_exon) ; 
		features_interest.append(exon) ;
	elif classe >= 1 :
		nb_introns = nb_exons -1 ;
		if nb_introns <= classe :
			retained_introns = range(nb_introns) ;
			classe = nb_introns ;
		else :
			retained_introns = rd.choice(range(nb_introns), classe,replace=False) ;
		# We'll take each retained intron in decreasing order to prevent exon index error
		for intron in sorted(retained_introns,reverse=True) :
			coord_start = calcCoord(mRNA[0:intron+1])
			mRNA[intron],feature_interest = keepIntron(mRNA[intron],mRNA[intron+1],coord_start) ;
			del mRNA[intron+1] ;
			features_interest.append(feature_interest) ;
	else :
		classe = 0 ;
	return addClasse(mRNA,classe),features_interest,classe ;


def keepIntron(exon1 : str, exon2 : str , coord_start : int) :
	items1 = exon1.split("\t") ;
	items2 = exon2.split("\t") ;
	
	# We check the strand 
	if items1[6] == "+" :
		# New exon's start is the first exon's start
		new_start = items1[3] ;
		intron_start = items1[4] ;
		# New exon's end is the second exon's end
		intron_end = items2[3] ;
		new_end = items2[4] ;
	else :
		# New exon's start is the second exon's start
		new_start = items2[3] ;
		intron_start = items2[4] ;
		# New exon's end is the first exon's end
		intron_end = items1[3] ;
		new_end = items1[4] ;
	
	intron_length = int(intron_end) - 1 - (int(intron_start) + 1) +1 ;
	coord_end = coord_start + intron_length -1

	# We update the new exon's attributes
	attributes1 = items1[-1].split(";") ;
	attributes2	= items2[-1].split(";") ;
	new_ID = attributes1[0] + "+" + attributes2[0].lstrip("ID=") ;
	new_attributes = ";".join([
		new_ID,
		attributes1[1],
		";".join(attributes1[2:])
		]) ;

	new_exon = "\t".join([
		items1[0],
		items1[1],
		"exon",
		new_start,
		new_end,
		".",
		items1[6],
		".",
		new_attributes
		]) ;
	
	intron = "\t".join([
		attributes1[1].lstrip("Parent="),
		".",
		"retained_intron",
		str(coord_start),
		str(coord_end),
		".",
		items1[6],
		".",
		";".join([
			"ID=intron" + attributes1[0].lstrip("ID=exon"),
			"length="+str(intron_length)
		])
		]) ;
	
	return new_exon, intron ; 


def spliceExon(mRNA : list, spliced_exon : int) :
	feature = mRNA.pop(spliced_exon) ;
	items = feature.split("\t") ;
	attributes = items[-1].split(";") ;

	new_start = calcCoord(mRNA[0:spliced_exon]) ;
	exon_length = int(items[4]) - int(items[3]) +1
	new_end = int(new_start) + exon_length -1 ;

	exon = "\t".join([
		attributes[1].lstrip("Parent="),
		".",
		"spliced_exon",
		str(new_start),
		str(new_end),
		".",
		".",
		".",
		";".join([attributes[0],"length="+str(exon_length)])
		]) ;
	return mRNA,exon ;


def calcCoord(features_list : list) :
	total_length = 0 ;
	for feature in features_list :
		items = feature.split("\t") ;
		start,end = items[3],items[4] ;
		total_length += int(end) - int(start) +1 ;

	return total_length + 1 ;
	


def annoToData(annotation : str, fasta : str, nb : int, output : str, grinder : bool, mix : bool) :



	outdir = output + "_annoToData" ; # Output directory name (where this function will write 
								 # all the results and tmp files).
	os.system("mkdir " + outdir) ;
	output_path = outdir + "/" + output ;

	choosen_genes = chooseGenes(annotation, int(nb)) ;
	generateTranscripts(annotation, choosen_genes, output_path, mix) ;
	extractFasta(fasta,output_path) ;

	if grinder :
		profile_path = os.path.abspath(os.path.dirname(sys.argv[0]) + "/../config/profile_file.txt")
		_grinder(output_path+"_library.fa", profile_path, output_path)
		split_read(output_path)
		os.system("rm {path}-reads.fa".format(path=output_path)) ;
		os.system("rm {path}_library.fa".format(path=output_path)) ;

