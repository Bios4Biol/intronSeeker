#!/usr/bin/env python3

# Modules
try :
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Blast import NCBIXML as bx ;
	from collections import defaultdict
	import pandas as pd
	import pysam
	import os
	import re
	import subprocess as sp
except ImportError as error :
	print(error) ;
	exit(1) ;

def change_key(dictionary: dict, key_to_clean: dict):
	"""
	Change key dict name

	:param dictionary: dictionary to change keys
	:param key_to_clean: dictionary of new keys associate to previous keys
	:return: clean dictionary
	"""
	for clef, value in key_to_clean.items():
		if clef in dictionary:
			dictionary[value] = dictionary.pop(clef)
	return dictionary


#######################
# Extract split reads #
#######################

def limit_from_cigar(cigar_list: list, start: int, ref_name: str, read_name: str, reference: str):
	"""
	Write the intron extract from cigar line in file_r

	:param reference: reference sequence to extract flanking sequence of split read
	:param cigar_list: list of tuple (equal to the cigar line)
	:param start: beginning of the read alignment
	:param ref_name: reference name on which the read is aligned
	:param read_name: read name
	:return:
	"""
	values = [1, 0, 1, 1, 0]  # [M, I, D, N, S]
	limit_from_start = [0, 0]
	i = 0
	cigar_tuple = cigar_list[i]
	# if there is a split, calculate its position based on cigar line tuple
	while cigar_tuple[0] != 3:
		limit_from_start[0] += values[cigar_tuple[0]] * cigar_tuple[1]
		i += 1
		cigar_tuple = cigar_list[i]
	# enf of the split, equal to the number of 'N'
	limit_from_start[1] = cigar_tuple[1]
	split_limit = [start + limit_from_start[0], start + limit_from_start[0] + limit_from_start[1]]
	length = limit_from_start[1]
	flank_left = reference[split_limit[0]: split_limit[0] + 2]
	flank_right = reference[split_limit[1] - 2: split_limit[1]]
	return list(map(str, [ref_name, "read"+read_name, split_limit[0], split_limit[1], "length="+str(length), "flank_left="+flank_left, "flank_right="+flank_right]))

def split_research(bamfilename, fastafilename, basename) :
	"""
	Search the split reads and write two output files : the first with all the spliced events and the second where all the identical spliced events are merged
	:param bamfilename: name of the input alignment file
	:param fastafilename: name of the reference fasta file (contains the contigs) 
	:param basename: prefix for the name of output files 
	:return: nothing
	"""
	outdir = basename + "_split" ; # Output directory name (where this function will write 
								 # all the results and tmp files).
	os.system("mkdir " + outdir) ;
	output_path = outdir + "/" + basename ;

	res_gff,res_tmp = find_split(fastafilename,bamfilename) ;

	with open(output_path + "_all-hits.gff", "w") as file_r:
		file_r.write("\n".join(res_gff)+"\n")

	events_gff=group_event(res_tmp)

	with open(output_path + "_introns.gff","w") as file_r :
		file_r.write("\n".join(events_gff)+"\n")


def find_split(fastafilename: str, bamfilename):
	"""
	For an alignment file, list all split reads.

	:param fastafilename: reference fasta file
	:param bamfilename: AlignmentFile object with all reads alignment information
	:return: list of split. For each split, save its reference, name, start, stop, length and flanking sequences.
	"""
	reads_list = [x.split("\t")[0] for x in pysam.idxstats(bamfilename).split("\n")[:-2]]
	length_bam = len(reads_list)
	liste_id = reads_list
	bamfile = pysam.AlignmentFile(bamfilename, "rb")
	try :
		ref_dict = pysam.FastaFile(fastafilename)
	except OSError as e :
		if str(e).startswith("error when opening file") :
			print("\nIndexingError : it's impossible to write in Reference fasta directory to indexing file.")
			print("Please move the reference file in a directory where it's possible to write in and retry.\n")
		else :
			print(e)
		exit(1)
	list_result_gff = []
	list_result_tmp = []
	for id in liste_id:
		aligned = bamfile.fetch(id, multiple_iterators=True)
		for read in aligned:
			reference = ref_dict.fetch(id)
			if read.cigartuples is not None:

				if '(3,' in str(read.cigartuples):
					res_tmp = limit_from_cigar(read.cigartuples, read.reference_start, read.reference_name,read.query_name, reference)
					res_gff = res_tmp.copy()
					if read.is_reverse :
						res_gff.insert(4,"-") ;
					else :
						res_gff.insert(4,"+") ;
					list_result_gff.append(GFFformating(res_gff,"splicing_event"))
					list_result_tmp.append(res_tmp)
	return list_result_gff,list_result_tmp

def GFFformating(feature : list, feature_type : str) :
	attributes = feature[5:]
	gff_feature = [
		feature[0],
		feature[1],
		feature_type,
		feature[2],
		feature[3],
		".",
		feature[4],
		".",
		";".join(attributes)
		]
	return "\t".join(gff_feature)

def group_event(all_hits) :
	"""
	From a list of spliced events, merged all the events considered identical events
	:param all_hits: a list of all found events
	:return: return a list of feature (GFF formated)
	"""
	deja_vu={}
	a_creer=True
	for hit in all_hits :
		for event in deja_vu.keys() :
			if hit[0] == event[0] and (int(hit[2])-int(event[1]) in range(-2,3)) and (int(hit[3])-int(event[2]) in range(-2,3)) :
				deja_vu[event][-1] += 1
				a_creer=False
				break
		if a_creer :
			temp = hit
			temp.append(1)
			deja_vu[(hit[0],hit[2],hit[3])] = temp
		a_creer=True
	gff = []
	
	count = 1 ;
	for ligne in deja_vu.values() :
		ligne[1]="." ; ligne[-1] = "depth="+str(ligne[-1])
		ligne.insert(4,".")
		ID="event"+str(count)+"."+ligne[0]+"."+ligne[2]+"-"+ligne[3]
		ligne.insert(5,"ID="+ID)
		gff.append(GFFformating(ligne,"splicing_event"))
		count += 1 ;
	return gff


#~ def find_split(liste_id, fasta_file: str, bamfile):
	#~ """
	#~ For each subpart of a bigger alignment file, list all split reads.

	#~ :param liste_id: list of aligned reads from the alignment file
	#~ :param fasta_file: reference fasta file
	#~ :param bamfile: AlignmentFile object with all reads alignment information
	#~ :return: list of split. For each split, save its reference, name, start, stop, length and flanking sequences.
	#~ """

	#~ ref_dict = pysam.FastaFile(fasta_file)
	#~ list_result = []
	#~ for id in liste_id:
		#~ aligned = bamfile.fetch(id, multiple_iterators=True)
		#~ for read in aligned:
			#~ reference = ref_dict.fetch(id)
			#~ if read.cigartuples is not None:
				#~ if '(3,' in str(read.cigartuples):
						#~ list_result.append(limit_from_cigar(
							#~ read.cigartuples, read.reference_start, read.reference_name,
							#~ read.query_name, reference))
	#~ return "\n".join(list_result)


#~ def parallel(cpu, bam_file, fasta_file, split):
	#~ """
	#~ parallelisation of the split research

	#~ :param cpu: number of thread for the parallel processing
	#~ :param bam_file: alignment file
	#~ :param fasta_file: reference fasta file use for the alignment
	#~ :param split: output name
	#~ :return: write the split file
	#~ """
	#~ reads_list = [x.split("\t")[0] for x in pysam.idxstats(bam_file).split("\n")[:-2]]
	#~ length_bam = int(len(reads_list) / cpu)

	#~ liste_id = []
	#~ bamfile = pysam.AlignmentFile(bam_file, "rb")
	#~ for i in range(0, len(reads_list), length_bam):
		#~ liste_id.append(reads_list[i: i + length_bam])
	#~ find_splitx = partial(find_split, fasta_file=fasta_file, bamfile=bamfile)
	#~ p = Pool(cpu)
	#~ results = p.map(find_splitx, liste_id)

	#~ with open(split, "w") as file_r:
		#~ file_r.write("\n".join(results))

	#~ # create a html page with split read analyse

	#~ with open("temp", "w") as file_r:
		#~ file_r.write("\n".join([fasta_file, split]))

	#~ fct_path = os.path.dirname(os.path.realpath(__file__))
	#~ ipynb = os.path.join(fct_path, "filter_parameters.ipynb")
	#~ template = os.path.join(fct_path, "full.tpl")

	#~ os.system("jupyter nbconvert --to html --template %s --execute %s --output-dir=%s" % (template, ipynb, "./temp"))
	#~ os.remove("./temp")

	#~ data = pd.read_csv(split, sep="\t", names=['ref', 'query', 'qstart', 'qend', 'length',
											   #~ 'flank_left', 'flank_right'])
	#~ introns_grouped = data.groupby(['ref', 'qstart', 'qend', 'length', 'flank_left',
									#~ 'flank_right']).size().reset_index().rename(columns={0: 'depth'})
	#~ introns_grouped.to_csv(split, header=None, sep="\t", index=False)


######################################
# aggregate closed spliced intervals #
######################################

#~ def categorie(data, bornes):
	#~ """
	#~ Class spliced interval depending on its splicing sites

	#~ :param bornes: list of selected splicing sites
	#~ :param data: a line about a spliced intervals
	#~ :return:
	#~ """
	#~ bornes_data = (data['flank_left'], data['flank_right'])
	#~ if bornes_data in bornes:
		#~ return 'can'
	#~ else:
		#~ return 'ncan'


#~ def read_parameters(parameters: str):
	#~ """
	#~ Extract the parameters from a parameters yaml file

	#~ :param parameters: parameter yaml file
	#~ :return: minimal depth and splicing site list
	#~ """
	#~ with open(parameters, "r") as data:
		#~ param = data.read()
	#~ param = yaml.load(param)
	#~ bornes = []
	#~ for i in param['splicing_sites']:
		#~ pair = i.split("_")
		#~ pair = tuple(pair)
		#~ bornes.append(pair)
	#~ return param['depth'], bornes


#~ def aggregate(intron_file: str, distance: int, parameters: str):
	#~ """
	#~ Merged spliced intervals (5bp)

	#~ :param parameters:
	#~ :param distance:
	#~ :param intron_file:
	#~ :return: concatenate contigs file
	#~ """

	#~ depth, bornes = read_parameters(parameters)

	#~ # load and prepare data
	#~ split = pd.read_csv(intron_file, sep="\t", names=['ref', 'qstart', 'qend', 'length',
													  #~ 'flank_left', 'flank_right', 'count'])
	#~ split['borne'] = split.apply(categorie, axis=1, bornes=bornes)
	#~ split['used'] = True

	#~ # create the contig dictionary
	#~ contig_list = list(set(split[split['borne'] == 'can'].ref.values))
	#~ contig_dict = defaultdict(list)

	#~ for i in range(0, len(split)):
		#~ contig = split.iloc[i]
		#~ if contig['ref'] in contig_list:
			#~ contig_dict[contig['ref']].append(list(contig))

	#~ # merged closed intervals (less than 5bp distance)
	#~ resume = pd.DataFrame()
	#~ for i in contig_list:
		#~ sub_split = pd.DataFrame(contig_dict[i], columns=['ref', 'qstart', 'qend', 'length', 'flank_left',
														  #~ 'flank_right', 'count', 'borne', 'used'])
		#~ sub_split = sub_split.sort_values(by=['borne', 'count'], ascending=[True, False]).reset_index()
		#~ if len(sub_split) > 1:
			#~ for z in range(0, len(sub_split) - 1):
				#~ reference = sub_split.iloc[z].copy()
				#~ if reference['used'] and reference['borne'] == 'can':
					#~ subData = sub_split[sub_split['used']]
					#~ left = subData.qstart - reference.qstart
					#~ right = subData.qend - reference.qend
					#~ expression = (abs(left) <= distance) & (abs(right) <= distance)
					#~ to_concat = subData[expression]
					#~ list_index = list(to_concat.index)
					#~ reference['count'] = sum(to_concat['count'])
					#~ resume = resume.append(reference, ignore_index=True)
					#~ already_use = list(sub_split['used'])
					#~ for y in list_index:
						#~ already_use[y] = False
					#~ sub_split['used'] = already_use
			#~ subData = sub_split[(sub_split['borne'] == 'can') & (sub_split['used'])]
			#~ resume = resume.append(subData, ignore_index=True)
		#~ else:
			#~ resume = resume.append(sub_split, ignore_index=True)

	#~ # write the merged spliced intervals into a file
	#~ resume[['count', 'length', 'qstart', 'qend']] = resume[['count', 'length', 'qstart', 'qend']].astype(int)
	#~ resume = resume[['ref',  'qstart', 'qend', 'length', 'flank_left', 'flank_right', 'count']]
	#~ resume.to_csv("merged_" + intron_file, header=False, index=False, sep="\t")


#########################
# Overlapping intervals #
#########################

#~ def recouvrement(data: list):
	#~ """
	#~ Merge overlapping intervals into longer intervals

	#~ :param data: list of all intervals starts and ends on a contig
	#~ :return: list of extend intervals
	#~ """
	#~ list_debut = [interv[0] for interv in data]
	#~ list_fin = [interv[1] for interv in data]

	#~ list_debut.sort()
	#~ list_fin.sort()

	#~ list_intervalle_final = []
	#~ nb_superposition = 0
	#~ debut_intervalle_courant = 0

	#~ while list_debut:
		#~ ordre_debut_fin = (list_debut[0] > list_fin[0]) - (list_debut[0] < list_fin[0])
		#~ if ordre_debut_fin == -1:
			#~ pos_debut = list_debut.pop(0)
			#~ if nb_superposition == 0:
				#~ debut_intervalle_courant = pos_debut
			#~ nb_superposition += 1
		#~ elif ordre_debut_fin == +1:
			#~ pos_fin = list_fin.pop(0)
			#~ nb_superposition -= 1
			#~ if nb_superposition == 0:
				#~ nouvel_intervalle = (debut_intervalle_courant, pos_fin)
				#~ list_intervalle_final.append(nouvel_intervalle)
		#~ else:
			#~ list_debut.pop(0)
			#~ list_fin.pop(0)

	#~ if list_fin:
		#~ pos_fin = list_fin[-1]
		#~ nouvel_intervalle = (debut_intervalle_courant, pos_fin)
		#~ list_intervalle_final.append(nouvel_intervalle)

	#~ return list_intervalle_final


#~ def recoupement(intron_file: str, reference_file: str, parameters: str, intervals_name: str):
	#~ """
	#~ Estimate the complexity of spliced events by looking at the overlapping events

	#~ :param intervals_name: output file name
	#~ :param parameters: file of parameters for spliced events selection
	#~ :param intron_file: file of spliced events (query, limits, depth ...)
	#~ :param reference_file: reference file used for alignment
	#~ :return: dataframe of all spliced intervals, depth and how many different
	#~ spliced events are include in each intervals
	#~ """

	#~ depth, bornes = read_parameters(parameters)

	#~ data = pd.read_csv(intron_file, sep="\t", names=['ref', 'qstart', 'qend', 'length',
													 #~ 'flank_left', 'flank_right', 'count'])
	#~ data['borne'] = data.apply(categorie, axis=1, bornes=bornes)

	#~ data_can = data[(data['borne'] == 'can') & (data['count'] > depth)]
	#~ contig_list = list(set(data_can[data_can['borne'] == 'can'].ref.values))

	#~ contig_dict = defaultdict(list)
	#~ for i in range(0, len(data_can)):
		#~ contig = data_can.iloc[i]
		#~ if contig['ref'] in contig_list:
			#~ contig_dict[contig['ref']].append(list(contig))

	#~ intervals = pd.DataFrame()
	#~ for i in contig_list:
		#~ contig_split = pd.DataFrame(contig_dict[i], columns=['ref', 'qstart', 'qend', 'length', 'flank_left',
															 #~ 'flank_right', 'count', 'borne'])
		#~ split_list = contig_split[['qstart', 'qend']].values.tolist()
		#~ labels = ['contig', 'start', 'stop', 'number', 'overlap', 'depth']
		#~ if len(split_list) > 1:
			#~ list_intervals = recouvrement(split_list)
			#~ for interval in list_intervals:
				#~ sub_split = contig_split[
					#~ (interval[0] <= contig_split['qstart']) & (contig_split['qstart'] <= interval[1])]
				#~ number = len(sub_split)
				#~ overlap = int(sum(sub_split['length']) * 100 / (interval[1] - interval[0]))
				#~ depth = sum(sub_split['count'])
				#~ to_add = pd.DataFrame([[i, interval[0], interval[1], number, overlap, depth]], columns=labels)
				#~ intervals = intervals.append(to_add, ignore_index=True)
		#~ else:
			#~ to_add = pd.DataFrame([[i, contig_split['qstart'][0], contig_split['qend'][0], 1, 100,
									#~ contig_split['count'][0]]], columns=labels)
			#~ intervals = intervals.append(to_add, ignore_index=True)

	#~ intervals.to_csv(intervals_name, header=False, index=False, sep="\t")

	#~ fct_path = os.path.dirname(os.path.realpath(__file__))
	#~ ipynb = os.path.join(fct_path, "recoupement_graph.ipynb")
	#~ template = os.path.join(fct_path, "full.tpl")
	#~ temp = os.path.join(fct_path, "temp")

	#~ temp_path = os.path.dirname(os.path.realpath(intervals_name))
	#~ reference_path = os.path.realpath(reference_file)
	#~ with open(temp, 'w') as tmp:
		#~ to_write = [os.path.join(temp_path, os.path.basename(intervals_name)), reference_path]
		#~ tmp.write("\n".join(map(str, to_write)))

	#~ os.system("jupyter nbconvert --to html --template %s --execute %s --output-dir=%s" % (template, ipynb, temp_path))
	#~ os.remove(temp)


#########################
# write truncated fasta #
#########################

def _Truncate(fasta_file : str, gff_feature : str, output_path : str) :
	print("\nTruncate FASTA file : "+fasta_file+"\n");
	with open(fasta_file,"r") as input_f :
		list_seqs = input_f.read().lstrip(">").rstrip().split("\n>") ;
		seqs = { seq.split("\n",1)[0].split()[0] : seq.split("\n",1)[1].replace("\n","") for seq in list_seqs} ;
		with open(gff_feature,"r") as features_f :
			list_features = features_f.read().rstrip().split("\n") ;
			fasta_truncated = "" ;
			num_feature = 1 ;
			old_contig = "" ;
			for feature in list_features :
				items = feature.split("\t") ;
				if old_contig != items[0] :
					num_feature = 1 ;
				else :
					num_feature += 1 ;
				seq_truncated = seqs[items[0]][0:int(items[3])] + seqs[items[0]][int(items[4]):]+"\n"
				seq_head = ">"+items[0]+".e"+str(num_feature)+"\t"+items[0]+"\t"+items[3]+"-"+items[4]+"\t"+items[-1]+"\n"
				print(seq_head)
				fasta_truncated += seq_head + seq_truncated ;
				old_contig = items[0] ;
		with open(output_path + "-trq.fa","w") as fa :
			fa.write(fasta_truncated.rstrip()) ;

def truncate(fasta_file : str, gff_feature : str, output : str) :
	"""
	From a fasta file and a ggf file, write a new fasta file which contains expurgated sequences of provided features
	:param fasta_file: original sequences file
	:param gff_feature: gff file which contains the features to delete
	:param output: name of the output fasta file
	"""
	outdir = output + "_truncate" ;
	os.system("mkdir " + outdir) ;
	output_path = outdir + "/" + output ;
	
	_Truncate(fasta_file,gff_feature,output_path) ;



#######################
# search ORF on fasta #
#######################

def runTransDecoder(fasta_file : str, tmp_directory : str, refine : bool) :
	"""
	Run TransDecoder.
	:param fasta_file: Fasta file on which the ORF prediction will be perform
	:param output_directory: directory name where the TransDecoder files will be stored
	:return: nothing
	"""
	if refine :
		os.system("TransDecoder.LongOrfs -t {fasta} -O {tempdir} ; TransDecoder.Predict -t {fasta} -O {tempdir}".format(fasta=fasta_file,tempdir=tmp_directory)) ;
	else :
		os.system("TransDecoder.LongOrfs -t {fasta} -O {tempdir} ; TransDecoder.Predict -t {fasta} -O {tempdir} --no_refine_starts".format(fasta=fasta_file,tempdir=tmp_directory)) ;

def writeORFasTsv(transdecoder_orfs : str) : 
	"""
	From the TransDecoder output, write a tsv file which contains information about each predicted ORF.
	:param transdecoder_orfs: Transdecoder output which contains all the predicte ORFs (nucleic or proteic fasta format)
	:param output_tsv: name of the tsv ouput file
	:return: nothing
	"""
	print("Write TSV file from "+transdecoder_orfs) ;
	tsv="\t".join(["Sequence","Start","End","ORFid","length","type","strand","frame","score"]) ;
	grep_out = sp.getoutput("grep 'CDS' {filename} | cut -f 8,9".format(filename=transdecoder_orfs.replace(".pep",".gff3"))) ;
	frames = { cds.split("Parent=")[1] : cds.split("\t")[0] for cds in grep_out.split("\n") }
	with open(transdecoder_orfs,"r") as ORFfile :
		lines = [orf.split("\n",1)[0] for orf in ORFfile.read().lstrip(">").split("\n>")] ;
		for line in lines :
			vals = line.split() ;
			ORFid = vals[0] ; typ=vals[3].lstrip("type:") ; length=vals[4].lstrip("len:") ;
			strand= vals[5].split(",")[0].lstrip("(").rstrip(")") ;
			score = vals[5].split(",")[1].lstrip("score=") ;
			seq = vals[-1].split(":")[0] ;
			start = vals[-1].split(":")[1].split("-")[0] ; end = vals[-1].split(":")[1].split("-")[1].split("(")[0] ;
			tsv += "\n"+"\t".join([seq,start,end,ORFid,length,typ,strand,frames[ORFid],score]) ;
	with open(transdecoder_orfs+".cds.tsv","w") as tsvfile :
		tsvfile.write(tsv) ;

def predictORF(fasta_file, trunc, gff, rm, output : str,refine) :

	outdir = output + "_orf" ;
	os.system("mkdir " + outdir) ;
	output_path = outdir + "/" + output ;
	fasta_name = os.path.basename(fasta_file)
	

	runTransDecoder(fasta_file,output_path + "_intermediate",refine) ;
	os.system("mv *.transdecoder.* {outdir}".format(outdir=outdir)) ;
	writeORFasTsv(outdir +"/"+ fasta_name + ".transdecoder.pep") ;

	if trunc and gff:
		_Truncate(fasta_file, gff, output_path) ;
		runTransDecoder(output_path+"-trq.fa",output_path + "_intermediate-tronque",refine) ;
		os.system("mv *.transdecoder.* {outdir}".format(outdir=outdir)) ;
		writeORFasTsv(output_path+"-trq.fa.transdecoder.pep") ;

	if rm :
		os.system("rm pipeliner* ; rm -r {outdir}_intermediate* ".format(outdir=output_path))


#################################################
# Search long gap on contigs-proteins alignment #
#################################################

def runDiamond(fasta_file : str, protein_db : str, output : str,threads : int) :
	"""
	Run Diamond v.0.9.9 to align nucleic sequences against proteic database.
	:param fasta_file: Filename of the fasta which contains the nucleic sequences.
	:param protein_db: Name of the Diamond proteic database.
	:param output: Name of the output file.
	:param threads: Number of threads used for the alignment.
	"""
	os.system("diamond blastx -q {fasta} -d {proteindb} -o {output} -f 5 -p {threads} -e 0.05 --max-target-seqs 0".format(
		fasta=fasta_file,proteindb=protein_db, output=output+".tmp.xml",threads=threads)) ;

def searchGapsOfInterest(alignseq_subject : str , alignseq_query : str) :
	"""
	From query eand susject aligned sequences (extracted from Diamond XML output), find long gaps and return a list of them.
	:param alignseq_subject: String which contains the subject sequence.
	:param alignseq_query: String which contains the query sequence.
	:return list_gaps: List of long Gaps found in the alignment.
	"""
	gap_pattern = re.compile("-+") ;
	begin_idx = 0 ;
	gaps = gap_pattern.findall(alignseq_subject)
	list_gaps = [] ;
	if len(gaps) < 4 and len(gaps) > 0 :
		for gap_str in gap_pattern.findall(alignseq_subject) :
			if len(gap_str) > 5 :
				gap_start = alignseq_subject.index(gap_str,begin_idx) ;
				gap_length = len(gap_str) 
				
				# We want the gap coordinates on query sequence, so we calculate the total length of all the query gaps located before the gap of interest 
				translation = len("".join(gap_pattern.findall(alignseq_query,0,gap_length)))
				query_AA_start = gap_start-translation ;
				query_AA_end = gap_start-translation+gap_length-1 ;
				list_gaps.append((query_AA_start,query_AA_end)) ;
		return list_gaps

def parseDiamondXML(diamXml : str,output : str) :
	"""
	Parse Diamond XML output and write a GFF file which contains the Long gaps.
	:param diamXml: Name of the Diamond ouput file (XML format).
	:param output: Name of the ouput file.
	"""
	print("\nParsing of Diamond XML output to produce Long Gaps GFF \n")
	gff = [] ;
	align_gff = [] ;
	with open(diamXml) as resDiamond:
		blast_records = bx.parse(resDiamond) ;
		for record in blast_records :
			for align in record.alignments :
				num_hsp = 1 ;
				for hsp in align.hsps :
					contig = record.query ;
					align_start = hsp.query_start ;
					align_end = hsp.query_end ;
					frame = hsp.frame[0] ;
					if frame > 0 :
						strand = "+" ;
					else :
						strand = "-" ;
					prot_align = align.title.split()[0]+".hsp"+str(num_hsp) ;
					ID_align = contig.split()[0] +"_"+align.title.split()[0]+".hsp"+str(num_hsp) ;
					features_align = ";".join(["ID="+ID_align,"PROT="+prot_align,"LEN="+str(hsp.align_length),"BITSCORE="+str(hsp.bits),"EVAL="+str(hsp.expect)])
					align_gff.append("\t".join(map(str,[contig.split()[0],"Diamond","Alignment",align_start,align_end,".",strand,frame,features_align])))

					list_gaps = searchGapsOfInterest(hsp.sbjct,hsp.query)
					if list_gaps :
						num_gap = 1
						for gap in list_gaps :
							# We multiply the gap coordinates 3 times because we want the cordinates on nucleic sequence.
							start = gap[0]*3 ;
							end = gap[1]*3 ;
							prot = align.title.split()[0]+".gap"+str(num_gap) ;
							ID = contig.split()[0] +"_"+align.title.split()[0]+".gap"+str(num_gap) ;
							features = ";".join(["ID="+ID,"PROT="+prot,"LEN="+str(end-start),"BITSCORE="+str(hsp.bits),"EVAL="+str(hsp.expect)])
							gff.append("\t".join(map(str,[contig.split()[0],"Diamond","Gap",start,end,".",strand,frame,features])))
							num_gap += 1 ;
					num_hsp += 1
	with open(output+".long_gaps.gff","w") as ofile :
		ofile.write("\n".join(gff)+"\n") ;
	with open(output+".alignments.gff","w") as o2file :
		o2file.write("\n".join(align_gff)+"\n") ;

def searchProtein(fasta : str, dbprotein : str, output : str, threads : int, rm = False) :
	"""
	Run all the Protein analysis to search potential introns (i.e. long gaps beetween contigs and proteins
	:param fasta: Name of the fasta file which contains contigs.
	:param dbprotein: Name of the Diamond proteic database.
	:param output: Name of the ouput file.
	:param threads: Number of threads used for the alignment.
	:param rm : Boolean which rules the temporary files erasure.
	"""
	outdir = output + "_protein" ;
	os.system("mkdir " + outdir) ;
	path_to_output = outdir + "/" + output ;
	
	runDiamond(fasta, dbprotein, path_to_output, threads) ;
	parseDiamondXML(path_to_output+".tmp.xml", path_to_output) ;
	if not rm :
		os.system("rm {outdir}/*.tmp.*".format(outdir=outdir)) ;

#~ if __name__ == "__main__" :
	
	#~ parseDiamondXML("test_xml.tmp.xml", "test")
