#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import re      # To work on regular expression
import gzip    # To open gzip files R1 R2
from collections import OrderedDict   # To parse flagstat
from Bio import SeqIO   # To parse fasta file



# Return 3 dict : nb_distinct_features, nb_ctg_by_feature, ctg_descr
def stat_from_gtf(gtf):
    nb_distinct_features = dict()   # Number of distinct features from all GTF lines
    nb_ctg_by_feature    = dict()   # Number of ctg by feature from all GTF lines (Ex: "Exon" see in X ctg, "Intron" see in Y ctg, ...)
    ctg_descr            = dict()   # Number of features profiles by ctg (Ex: "1 Exon & 2 Intron" see in X ctg, "3 Introns" see in Y ctg, ...)

    tmp_array = []
    for line in open(gtf):
        if not line.startswith('#'):
            k = line.split()[0]
            if k in ctg_descr:
                ctg_descr[k][line.split()[2]] = ctg_descr[k].get(line.split()[2], 0) + 1
                if(line.split()[2] not in tmp_array):
                    tmp_array.append(line.split()[2])
                    nb_ctg_by_feature[line.split()[2]] = nb_ctg_by_feature.get(line.split()[2], 0) + 1
            else:
                ctg_descr[k] = dict()
                ctg_descr[k][line.split()[2]] = ctg_descr[k].get(line.split()[2], 0) + 1
                tmp_array = []
                tmp_array.append(line.split()[2])
                nb_ctg_by_feature[line.split()[2]] = nb_ctg_by_feature.get(line.split()[2], 0) + 1
                
            nb_distinct_features[line.split()[2]] = nb_distinct_features.get(line.split()[2], 0) + 1
    
    res = []
    for ctg in ctg_descr.values():
        tmpstr = ""
        for k, v in sorted(ctg.items(), key=lambda t: t[0]):
            if(tmpstr != ""):
                tmpstr += " and " 
            tmpstr += str(v)+" "+str(k)
        res.append(tmpstr)
    unique_elements, counts_elements = np.unique(res, return_counts=True)
    ctg_descr = dict()
    for i, e in enumerate(unique_elements):
        ctg_descr[e] = ctg_descr.get(e, counts_elements[i])
    
    return nb_distinct_features, nb_ctg_by_feature, ctg_descr


# Return : 1 array of array with each len of each features from GTF file
#          1 array of feature type
def len_dist_from_gtf(gtf):
    feature_names = []
    len_by_features = []
    nbf = 0
    for line in open(gtf):
        if not line.startswith('#'):
            feature = line.split()[2]
            start   = line.split()[3]
            end     = line.split()[4]
            if(feature not in feature_names):
                feature_names.append(feature)
                len_by_features.append([int(end)-int(start)+1])
                nbf += 1
            else:
                len_by_features[feature_names.index(feature)].append(int(end)-int(start)+1)        
    return len_by_features, feature_names

# Parse fasta file and return pandas.DataFrame
def parse_fasta(fastafile, save_seq) :
    with open(fastafile,"r") as ff :
        if(save_seq) :
            fasta = {record.id : pd.Series({
                'length':len(record),
                'sequence':record.seq,
                **{a.split("=")[0]:a.split("=")[1] for a in record.description.split() if a.startswith("class")}
            })
            for record in SeqIO.parse(ff, "fasta")}
        else :
            fasta = {record.id : pd.Series({
                'length':len(record),
                **{a.split("=")[0]:a.split("=")[1] for a in record.description.split() if a.startswith("class")}
            })
            for record in SeqIO.parse(ff, "fasta")}
        df = pd.DataFrame.from_dict(fasta,orient='index')
        df.index.name='contig'
    return df


# Parse library R1 and R2 return a pandas.DataFrame named library where each line is a read description
def parse_library(r1, r2=0) :
    if r1.endswith('.gz') :
        my_open = gzip.open
    else :
        my_open = open
    lectures=[]
    with my_open(r1,"rt") as file1 :
        for record in SeqIO.parse(file1, "fastq") :
            reference = record.description.split()[1].lstrip("reference=")
            id = record.id
            start,end,complement =  parse_positions(record.description.split()[2])
            lectures.append([reference,id,start,end,complement])
    if r2 :
        with my_open(r2,"rt") as file2 :
            for record in SeqIO.parse(file2, "fastq") :
                reference = record.description.split()[1].lstrip("reference=")
                id = record.id
                start,end,complement =  parse_positions(record.description.split()[2])
                lectures.append([reference,id,int(start),int(end),complement])
    return pd.DataFrame(lectures,columns=["contig","lecture","start","end","complement"]).sort_values(["contig","start","end"]).set_index('lecture') 


# Return read description (start, end, complement)
def parse_positions(fastq_pos) :
    pos = fastq_pos.lstrip("position=").split("..")
    complement = ('complement(' in pos[0])
    start = int(pos[0].lstrip("complement("))-1
    end = int(pos[1].rstrip(")"))
    return start,end,complement


# Return panda which contains gtf features desc (seqref feature start end)
def parse_gtf(gtf) :
    t = pd.read_table(gtf, usecols=[0,2,3,4], names=['contig','feature','start', 'end'], header=None)
    t["length"] = t["end"]-t["start"]
    t['features'] = t.apply(lambda df : "|".join([df.contig,str(df.start),str(df.end)]),axis=1)
    return t.set_index('features')

def parse_control_introns(introns_coord_file) :
    table = pd.read_table(introns_coord_file, usecols=[0,3,4], names=['contig','start', 'end'], header=None)
    table["length"] = table["end"]-table["start"]
    table['intron'] = table.apply(lambda df : "|".join([df.contig,str(df.start),str(df.end)]),axis=1)
    return table.set_index('intron')    
    
# Return panda which contains candidats desc 
def parse_candidat(candidat) :
    t = pd.read_table(candidat, usecols=[0,2,3,4,5,6], names=['ID', 'start', 'end', 'depth','split_borders', 'filter'],  header=0)   #header=0 to remove commented header
    print(type(t))
    print(t.dtypes)
    return t.set_index('ID')

# Return panda which contains split desc 
def parse_split(split):
    t = pd.read_table(split, usecols=[0,1,4,5,6], names=['reference', 'read', 'split_length', 'split_borders', 'strand'],  header=0)
    print(t.head(5))
    return t.set_index('read')


# Return int : nbreads, mapped, paired, proper
def parse_flagstat2(flagstat) :
    with open(flagstat) as f:
        mylist = [line.rstrip('\n') for line in f]
        for i in range(0, 12):
            line=mylist[i]
            #pos1 = line.find('\D\s')
            pos2 = line.find('+')  
            if "QC-passed reads" in line:
                nbreads=line[0:pos2]
            if "mapped (" in line:
                mapped=line[0:pos2]
            if "paired in sequencing" in line:
                paired=line[0:pos2]
            if "properly paired" in line:
                proper=line[0:pos2]
    return nbreads, mapped, paired, proper



def parse_flagstat(filename : str , lib_size : int, name : str) :
    flagstat = OrderedDict({})
    with open(filename,"r") as f :
        
        first = f.readline().split(" ",3)
        if int(first[2]) == 0 :
            qc_failed = False
            flagstat["Total count"] = int(first[0])
        else :
            qc_failed = True
            flagstat["Total count"] = [int(first[0]),int(first[2])]
            
        items_of_interest=["secondary",
                           "supplementary",
                           "duplicates",
                           "mapped",
                           "properly paired",
                           "singletons",
                           "with mate mapped to a different chr"]
        
        for ligne in f :
            values = ligne.rstrip().split(" ",3)
            if not qc_failed and not int(values[0]) == 0 :
                item = values[-1].split(" (")[0]
                if item in items_of_interest :
                    if item in ["mapped","properly paired","singletons"] :
                        flagstat[item] = "{value} ({percentage}%)".format(
                            value=values[0],
                            percentage=round((int(values[0])/(lib_size+flagstat["secondary"]))*100,2))
                    else :
                        flagstat[item] = int(values[0])
                    
            elif qc_failed and (not int(values[0]) == 0 or not int(values[2]) == 0) :
                item = values[-1].split(" (")[0]
                if item in items_of_interest :
                    if item in ["mapped","properly paired","singletons"] :
                        flagstat[item] = [
                            "{value} ({percentage}%)".format(
                                value=values[0],
                                percentage=round((int(values[0])/(lib_size+flagstat["secondary"]))*100,2)
                            ),
                            "{value} ({percentage}%)".format(
                                value=values[2],
                                percentage=round((int(values[2])/(lib_size+flagstat["secondary"]))*100,2)
                            )]
                    else :
                        flagstat[item] = [int(values[0]),int(values[2])]
            if item == "with mate mapped to a different chr" :
                break
        if not qc_failed :
            return pd.DataFrame.from_dict(flagstat,"index",columns=[name])
        else :
            return pd.DataFrame.from_dict(flagstat,"index",columns=pd.MultiIndex.from_tuples([(name,"QC-passed"),(name,"QC-failed")]))


def compute_tr_length(df_mfasta, df_features) :
    return df_mfasta.length - df_features.loc[lambda df : df.contig == df_mfasta.name,"length" ].sum()


def compute_pos_on_mfasta(df_features, df_mfasta) :
    #print(df_mfasta.at[df_features.contig,"short_length"]) #608
    pos_on_contig = df_features.start/df_mfasta.at[df_features.contig,"short_length"]*100
    #print('pos_on_contig',pos_on_contig) #pos_on_contig 42.26973684210527
    c_seq = str(df_mfasta.at[df_features.contig,'sequence'])
    #print('c_seq',c_seq)
    #c_seq TCAGGGCTCGAATAAACAGGCAAGCGGCTCGTAGATGGTGCTATCTTAACAACAAGGAAACGGCCCTGGATCGCCAGTTATACAAGGCGGAG...
    flanks = str(c_seq[df_features.start:df_features.start+2])+"_"+str(c_seq[df_features.end-2:df_features.end])
    #print('flanks',flanks) #flanks CT_AC
    
    return pd.Series([flanks,pos_on_contig],index=["flanks","pos_on_contig"])
        
# Parse ranks file
def parse_rank_file(rank_file) :
    with open(rank_file,"r") as rf :
        for line in rf.read().rstrip().split("\n") :
            if line.startswith('#') :
                names = line.lstrip("# ").split("\t")
                names[1] = 'contig'
                ranks = []
            else :
                ranks.append(line.split("\t"))
    return pd.DataFrame(data=ranks, columns=names).set_index(names[1]).sort_index()

#Useful functions :
def parsing_test(items) :
    items_of_interest = ["Number of contigs","Total size of contigs","Longest contig","Shortest contig","Number of contigs > 1K nt","N50 contig length","L50 contig count"]
    if len(items) == 2 and items[0] in items_of_interest :
        return True
    else :
        return False

# Parse assemblathon files to compare assembly with Assemblathon.pl statistics
def parse_assemblathon(filename : str, name : str ) :
    with open(filename,"r") as f :
        assemblathon = { re.split("\s\s+",line.strip(),1)[0] : re.split("\s\s+",line.strip(),1)[1] for line in f if parsing_test(re.split("\s\s+",line.strip(),1))}
    return pd.DataFrame(data=assemblathon.values(),index=assemblathon.keys(),columns=[name])    

# # Parse Alignment BAM files
# def parse_BAM(BamPath:str):    
#     bamfile = pysam.AlignmentFile(BamPath, "rb")
#     alignments = [{
#         'query_name' : record.query_name,
#         'reference_name' : record.reference_name,
#         'reference_start' : record.reference_start,
#         'reference_end' : record.reference_end,
#         'cigartuples' : record.cigartuples,
#         'is_secondary' : record.is_secondary,
#         'is_supplementary' : record.is_supplementary,
#         'mapping_quality' : record.mapping_quality
#         } for record in bamfile.fetch(until_eof=True)]   
    
#     return alignments

# def limit_from_cigar(cigar_list: list, start: int, ref_seq: str):
#     """
#     Write the intron extract from cigar line in file_r

#     :param reference: reference sequence to extract flanking sequence of split read
#     :param cigar_list: list of tuple (equal to the cigar line)
#     :param start: beginning of the read alignment
#     :param ref_name: reference name on which the read is aligned
#     :param read_name: read name
#     :return:
#     """
#     values = [1, 0, 1, 1, 0]  # [M, I, D, N, S]
#     limit_from_start = [0, 0]
#     i = 0
#     cigar_tuple = cigar_list[i]
#     # if there is a split, calculate its position based on cigar line tuple
#     while cigar_tuple[0] != 3:
#         limit_from_start[0] += values[cigar_tuple[0]] * cigar_tuple[1]
#         i += 1
#         cigar_tuple = cigar_list[i]
#     # enf of the split, equal to the number of 'N'
#     limit_from_start[1] = cigar_tuple[1]
#     split_start = start + limit_from_start[0]
#     split_end = start + limit_from_start[0] + limit_from_start[1]
#     length = limit_from_start[1]
#     flank_left = ref_seq[split_start: split_start + 2]
#     flank_right = ref_seq[split_end - 2: split_end]
#     return pd.Series([int(split_start),int(split_end),length,flank_left+"_"+flank_right],
#                     index = ["start_split","end_split","split_length","split_flanks"])

# #return: list of split. For each split, save its reference, name, start, stop, length and flanking sequences.
# def process_bam(alignments, df_mfasta, df_features, df_library):
#     """
#     For an alignment file, list all split reads.

#     :param fastafilename: reference fasta file
#     :param bamfilename: AlignmentFile object with all reads alignment information
#     :return: list of split. For each split, save its reference, name, start, stop, length and flanking sequences.
#     """
#     print("Enter process bam")
#     rows = []
#     for record in alignments :
#         begin = pd.Series([
#             record['query_name'],                        # ID of the read
#             record['reference_name'],                    # ID of the contig where the read is mapped
#             record['reference_start'],                   # Start of the alignment on the contig
#             record['reference_end'],                     # End of the alignment on the contig1
#             df_library.at[record['query_name'],"covering"], # Bool if the read normally covers an intron (i.e. should be split)
#             record['cigartuples'] is not None,           # Bool if the read is mapped
#             not record['reference_name'] == df_library.at[record['query_name'],"contig"],  # Bool if the read is mapped on right contig
#             (record['cigartuples'] is not None) and ('(3,' in str(record['cigartuples'])), # Bool if the read is split by aligner
#             record['is_secondary'],
#             record['is_supplementary'],
#             record['mapping_quality']]
#             ,
#             index = ["read","contig","align_start","align_end",'covering','mapped',"mismapped",'split','second','suppl','score']
#         )
        
#         if record['cigartuples'] is not None and '(3,' in str(record['cigartuples']) :
#             row = begin.append(limit_from_cigar(
#                 record['cigartuples'], 
#                 record['reference_start'], 
#                 str(df_mfasta.at[record['reference_name'],"sequence"])
#             ))
#             introns_to_check = df_features.loc[lambda df : df.contig == record['reference_name'],:]
#             for limits in zip(introns_to_check["start"],introns_to_check["end"]) :
#                 row["missplit"] = not (limits[0]-row.start_split in range(-3,4) and limits[1]-row.end_split in range(-3,4))
#         else :
#             row = begin.append(pd.Series([None,None,None,None,None],
#                                          index=["start_split","end_split","split_length","split_flanks","missplit"]))
#         rows.append(row)
    
#     return pd.DataFrame(rows)

# def compute_pos_on_read(cov_lect,intron_start):
#     if not cov_lect.complement :
#         return (intron_start - cov_lect.start)/(cov_lect.end-cov_lect.start)*100
#     else :
#         return (cov_lect.end - intron_start)/(cov_lect.end-cov_lect.start)*100

# # Return df_cov_lect, a new Dataframe which contains reads from df_library if INTRON cov (lecture,
# # contig, start, end, complement) join with 3 news columns:
# #   "covering"      : one for the intron covering reads (True/False)
# #   "intron (name)" : another for the covered intron id (if True) 
# #   "pos_on_read"   : intron insertion position in read (if True - in term of read length percentage)
# def process_intron(df_features,df_library) :
#     df_cov_lect = pd.DataFrame(df_library.loc[lambda df : 
#                          (df.contig+'.modif'== str(df_features.contig))
#                          & (df_features.start > df.start)
#                          & (df_features.start < df.end)
#                           ])
#     if df_cov_lect.shape[0] > 0 :
#         df_cov_lect['covering'] = True
#         df_cov_lect['intron'] = df_features.name
#         df_cov_lect['pos_on_read'] = df_cov_lect.apply(
#                 compute_pos_on_read,
#                 axis=1,
#                 intron_start=df_features.start
#                 )             
#     return df_cov_lect

# # Return DataFrame Reads
# def prlz_process_intron(df_features,df_library) :
#     df_reads = pd.concat(
#         df_features.apply(
#             process_intron,
#             axis=1,
#             df_library=df_library
#             ).values
#         )   
#     return df_reads

# Return int formatted by 3 numbers. Example : 1 234 instead of 1234
def split_int(number, separator=' ', count=3):
    return separator.join(
        [str(number)[::-1][i:i+count] for i in range(0, len(str(number)), count)]
    )[::-1]

# Return string split by 3 characters
def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]