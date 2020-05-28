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
    t = pd.read_table(candidat, usecols=[0,1,2,3,4,5,6], names=['ID', 'reference', 'start', 'end', 'depth','split_borders', 'filter'],  header=0)   #header=0 to remove commented header
    t['key'] = t['ID']
    return t.set_index('key')

# Return comparison between dataframes df_candidat and df_features
# Return 5 stats:
# nbTotCandidatsIncludingFeatures : Total number of candidats including features
# nbSameStartEnd : Number of features with the same start and end than candidats
# nbLen : Features length >= 80 (default value in Split Read Search)
# minDepth : Features with depth inf or equals to 1 (value by default)
# nonCanonical : Number of features without canonical junctions
def candidatsVsFeatures(df_candidat, df_features):
    # Join candidat and features dataframe
    df_candidat = df_candidat.join(df_features,  lsuffix='_candidat', rsuffix='_features')
   
    # Total number of candidats including features
    nbTotCandidatsIncludingFeatures = df_candidat.shape[0]

    # Number of features with the same start and end than candidats
    conditionStartEnd = ((df_candidat['start_candidat'] == df_candidat['start_features']) & (df_candidat['end_candidat'] == df_candidat['end_features']))
    nbSameStartEnd=len(df_candidat.loc[conditionStartEnd]) - 1 # -1 for header
   
    # Features length >= 80 (default value in Split Read Search)
    condLen = ((((df_candidat['start_features'] - df_candidat['end_features'])/df_candidat['length'])*100) >= 80)
    nbLen   = len(df_candidat.loc[condLen]) 

    # Features with depth inf or equals to 1 (value by default)
    condDepth = (df_candidat['depth'] <= 1)
    minDepth  = len(df_candidat.loc[condDepth]) 

    # Number of features without canonical junctions
    condNonCanonical = ((df_candidat['split_borders'] != 'CT_AC') & (df_candidat['split_borders'] != 'GT_AG'))
    nonCanonical = len(df_candidat.loc[condNonCanonical])
    
    return nbTotCandidatsIncludingFeatures, nbSameStartEnd, nbLen, minDepth, nonCanonical

# Return panda which contains split desc 
def parse_split(split):
    t = pd.read_table(split, usecols=[0,1,4,5,6], names=['reference', 'read', 'split_length', 'split_borders', 'strand'],  header=0)
    return t.set_index('read')

# Return int : nbreads, mapped, paired, proper
def parse_flagstat(flagstat) :
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
            if "secondary" in line:
                secondary=line[0:pos2]
            if "singletons (" in line:
                singletons=line[0:pos2]   
    return nbreads, mapped, paired, proper, secondary, singletons

def compute_tr_length(df_mfasta, df_features) :
    return df_mfasta.length - df_features.loc[lambda df : df.contig == df_mfasta.name,"length" ].sum()


def compute_pos_on_mfasta(df_features, df_mfasta) :
    pos_on_contig = df_features.start/df_mfasta.at[df_features.contig,"short_length"]*100
    c_seq = str(df_mfasta.at[df_features.contig,'sequence'])
    flanks = str(c_seq[df_features.start:df_features.start+2])+"_"+str(c_seq[df_features.end-2:df_features.end])
    
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

# Return int formatted by 3 numbers. Example : 1 234 instead of 1234
def split_int(number, separator=' ', count=3):
    return separator.join(
        [str(number)[::-1][i:i+count] for i in range(0, len(str(number)), count)]
    )[::-1]

# Return string split by 3 characters
def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]