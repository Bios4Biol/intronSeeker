#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import subprocess as sp # To run subprocess
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
    
# Return :
# panda which contains candidats desc 
# mindepth (int)
# maxlen (int)
def parse_candidat(candidat) :
    mindepth  = 0 # Extract min depth from candidates.txt file
    maxlen    = 0 # Extract max len from candidates.txt file
    skip_rows = 0 # Remove commented first 2 lines with mindepth and maxlen
    with open(candidat,"r") as fi:
        for ln in fi:
            if ln.startswith("##mindepth:"):
                mindepth=ln.split(":")[1].rstrip()
            elif ln.startswith("##maxlen:"):
                maxlen=ln.split(":")[1].rstrip()
            if ln.startswith("#"):
                skip_rows += 1
            else:
                break
    
    t = pd.read_table(candidat, usecols=[0,1,2,3,4,5,6], names=['ID', 'reference', 'start', 'end', 'depth','split_borders', 'filter'], skiprows=skip_rows)   #header=0 to remove commented header
    t['key'] = t['ID']
    return t.set_index('key'), mindepth, maxlen

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
# Return ranks dataframe
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

# Run assemblathon_stats.pl script then
# Return statistics (INT format) from assemblathon file(s) :
# nbContigs      : Number of contigs
# totContigSize  : Total size of contigs
# longestContig  : Longest contig
# shortestContig : Shortest contig
# nbContigsSup1K : Number of contigs > 1K nt
# n50            : N50 contig length
# l50            : L50 contig count
# meanContigSize : Mean contig size
def run_assemblathon(fasta : str ) :
    popen = sp.Popen('assemblathon_stats.pl '+ fasta, stdout = sp.PIPE, shell = True, encoding = 'utf8')
    reader = popen.stdout.read()
    nbContigs=0
    totContigSize=0
    longestContig=0
    shortestContig=0
    nbContigsSup1K=0
    n50=0
    l50=0
    meanContigSize=0
    res = [x.strip() for x in reader.split('\n')]
    for i in res:
        if "Number of contigs       " in  i:
            nbContigs=int(i.split("       ")[1].rstrip())
        elif "Total size of contigs    " in i:
            totContigSize=int(i.split("    ")[1].rstrip())
        elif "Longest contig" in i:
            longestContig=int(i.split("       ")[1].rstrip())
        elif "Shortest contig" in i:
            shortestContig=int(i.split("        ")[1].rstrip())
        elif "Number of contigs > 1K nt        " in i:  
            nbContigsSup1Ktmp=i.split("        ")[1].rstrip()
            nbContigsSup1K=int(re.sub(r'([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%)', r" ", nbContigsSup1Ktmp))
        elif "N50 contig length" in i:             
            n50=int(i.split("       ")[1].rstrip())
        elif "L50 contig count" in i:      
            l50=int(i.split("        ")[1].rstrip())
        elif "Mean contig size" in i:               
            meanContigSize=int(i.split("       ")[1].rstrip())
        popen.wait()
    if popen.returncode != 0:
        raise RuntimeError('Error: Assemblathon output is empty.')

    return nbContigs, totContigSize, longestContig, shortestContig, nbContigsSup1K, n50, l50, meanContigSize


# Return int formatted by 3 numbers. Example : 1 234 instead of 1234
def split_int(number, separator=' ', count=3):
    return separator.join(
        [str(number)[::-1][i:i+count] for i in range(0, len(str(number)), count)]
    )[::-1]

# Return string split by 3 characters
def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]

# Add DP to df_features using df_library
def compute_dp(df_features, df_library) :
    ref = df_features['contig'].replace('.modif','')
    return len(df_library.loc[lambda df : (df['contig'] == ref) & (df_features['start'] > df['start']) & (df_features['start'] < df['end'])])

# Add contig length to df_features using df_fasta
def compute_len(df_features, df_fasta) :
    ref = df_features['contig'].replace('.modif','')
    return (df_fasta.loc[ref, 'length'])


# Process df_features (dict), df_library (dict) and meanCoverage (int)
# Return :
# detectableIntrons (int)            : Number of dectable introns.
# TD - True Detectable (int)         : Number of same features found in candidats dataframe and df_cov_lect dataframe.
# detectablePreditNeg (int)          : Number of detectable introns predits negatives.
# nbFeaturesWithoutReads (int)       : Number of features covered by no read.
# nbIntronsWithReadsBelowCov (int)   : Number of reads below coverage
# https://stackoverflow.com/questions/23508351/how-to-do-a-conditional-join-in-python-pandas  
def process_intron(df_features : dict, df_library: dict, df_candidat:dict,meanCoverage :int): 
  
    # Add one column in df_feature with nb reads overlapping each intron/feature
    
    # library head :              contig  start  end  complement
    # lecture                                    
    # 80032/1   SEQUENCE1      0  101       False
    # 110561/1  SEQUENCE1      0  101       False 

    # features head :                                     contig          feature  start   end  length flanks  pos_on_contig
    # features                                                                                              
    # SEQUENCE1.modif|635|962    SEQUENCE1.modif  retained_intron    635   962     327  GT_AG      63.310070
    # SEQUENCE2.modif|941|1318   SEQUENCE2.modif  retained_intron    941  1318     377  GT_AG      76.317924
    
    # Add a new column in df_library with SEQ.modif
    df_library.columns = ['old_contig','start','end','complement']
    df_library['contig'] = (df_library['old_contig']+'.modif')
    # Select only columns we need
    df_library = df_library[['contig', 'old_contig','start', 'end']]
    df_features = df_features[['contig', 'feature', 'start', 'end']]

    # Merge 
    df_cov_lect = df_library.merge(df_features,on='contig', suffixes=('_lect', '_features'))
    # df_cov_lect                    contig   old_contig  start_lect  end_lect          feature  start_features  end_features
    # 0         SEQUENCE1.modif    SEQUENCE1           0       101  retained_intron             635           962
    # 1         SEQUENCE1.modif    SEQUENCE1           0       101  retained_intron             635           962
    # ...                   ...          ...         ...       ...              ...             ...           ...
    # 208696  SEQUENCE999.modif  SEQUENCE999         328       429  retained_intron             290          1017
    # 208697  SEQUENCE999.modif  SEQUENCE999         328       429  retained_intron             290          1017
    
    # And then select overlapping reads only
    df_cov_lect = df_cov_lect[((df_cov_lect.start_features > df_cov_lect.start_lect)  & (df_cov_lect.start_features < df_cov_lect.end_lect))]
    # 	contig	old_contig	start_lect	end_lect	feature	start_features	end_features
    # 333	SEQUENCE1.modif	SEQUENCE1	536	637	retained_intron	635	962
    # 334	SEQUENCE1.modif	SEQUENCE1	537	638	retained_intron	635	962

    # group row by contig in order to count nb reads by contig, filter with mean coverage.
    df_cov_lect['ID'] = df_cov_lect.apply(lambda df : "|".join([df.contig,str(df.start_features),str(df.end_features)]),axis=1)  # Add ID like df_candidat dataframe
    df_coverage = (df_cov_lect['ID'].value_counts()).to_frame()
    df_coverage.columns = ['coverage']
    df_coverage.index.name='ID'
    # ID	coverage
    # SEQUENCE564.modif|98|408	114
    # SEQUENCE703.modif|442|829	107
    # SEQUENCE763.modif|356|790	104print()
  
    # Filter on mean coverage 
    cond = (df_coverage['coverage'] > meanCoverage)
    detectableIntrons  = len(df_coverage.loc[cond]) 
    
    nbFeaturesWithoutReads = df_features.shape[0]- df_coverage.shape[0] # le nombre de "features" qui ne sont pas couvertes par des lectures 
    
    condBelowCov = (df_coverage['coverage'] <= meanCoverage)  # nombre de "features" qui sont couvertes par un nombre de lectures sous le seuil de profondeur.
    nbIntronsWithReadsBelowCov = len(df_coverage.loc[condBelowCov])
    
    df_detectables=df_coverage.loc[cond] # Dataframe of detectables introns
    # ID	coverage
    # SEQUENCE564.modif|98|408	114
    # SEQUENCE703.modif|442|829	107
    # SEQUENCE763.modif|356|790	104
    # SEQUENCE245.modif|88|643	103

    #https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html
    df_detectables = pd.merge(df_detectables, df_candidat, on='ID', suffixes=('_lect', '_candidat'))
    TP=df_detectables.shape[0]  #TP = features strictly identical between df_candidat and df_detectables
    condNoPASS = (df_detectables['filter'] != 'PASS')
    detectablePreditNeg=len(df_detectables.loc[condNoPASS])
         
    return detectableIntrons, TP, detectablePreditNeg, nbFeaturesWithoutReads, nbIntronsWithReadsBelowCov



# Return comparison between dataframes df_candidat and df_features
def candidatsVsFeatures(df_candidat:dict, df_features:dict, mindepth:int, maxlen:int):
    # Join candidat and features dataframe
    df_candidat = df_candidat.join(df_features,  lsuffix='_candidat', rsuffix='_features')
   
    # Total number of candidats including features
    nbTotCandidatsIncludingFeatures = df_candidat.shape[0]

    # Number of features with the same start and end than candidats
    conditionStartEnd = ((df_candidat['start_candidat'] == df_candidat['start_features']) & (df_candidat['end_candidat'] == df_candidat['end_features']))
    nbSameStartEnd=len(df_candidat.loc[conditionStartEnd]) - 1 # -1 for header
   
    # Features length >= 80 (default value in Split Read Search)
    condLen = ((((df_candidat['start_features'] - df_candidat['end_features'])/df_candidat['length'])*100) >= int(maxlen))
    nbLen   = len(df_candidat.loc[condLen]) 

    # Features with depth inf or equals to 1 (value by default)
    condDepth = (df_candidat['depth'] <= int(mindepth))
    minDepth = len(df_candidat.loc[condDepth])   

    # Number of features without canonical junctions
    condNoCanonical = ((df_candidat['split_borders'] != 'CT_AC') & (df_candidat['split_borders'] != 'GT_AG'))
    noCanonical = len(df_candidat.loc[condNoCanonical])
    
    return nbTotCandidatsIncludingFeatures, nbSameStartEnd, nbLen, minDepth, noCanonical