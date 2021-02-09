#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import subprocess as sp # To run subprocess
import re      # To work on regular expression
import gzip    # To open gzip files R1 R2
from collections import OrderedDict   # To parse flagstat
from Bio import SeqIO   # To parse fasta file


def stat_from_gtf(gtf):
    """
    stat_from_gtf function open input GTF file than split lines in order to return 3 dictionnaries:
        nb_distinct_features : Number of distinct features from all GTF lines
        nb_ctg_by_feature    : Number of ctg by feature from all GTF lines (Ex: "Exon" see in X ctg, "Intron" see in Y ctg, ...)
        ctg_descr            : Number of features profiles by ctg (Ex: "1 Exon & 2 Intron" see in X ctg, "3 Introns" see in Y ctg, ...)
    """
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


def len_dist_from_gtf(gtf):
    """
    len_dist_from_gtf function open input GTF file than split lines in order to return 2 arrays:
        1 array of array with each len of each features from GTF file
        1 array of feature type
    """
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


def parse_fasta(fastafile, save_seq) :
    """
    parse_fasta function parse input fasta file and return pandas.DataFrame where each line is a seq description from FASTA file
    If save_seq is True, sequence is also returned in the pandas.DataFrame
    """
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

# Return read description (start, end, complement)
def parse_positions(fastq_pos) :
    pos = fastq_pos.lstrip("position=").split("..")
    complement = ('complement(' in pos[0])
    start = int(pos[0].lstrip("complement("))
    end = int(pos[1].rstrip(")"))
    return start,end,complement


def parse_library(r1, r2=0, mfasta=0) :
    """
    parse_library function parse library R1, and, if present, library R2 to return a pandas.DataFrame named library where each line is a read description.
    """
    if r1.endswith('.gz') :
        my_open = gzip.open
    else :
        my_open = open
    lectures=[]
    with my_open(r1,"rt") as file1 :
        for record in SeqIO.parse(file1, "fastq") :
            reference = record.description.split()[1].lstrip("reference=")
            id = record.id
            if mfasta:
                start,end,complement = parse_positions(record.description.split()[2])
            else:
                start = 0
                end = 0
                complement = 0
            lectures.append([reference, id, len(record), start, end, complement])
    if r2 :
        with my_open(r2,"rt") as file2 :
            for record in SeqIO.parse(file2, "fastq") :
                reference = record.description.split()[1].lstrip("reference=")
                id = record.id
                if mfasta:
                    start,end,complement = parse_positions(record.description.split()[2])
                else:
                    start = 0
                    end = 0
                    complement = 0
                lectures.append([reference, id, len(record), start, end, complement])
    return pd.DataFrame(lectures,columns=["contig","lecture","length","start","end","complement"]).sort_values(["contig","start","end"]).set_index('lecture') 

def parse_gtf(gtf) :
    """
    parse_gtf function parse input GTF file to return a panda which contains gtf features desc (seqref feature start end)
    """
    t = pd.read_table(gtf, usecols=[0,2,3,4], names=['contig','feature','start', 'end'], header=None)
    t["length"] = t["end"]-t["start"]
    t['features'] = t.apply(lambda df : "|".join([str(df.contig),str(df.start),str(df.end)]),axis=1)
    return t.set_index('features')

# def parse_control_introns(introns_coord_file) :
#     """
#     not used ?
#     """
#     table = pd.read_table(introns_coord_file, usecols=[0,3,4], names=['contig','start', 'end'], header=None)
#     table["length"] = table["end"]-table["start"]
#     table['intron'] = table.apply(lambda df : "|".join([df.contig,str(df.start),str(df.end)]),axis=1)
#     return table.set_index('intron')    
    

def parse_candidat(candidat) :
    """
    parse_candidat function parse input candidat file to return 1 dataframe and 2 int values:
        panda which contains candidats desc 
        mindepth (int)
        maxlen (int)
    """
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
    
    t = pd.read_table(candidat, usecols=[0,1,2,3,4,5,6], names=['ID', 'reference', 'start', 'end', 'depth','split_borders', 'filter'], skiprows=skip_rows) 
    t['key'] = t['ID']
    return t.set_index('key'), mindepth, maxlen


def parse_split(split):
    """
    parse_split function parse input split file to return a panda dataframe which contains split desc 
    """
    t = pd.read_table(split, usecols=[0,1,4,5,6], names=['reference', 'read', 'split_length', 'split_borders', 'strand'],  header=0)
    return t.set_index('read')

# Return int : nbreads, mapped, paired, proper
def parse_flagstat(flagstat) :
    """
    parse_flagstat function parse input flagstat file to return 8 int values:
        nbreads       : Number of lines in the BAM
        mapped        : Number of mapped
        mappercent    : Percentage of mapped
        paired        : Number of reads
        proper        : Number of properly paired
        properpercent : Percentage of properly paired
        secondary     : Secondary
        singletons    : Singletons
    """
    with open(flagstat) as f:
        mylist = [line.rstrip('\n') for line in f]
        for line in mylist:
            v = re.search('^(\d+)', line)
            if v:   
                if "QC-passed reads" in line:
                    nbreads = v.group(1)
                if "mapped (" in line:
                    mapped = v.group(1)
                    v = re.search('mapped \((\d+.\d+)\%', line)
                    mappercent = 0
                    if v:
                        mappercent = v.group(1)
                if "paired in sequencing" in line:
                    paired = v.group(1)
                if "properly paired" in line:
                    proper = v.group(1)
                    v = re.search('properly paired \((\d+.\d+)\%', line)
                    properpercent = 0
                    if v:
                        properpercent = v.group(1)
                if "secondary" in line:
                    secondary = v.group(1)
                if "singletons (" in line:
                    singletons = v.group(1)  
    return nbreads, mapped, mappercent, paired, proper, properpercent, secondary, singletons


def compute_tr_length(df_mfasta, df_features) :
    """
    compute_tr_length function extract "fasta" length from df_features dataframe 
    to add a column to df_fasta dataframe (without any simulated features).
    """
    return df_mfasta.length - df_features.loc[lambda df : df.contig == df_mfasta.name,"length" ].sum()


def compute_pos_on_mfasta(df_features, df_mfasta) :
    """
    compute_pos_on_mfasta function extract from df_mfasta dataframe 2 informations:
        1- the true insertion position of the simulated feature (in term of mfasta length percentage)
        2- the borders of the simulated features (in term of nucleotides)
    as 2 columns added to df_features dataframe.
    """
    pos_on_contig = df_features.start/df_mfasta.at[df_features.contig,"short_length"]*100
    c_seq = str(df_mfasta.at[df_features.contig,'sequence'])
    flanks = str(c_seq[df_features.start-1:df_features.start+1])+"_"+str(c_seq[df_features.end-2:df_features.end])
    
    return pd.Series([flanks,pos_on_contig],index=["flanks","pos_on_contig"])
        

def parse_rank_file(rank_file) :
    """
    parse_rank_file function parse input rank_file to return ranks dataframe
    """
    with open(rank_file,"r") as rf :
        for line in rf.read().rstrip().split("\n") :
            if line.startswith('#') :
                names = line.lstrip("# ").split("\t")
                names[1] = 'contig'
                ranks = []
            else :
                ranks.append(line.split("\t"))
    return pd.DataFrame(data=ranks, columns=names).set_index(names[1]).sort_index()


def run_assemblathon(fasta : str ) :
    """
    run_assemblathon function run assemblathon_stats.pl script then parse assemblathon result(s) to return statistics (INT format) :
        nbContigs      : Number of contigs
        totContigSize  : Total size of contigs
        longestContig  : Longest contig
        shortestContig : Shortest contig
        nbContigsSup1K : Number of contigs > 1K nt
        n50            : N50 contig length
        l50            : L50 contig count
        meanContigSize : Mean contig size
    """
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
        if "Number of contigs  " in  i:
            nbContigs = [int(w) for w in i.split() if w.isdigit()][0]
        elif "Total size of contigs" in i:
            totContigSize = [int(w) for w in i.split() if w.isdigit()][0]
        elif "Longest contig" in i:
            longestContig = [int(w) for w in i.split() if w.isdigit()][0]
        elif "Shortest contig" in i:
            shortestContig = [int(w) for w in i.split() if w.isdigit()][0]
        elif "Number of contigs > 1K nt" in i:
            nbContigsSup1K = [int(w) for w in i.split() if w.isdigit()][0]
        elif "N50 contig length" in i:             
            n50 = [int(w) for w in i.split() if w.isdigit()][0]
        elif "L50 contig count" in i:      
            l50 = [int(w) for w in i.split() if w.isdigit()][0]
        elif "Mean contig size" in i:               
            meanContigSize = [int(w) for w in i.split() if w.isdigit()][0]
        popen.wait()
    if popen.returncode != 0:
        raise RuntimeError('Error: Assemblathon output is empty.')

    return nbContigs, totContigSize, longestContig, shortestContig, nbContigsSup1K, n50, l50, meanContigSize


def split_int(number, separator=' ', count=3):
    """
    split_int function format int number to return a str formatted by 3 numbers. Example : 1 234 instead of 1234.
    Not apply for float numbers.
    """
    l = str(number).find('.')
    if l == -1:
        return str(separator.join(
            [str(number)[::-1][i:i+count] for i in range(0, len(str(number)), count)]
            )[::-1])
    else:
        return str(number)


def compute_dp(df_features, df_library) :
    """
    compute_dp function add DP to df_features using df_library
    """
    ref = df_features['contig'].replace('.modif','')
    return len(df_library.loc[lambda df : (df['contig'] == ref) & (df_features['start'] > df['mstart']) & (df_features['start'] < df['mend'])])


def compute_len(df_features, df_mfasta) :  
    """
    compute_len function add contig length to df_features using df_mfasta dataframe
    """
    ref = df_features['contig']
    return (df_mfasta.loc[ref, 'length'])
