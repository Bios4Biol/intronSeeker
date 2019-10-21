#!/usr/bin/env python3

import re
import pickle
import os
import pandas as pd
import pysam
import gzip
import time
import numpy as np
from pprint import pprint
from Bio import SeqIO
from collections import OrderedDict
import configparser
import concurrent.futures as prl
from itertools import repeat

def parse_config() :
    config = configparser.RawConfigParser()
    config.read('jupyter.properties')
    return config
    
def parse_positions(fastq_pos) :
    pos = fastq_pos.lstrip("position=").split("..")
    complement = ('complement(' in pos[0])
    start = int(pos[0].lstrip("complement("))-1
    end = int(pos[1].rstrip(")"))
    return start,end,complement

def parse_library(library_prefix) :
    lectures=[]
    with gzip.open(library_prefix+"read_1.fastq.gz","rt") as file1 :
        for record in SeqIO.parse(file1, "fastq") :
            reference = record.description.split()[1].lstrip("reference=")
            id = record.id
            start,end,complement =  parse_positions(record.description.split()[2])
            lectures.append([reference,id,start,end,complement])
    with gzip.open(library_prefix+"read_2.fastq.gz","rt") as file2 :
        for record in SeqIO.parse(file2, "fastq") :
            reference = record.description.split()[1].lstrip("reference=")
            id = record.id
            start,end,complement =  parse_positions(record.description.split()[2])
            lectures.append([reference,id,int(start),int(end),complement])
    return pd.DataFrame(lectures,columns=["contig","lecture","start","end","complement"]).sort_values(["contig","start","end"]).set_index('lecture') 

def parse_fasta(fastafile) :
    
    with open(fastafile,"r") as ff :
        fasta = {record.id : pd.Series({
            'length':len(record),
            'sequence' : record.seq,
            **{a.split("=")[0]:a.split("=")[1] for a in record.description.split() if a.startswith("classe")}
            })
        for record in SeqIO.parse(ff, "fasta")}
        contigs = pd.DataFrame.from_dict(fasta,orient='index')
        contigs.index.name='contig'
    return contigs

def parse_control_introns(introns_coord_file) :
    table = pd.read_table(introns_coord_file)
    table["length"] = table["end"]-table["start"]
    table['intron'] = table.apply(lambda df : "|".join([df.contig,str(df.start),str(df.end)]),axis=1)
    return table.set_index('intron')

def compute_tr_length(contig,features) :
    return contig.length -  features.loc[lambda df : df.contig == contig.name,"length" ].sum()

def compute_pos_on_contig(intron,contigs) :
    pos_on_contig = intron.start/contigs.at[intron.contig,"transcript_length"]*100
    c_seq = str(contigs.at[intron.contig,'sequence'])
    flanks = str(c_seq[intron.start:intron.start+2])+"_"+str(c_seq[intron.end-2:intron.end])
    return pd.Series([flanks,pos_on_contig],index=["flanks","pos_on_contig"])

def compute_pos_on_read(cov_lect,intron_start):
    if not cov_lect.complement :
        return (intron_start - cov_lect.start)/(cov_lect.end-cov_lect.start)*100
    else :
        return (cov_lect.end - intron_start)/(cov_lect.end-cov_lect.start)*100

def process_intron(intron,lectures) :
    
    cov_lect = pd.DataFrame(lectures.loc[lambda df : 
                         (df.contig == str(intron.contig)+".ori")
                         & (intron.start > df.start)
                         & (intron.start < df.end)
                          ])
    cov_lect['covering'] = True
    cov_lect['intron'] = intron.name
    cov_lect['pos_on_read'] = cov_lect.apply(
            compute_pos_on_read,
            axis=1,
            intron_start=intron.start
            )
    return cov_lect

def prlz_process_intron(df_introns,library) :
    reads = pd.concat(
        df_introns.apply(
            process_intron,
            axis=1,
            lectures=library
            ).values
        )
    return reads

def limit_from_cigar(cigar_list: list, start: int, ref_seq: str):
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
    split_start = start + limit_from_start[0]
    split_end = start + limit_from_start[0] + limit_from_start[1]
    length = limit_from_start[1]
    flank_left = ref_seq[split_start: split_start + 2]
    flank_right = ref_seq[split_end - 2: split_end]
    return pd.Series([int(split_start),int(split_end),length,flank_left+"_"+flank_right],
                    index = ["start_split","end_split","split_length","split_flanks"])

def process_bam(alignements, contigs, introns, library):
    """
    For an alignment file, list all split reads.

    :param fastafilename: reference fasta file
    :param bamfilename: AlignmentFile object with all reads alignment information
    :return: list of split. For each split, save its reference, name, start, stop, length and flanking sequences.
    """
    rows = []
    for record in alignements :
        begin = pd.Series([
            record['query_name'],                        # ID of the read
            record['reference_name'],                    # ID of the contig where the read is mapped
            record['reference_start'],                   # Start of the alignement on the contig
            record['reference_end'],                     # End of the alignement on the contig
            library.at[record['query_name'],"covering"], # Bool if the read normally covers an intron (i.e. should be split)
            record['cigartuples'] is not None,           # Bool if the read is mapped
            not record['reference_name'] == library.at[record['query_name'],"contig"].rstrip(".ori"), # Bool if the read is mapped on right contig
            (record['cigartuples'] is not None) and ('(3,' in str(record['cigartuples'])) # Bool if the read is split by aligner
            ],
            index = ["read","contig","align_start","align_end",'covering','mapped',"mismapped",'split']
        )
        
        
        if record['cigartuples'] is not None and '(3,' in str(record['cigartuples']) :
            
            row = begin.append(limit_from_cigar(
                record['cigartuples'], 
                record['reference_start'], 
                str(contigs.at[record['reference_name'],"sequence"])
            ))
            introns_to_check = introns.loc[lambda df : df.contig == record['reference_name'],:]
            for limits in zip(introns_to_check["start"],introns_to_check["end"]) :
                row["missplit"] = not (limits[0]-row.start_split in range(-3,4) and limits[1]-row.end_split in range(-3,4))
        else :
            row = begin.append(pd.Series([None,None,None,None,None],
                                         index=["start_split","end_split","split_length","split_flanks","missplit"]))
        if row.covering and row.mapped and not row.mismapped and row.split and not row.missplit :
            row["correct"] = True
        else :
            row["correct"] = None
        
        rows.append(row)
    
    return pd.DataFrame(rows)

def class_read(mapping) :
    if (mapping.covering and mapping.correct) or (mapping.split and not mapping.missplit) :
        return "TP"
    elif not mapping.covering and not mapping.split :
        return "TN"
    elif mapping.covering and not mapping.correct :
        return "FN"
    elif (not mapping.covering and mapping.split) or (mapping.covering and mapping.missplit):
        return "FP"

def class_multi(group) :
    items = set(group['classe'].values)
    if 'TP' in items :
        return 'TP'
    elif 'FP' in items :
        return 'FP'
    elif 'FN' in items :
        return 'FN'
    elif 'TN' in items :
        return 'TN'

def merge_split(contig_reads) :
    split_signal = contig_reads.sort_values(by=['start_split','end_split']).to_dict('index',OrderedDict)
    
    candidates=[]
    
    first = split_signal.popitem(last=False)
    current_candidate =[first[1]]
    curr_cand_starts = [first[1]['start_split']]
    curr_cand_ends = [first[1]['end_split']]
    to_check=OrderedDict()
    while split_signal or to_check :
        idx, split_read = split_signal.popitem(last=False)
        if split_read['start_split']-np.mean(curr_cand_starts) in range(-5,5) and split_read['end_split']-np.mean(curr_cand_ends) in range(-5,5) :
             print(split_read)
        

def alignement_analysis(args_dict) :

    print("Analysis "+args_dict["name"])
    
    path_rmpg_pick = "/".join([
        results_dir,
        "_".join(["reads_mapping",*args_dict["name"].split()])
        ])
    path_resreads_pick = "/".join([
        results_dir,
        "_".join(["res_lectures",*args_dict["name"].split()])
        ])
    if not os.path.exists(path_rmpg_pick) or not os.path.exists(path_resreads_pick) :
        # Reads library parsing and pickling 
        print("Lib parsing")
        path_lib_pick = '/'.join([
                results_dir,
                os.path.basename(args_dict['library'])
                ])+'.gz'
        if not os.path.exists(path_lib_pick) :
            library = parse_library(args_dict['library'])
            library.to_pickle(path_lib_pick)
        else :
            library = pd.read_pickle(path_lib_pick)
        
        # Reference Fasta parsing and pickling 
        print("contig parsing")
        path_ref_pick = '/'.join([
                results_dir,
                os.path.basename(args_dict['reference'])
                ])+'.gz'
        if not os.path.exists(path_ref_pick) :
            contigs = parse_fasta(args_dict['reference'])
            contigs.to_pickle(path_ref_pick)
        else :
            contigs = pd.read_pickle(path_ref_pick)
        
        # Introns control file parsing and pickling
        print("Control parsing")
        path_ctrl_pick = '/'.join([
                results_dir,
                os.path.basename(args_dict['control'])
                ])+'.gz'
        if not os.path.exists(path_ctrl_pick) :
            control = parse_control_introns(args_dict['control'])
            control.to_pickle(path_ctrl_pick)
        else :
            control = pd.read_pickle(path_ctrl_pick)
        
        # For each contig, we calculate the transcript length (i.e. only the exons total length)
        print("contig computation")
        contigs["transcript_length"] = contigs.apply(
                compute_tr_length,
                axis = 1,
                features=control
                )
        # For each intron, we calculate the insertion position in contig (in term of percentage of transcript length)
        # and also the "real" start of intron (i.e. insertion position in term of transcript coordinates - without introns
        print("introns computation")
        introns = control.join(
                    other = control.apply(
                        compute_pos_on_contig,
                        axis=1,
                        contigs=contigs
                    )
                )
        # For each intron, we determine all the reads which cover the insertion locus and the position in the read of this insertion locus
        # (precision : the function is called on intronns DataFrame but it returns a library-like DataFrame)
        print("library computation")
        with prl.ProcessPoolExecutor(max_workers=8) as ex :
            s_t = time.time()
            introns_split = np.array_split(introns,ex._max_workers)
            lectures = pd.concat(ex.map(prlz_process_intron,introns_split,repeat(library,ex._max_workers)))
            print(time.time()-s_t)
        lectures = library.join(lectures,lsuffix='',rsuffix='_cov').loc[:,lectures.columns]
        lectures.loc[lambda df : df.covering != True, "covering"] = False
        
        print("alignment computation")
        if not os.path.exists(path_rmpg_pick) :
            print("alignement parsing")
            bamfile = pysam.AlignmentFile(args_dict["bamfile"], "rb")
            alignements = [{
                        'query_name' : record.query_name,
                        'reference_name' : record.reference_name,
                        'reference_start' : record.reference_start,
                        'reference_end' : record.reference_end,
                        'cigartuples' : record.cigartuples
                        } for record in bamfile.fetch(until_eof=True)]
            
            print("alignement plrz")
            
            with prl.ProcessPoolExecutor(max_workers=8) as ex :
                s_t = time.time()
                align_split = np.array_split(alignements,ex._max_workers)
                reads_mapping = pd.concat(ex.map(process_bam,
                    align_split,
                    repeat(contigs,ex._max_workers),
                    repeat(introns,ex._max_workers),
                    repeat(lectures,ex._max_workers)))
                print(time.time()-s_t)
            print("duplicates computation")
            reads_mapping['multi_aligned'] = reads_mapping['read'].duplicated(keep=False)
            reads_mapping.to_pickle(path_rmpg_pick)
        else :
            reads_mapping = pd.read_pickle(path_rmpg_pick)
        
        print("lectures classification")
        reads_mapping['classe'] = reads_mapping.apply(class_read,axis=1)
        simple_aligned_classe = reads_mapping.loc[lambda df : df.multi_aligned != True,['read','classe']].set_index('read')
        multi_aligned_classe = pd.DataFrame(reads_mapping.loc[lambda df : df.multi_aligned == True,:].groupby('read').apply(class_multi),columns=['classe'])
        lectures = lectures.join(pd.concat([simple_aligned_classe,multi_aligned_classe]))
        
        lectures.to_pickle(path_resreads_pick)
    else :
        reads_mapping = pd.read_pickle(path_rmpg_pick)
        lectures = pd.read_pickle(path_resreads_pick)
    # ~ print()
    # ~ print("####### Results  ######")
    # ~ print(args_dict["name"])
    # ~ print(len(reads_mapping))
    # ~ print(len(lectures))
    # ~ print(lectures['classe'].value_counts())
    # ~ print("###################")
    # ~ print()
    
    reads_mapping.loc[lambda df : df.split == True,:].groupby('contig').apply(merge_split)
    # ~ print(len(np.array_split(reads_mapping.groupby('contig'),8)))
    
    return 
    
            
if __name__ == '__main__' :
    
    print("config parsing")
    config = parse_config()
    
    print("Res dir creation")
    results_dir = config["Global"]["env_dir"]
    if not os.path.exists(results_dir) :
        os.mkdir(results_dir)
        
    args_dicts = [config[analysis] for analysis in config['Global']['analysis'].split(',')]
    
    alignement_analysis(args_dicts[0])

    #~ start_time = time.time()
    #~ with prl.ProcessPoolExecutor(max_workers=4) as ex :
        #~ ex.map(alignement_analysis,args_dicts)
    #~ print('la totale : '+str(time.time() - start_time))
    
    
