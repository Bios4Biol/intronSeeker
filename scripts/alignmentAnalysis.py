#!/usr/bin/env python3

import re
import pickle
import os
import pandas as pd
import pysam
import gzip
import time
import sys
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


def process_bam(alignments, contigs, introns, library):
    """
    For an alignment file, list all split reads.

    :param fastafilename: reference fasta file
    :param bamfilename: AlignmentFile object with all reads alignment information
    :return: list of split. For each split, save its reference, name, start, stop, length and flanking sequences.
    """
    rows = []
    for record in alignments :
        begin = pd.Series([
            record['query_name'],                        # ID of the read
            record['reference_name'],                    # ID of the contig where the read is mapped
            record['reference_start'],                   # Start of the alignment on the contig
            record['reference_end'],                     # End of the alignment on the contig1
            library.at[record['query_name'],"covering"], # Bool if the read normally covers an intron (i.e. should be split)
            record['cigartuples'] is not None,           # Bool if the read is mapped
            not record['reference_name'] == library.at[record['query_name'],"contig"].rstrip(".ori"), # Bool if the read is mapped on right contig
            (record['cigartuples'] is not None) and ('(3,' in str(record['cigartuples'])), # Bool if the read is split by aligner
            record['is_secondary'],
            record['is_supplementary'],
            record['mapping_quality']]
            ,
            index = ["read","contig","align_start","align_end",'covering','mapped',"mismapped",'split','second','suppl','score']
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
        rows.append(row)
    
    return pd.DataFrame(rows)

def class_read(mapping) :
    if mapping.covering and mapping.split and not mapping.missplit :
        return "TP"
    if not mapping.covering and (not mapping.split or mapping.missplit) :
        return "TN"
    if mapping.covering and (not mapping.split or mapping.missplit) :
        return "FN"
    if not mapping.covering and mapping.split and not mapping.missplit :
        return "FP"

def merge_split(contig_reads) :
    candidates = []
    split_alignments = contig_reads.sort_values(by=['start_split','end_split']).reset_index(drop=True)
    while not split_alignments.empty :
        current = split_alignments.loc[0,['start_split','end_split']] 
        split_alignments['start_split']=pd.to_numeric(split_alignments['start_split'])
        split_alignments['end_split']=pd.to_numeric(split_alignments['end_split'])
        selected_reads = split_alignments.query(
            '(-5 <= start_split-{current_start} <= 5) and (-5 <= end_split-{current_end} <= 5)'.format(
                current_start=float(current.start_split),
                current_end=float(current.end_split)
                )
            )
        old_size = 1
        while len(selected_reads) != old_size :
            current.start_split = round(selected_reads['start_split'].mean())
            current.end_split = round(selected_reads['end_split'].mean())
            old_size = len(selected_reads)
            selected_reads = split_alignments.query(
                '(-5 <= start_split-{current_start} <= 5) and (-5 <= end_split-{current_end} <= 5)'.format(
                    current_start=current.start_split,
                    current_end=current.end_split
                    )
                )
        candidates.append(pd.Series(
                name= '|'.join([
                    contig_reads.name,
                    str(int(current.start_split)),
                    str(int(current.end_split))
                    ]),
                data= [contig_reads.name,current.start_split,current.end_split,len(selected_reads)],
                index=['contig','start','end','depth']
                ))
        split_alignments = split_alignments.drop(selected_reads.index).reset_index(drop=True)
    return pd.DataFrame(candidates)

def alignment_analysis(args_dict) :

    print("Analysis "+args_dict["name"])
    
    path_rmpg_pick = "/".join([
        results_dir,
        "_".join(["reads_mapping",*args_dict["name"].split()])
        ])
    path_resreads_pick = "/".join([
        results_dir,
        os.path.basename(args_dict['library'])+"res"
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
        
        
        path_ref_pick = '/'.join([
                results_dir,
                os.path.basename(args_dict['reference'])
                ])+'.gz'
        path_ctrl_pick = '/'.join([
                results_dir,
                os.path.basename(args_dict['control'])
                ])+'.gz'
        if not os.path.exists(path_ref_pick) or not os.path.exists(path_ctrl_pick):
            print("contig parsing")
            contigs = parse_fasta(args_dict['reference'])
            print("Control parsing")
            control = parse_control_introns(args_dict['control'])
            
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
            contigs.to_pickle(path_ref_pick)
            introns.to_pickle(path_ctrl_pick)
        else :
            contigs = pd.read_pickle(path_ref_pick)
            introns = pd.read_pickle(path_ctrl_pick)

        # For each intron, we determine all the reads which cover the insertion locus and the position in the read of this insertion locus
        # (precision : the function is called on intronns DataFrame but it returns a library-like DataFrame)
        if not os.path.exists(path_rmpg_pick) :
            print("library computation")
            with prl.ProcessPoolExecutor(max_workers=8) as ex :
                s_t = time.time()
                introns_split = np.array_split(introns,ex._max_workers)
                lectures = pd.concat(ex.map(prlz_process_intron,introns_split,repeat(library,ex._max_workers)))
                print(time.time()-s_t)
            lectures = library.join(lectures,lsuffix='',rsuffix='_cov').loc[:,lectures.columns]
            lectures.loc[lambda df : df.covering != True, "covering"] = False
            lectures.to_pickle(path_resreads_pick)
        else :
            lectures = pd.read_pickle(path_resreads_pick)
        
        if not os.path.exists(path_rmpg_pick) :
            
            print("alignment parsing")
            bamfile = pysam.AlignmentFile(args_dict["bamfile"], "rb")
            alignments = [{
                        'query_name' : record.query_name,
                        'reference_name' : record.reference_name,
                        'reference_start' : record.reference_start,
                        'reference_end' : record.reference_end,
                        'cigartuples' : record.cigartuples,
                        'is_secondary' : record.is_secondary,
                        'is_supplementary' : record.is_supplementary,
                        'mapping_quality' : record.mapping_quality
                        } for record in bamfile.fetch(until_eof=True)]
            
            print("alignment computation")
            
            with prl.ProcessPoolExecutor(max_workers=8) as ex :
                s_t = time.time()
                align_split = np.array_split(alignments,ex._max_workers)
                reads_mapping = pd.concat(ex.map(process_bam,
                    align_split,
                    repeat(contigs,ex._max_workers),
                    repeat(introns,ex._max_workers),
                    repeat(lectures,ex._max_workers)))
                print(time.time()-s_t)
            print("duplicates computation")
            reads_mapping['multi_aligned'] = reads_mapping['read'].duplicated(keep=False)
            print("alignements classification")
            reads_mapping['classe'] = reads_mapping.apply(class_read,axis=1)
            reads_mapping.to_pickle(path_rmpg_pick)
        else :
            reads_mapping = pd.read_pickle(path_rmpg_pick)
        
    else :
        reads_mapping = pd.read_pickle(path_rmpg_pick)
        lectures = pd.read_pickle(path_resreads_pick)
    print()
    print("####### Results  ######")
    print(args_dict["name"])
    print(len(reads_mapping))
    print(reads_mapping['classe'].value_counts())
    print("###################")
    print()
    
    path_candidates_pick = "/".join([
        results_dir,
        os.path.basename("_".join(["candidates",*args_dict["name"].split()]))
        ])
    
    if not os.path.exists(path_candidates_pick) :
        candidates = reads_mapping.loc[lambda df : df.split == True,:].groupby('contig').apply(merge_split).droplevel(0)
        candidates.to_pickle(path_candidates_pick)
    else :
        candidates = pd.read_pickle(path_candidates_pick)
        path_ctrl_pick = '/'.join([
                results_dir,
                os.path.basename(args_dict['control'])
                ])+'.gz'
        control = pd.read_pickle(path_ctrl_pick)
    
    print('Number of candidates')
    print(len(candidates))
    
    print('Number of contigs with/without split reads')
    print(reads_mapping.groupby('contig').apply(lambda df : not df['split'].any()).value_counts())
    
    print(candidates)
    print(control)
    
    return 
    
            
if __name__ == '__main__' :
    
    print("config parsing")
    config = parse_config()
    
    print("Res dir creation")
    results_dir = config["Global"]["env_dir"]
    if not os.path.exists(results_dir) :
        os.mkdir(results_dir)
        
    args_dicts = [config[analysis] for analysis in config['Global']['analysis'].split(',')]
    alignment_analysis(args_dicts[int(sys.argv[-1])])
    
