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
from intronSeekerPlot import * 

def parse_positions(fastq_pos) :
    pos = fastq_pos.lstrip("position=").split("..")
    complement = ('complement(' in pos[0])
    start = int(pos[0].lstrip("complement("))-1
    end = int(pos[1].rstrip(")"))
    return start,end,complement

def parse_library(r1,r2) :
    
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

def compute_tr_length(contig,control) :
    return contig.length -  control.loc[lambda df : df.contig == contig.name,"length" ].sum()

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
                         (df.contig+'.modif'== str(intron.contig))
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

def simulationReport(r1,r2,reference,control,threads,bamfile,splitfile) :
    # library : pandas.DataFrame where each line is a read description
    # contigs : pandas.DataFrame where each line is a contig description 
    # control : pandas.DataFrame where each line is a simulated intron description 
    if r2 :
        library = parse_library(r1.name,r2.name)
    else :
        library = parse_library(r1.name)
    
    contigs = parse_fasta(reference.name)
    control = parse_control_introns(control.name)
    
    # Add a column to contigs with the true transcript length (without any simulated introns)
    contigs["transcript_length"] = contigs.apply(
        compute_tr_length,
        axis = 1,
        control=control
        )
    
    # Add two columns to control : one with the true insertion position of the simulated intron (in term of contig length percentage)
    # and another one with the borders of the simulated introns (in term of nucleotides)
    control = control.join(
        other = control.apply(
            compute_pos_on_contig,
            axis=1,
            contigs=contigs
        )
    )
    
    plot_insertion_in_contig(control['pos_on_contig'])
    
    # Add three columns on library : one for the intron covering reads (True/False), another for the covered intron id (if True) and the last
    # for the intron insertion position in read (if True - in term of read length percentage)
    # (precision : the function is called on control DataFrame but it returns a library-like DataFrame)
    with prl.ProcessPoolExecutor(max_workers=threads) as ex :
        introns_split = np.array_split(control,ex._max_workers)
        library_cov = pd.concat(ex.map(prlz_process_intron,introns_split,repeat(library,ex._max_workers)))
    
    library = library.join(library_cov,lsuffix='',rsuffix='_cov').loc[:,library_cov.columns]
    library.loc[lambda df : df.covering != True, "covering"] = False
    
    bamfile = pysam.AlignmentFile(bamfile.name, "rb")

    alignments = pd.DataFrame(
        [pd.Series([
            record.query_name,                        # ID of the read
            record.reference_name,                    # ID of the contig where the read is mapped
            record.reference_start,                   # Start of the alignment on the contig
            record.reference_end,                     # End of the alignment on the contig1
            library.at[record.query_name,"covering"], # Bool if the read normally covers an intron (i.e. should be split)
            record.cigartuples is not None,           # Bool if the read is mapped
            not (record.cigartuples is not None and record.reference_name.rstrip(".modif") == library.at[record.query_name,"contig"]), # Bool if the read is mapped on right contig
            record.is_secondary,
            record.is_supplementary]
            ,
            index = ["read","contig","align_start","align_end",'covering','mapped',"mismapped",'second','suppl']
            )
        for record in bamfile.fetch(until_eof=True) ])

    split_aln = pd.read_csv(splitfile.name,sep='\t').rename(columns={'#start_split':'start_split'})
    # ~ print(control.head(5))
    # ~ print(contigs.head(5))
    # ~ print(library.head(5))
    # ~ print(alignments.head(5))
    # ~ print(split_aln.head(5))

if __name__ == '__main__' :
    
    import argparse 
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-1','--R1', type=argparse.FileType('r'), required=True, dest='r1')
    parser.add_argument('-2','--R2', type=argparse.FileType('r'), required=False, dest='r2')
    parser.add_argument('-r','--reference', type=argparse.FileType('r'), required=True, dest='reference')
    group_ctrl = parser.add_mutually_exclusive_group(required=True) 
    group_ctrl.add_argument('--gbs', type=argparse.FileType('r'), dest='control')
    group_ctrl.add_argument('--frs', type=argparse.FileType('r'), dest='control')
    parser.add_argument('-t','--threads', type=int, default=1, required=False, dest='threads')
    parser.add_argument('-a','--alignment', type=argparse.FileType('r'), required=True, dest='bamfile') 
    parser.add_argument('-s','--split-alignments', type=argparse.FileType('r'), required=True, dest='splitfile') 

    args = vars(parser.parse_args())
    
    simulationReport(**args)
