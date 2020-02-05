#!/usr/bin/env python3

import os
import configparser
import numpy as np
import pandas as pd
import plotly as py
import plotly.graph_objects as go



import re
import pickle
import pysam
import gzip
import time
import sys
import concurrent.futures as prl
import tempfile
#import matplotlib.pyplot as plt
#from scipy.stats import gaussian_kde
from IPython.display import HTML
from pprint import pprint
from collections import OrderedDict
from itertools import repeat
from Bio import SeqIO
from intronSeekerPlot import *

from json import JSONEncoder
import json


#step 1  full random simulation : intronSeeker fullRandomSimulation -r -o FRS/ 
#$ conda deactivate
#module load system/Miniconda3-4.7.10;
#source activate ISeeker_environment;
#cd scripts/; 
#python3 simulation2HTML.py -g ../../archives_intronSeeker/FRS/frs_modifications.gtf -f ../../archives_intronSeeker/FRS/frs_contigs-modified.fa -o HTML
#python3 simulation2HTML.py -f../../archives_intronSeeker/IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/GBS/gbs_Cele_test1_transcripts-modified.fa -g ../../archives_intronSeeker/IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/GBS/gbs_Cele_test1_transcripts-modified.gtf -o HTML

#python3 simulation2HTML.py -g /work/project/sigenae/sarah/intronSeeker/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -r /work/project/sigenae/sarah/intronSeeker/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /work/project/sigenae/sarah/intronSeeker/FRS/CAS-A/sample1/frs_sample1_contigs.fa -o HTML/
#simulate reads : config/grinder_frs_testA.cfg
#Tests 
#../../IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/FRS/
#simulationReport.py --R1 ../../IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/FRS/sr_test1_R1.fastq.gz --R2 ../../IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/FRS/sr_test1_R2.fastq.gz --reference ../../IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/FRS/frs_test1_contigs.fa --alignment ../../IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/FRS/hisat2_test1.sort.bam --split-alignments ../../IntronSeeker_1rst_tests/Tests_Janv2020/intronSeeker/FRS/srs_test1_split_alignments.txt  --frs ???? TODO 



    
#Tableau - https://plot.ly/python/table/
  
   
#return tab1
def tabInput(fasta, gtf):
    #files and paths
    tab = go.Figure(data=[go.Table(header=dict(values=['File name/path'], fill_color='paleturquoise', align='left'),
                cells=dict(values=[[gtf, fasta]], fill_color='lavender', align='left'))
                    ])
    tab.update_layout(
        autosize=False,
        width=900,
        height=300,
        paper_bgcolor='rgba(0,0,0,0)',
    )
    return tab

#return tab2
def tabGlobal(fasta, gtf):
    #table contig/transcript nb and nb tot of contig/transcript modified
    tab = go.Figure(data=[go.Table(header=dict(values=['Nb contig(s)/transcript(s)', 'Nb contig(s)/transcript(s) modified'], fill_color='paleturquoise', align='left'),
                cells=dict(values=[[count_seq_from_fa(fasta)], [nbContigsModified(gtf)]], fill_color='lavender', align='left'))
                    ])
    tab.update_layout(
        autosize=False,
        height=400,
        paper_bgcolor='rgba(0,0,0,0)',
    )
    return tab

#return tab3
def tabByContig(gtf):
    nb_features, nb_features_by_ctg, ctg_descr=stat_from_gtf(gtf)

    #https://stackoverflow.com/questions/18837262/convert-python-dict-into-a-dataframe
    df=pd.DataFrame(list(ctg_descr.items()), columns=['Features', 'nbContigsByFeature'])

    tab3 = go.Figure(data=[go.Table(
    header=dict(values=list(df.columns),
                fill_color='paleturquoise',
                align='left'),
    cells=dict(values=[df.Features, df.nbContigsByFeature],
               fill_color='lavender',
               align='left'))
    ])
    tab3.update_layout(
        autosize=False,
        height=400,
        paper_bgcolor='rgba(0,0,0,0)',
    )

    return tab3


#return int
def nbContigsModified(gtf):
    n=count_items_from_gtf(gtf)
    nb_contig_modif=0
    for ctg in n.values():
        nb_contig_modif += ctg
    #print ('nb_contig_modif', nb_contig_modif)

    return nb_contig_modif

#https://plot.ly/python/creating-and-updating-figures/



# Return int
def count_seq_from_fa(fasta):
    nb_seq = 0;
    for line in open(fasta):
        if line.startswith('>'):
            nb_seq += 1
    return nb_seq

# Return dict : item_name => int
def count_items_from_gtf(gtf):
    nb_items = dict();
    for line in open(gtf):
        if not line.startswith('#'):
            nb_items[line.split()[2]] = nb_items.get(line.split()[2], 0) + 1
    return nb_items

# Return two dict
def stat_from_gtf(gtf):
    ctg_descr          = dict()
    nb_features        = dict()
    nb_features_by_ctg = dict()
    tmp_array = []
    for line in open(gtf):
        if not line.startswith('#'):
            k = line.split()[0]
            if k in ctg_descr:
                ctg_descr[k][line.split()[2]] = ctg_descr[k].get(line.split()[2], 0) + 1
                if(line.split()[2] not in tmp_array):
                    tmp_array.append(line.split()[2])
                    nb_features_by_ctg[line.split()[2]] = nb_features_by_ctg.get(line.split()[2], 0) + 1
            else:
                ctg_descr[k] = dict()
                ctg_descr[k][line.split()[2]] = ctg_descr[k].get(line.split()[2], 0) + 1
                tmp_array = []
                tmp_array.append(line.split()[2])
                nb_features_by_ctg[line.split()[2]] = nb_features_by_ctg.get(line.split()[2], 0) + 1
                
            nb_features[line.split()[2]] = nb_features.get(line.split()[2], 0) + 1
    
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
    
    return nb_features, nb_features_by_ctg, ctg_descr
   
def simulationReport(fasta : str, gtf : str, output : str, prefix : str) :
    output_path = output + "/html";
    if prefix:
        output_path += "_" + prefix;
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.mkdir(output)
    
    print(stat_from_gtf(gtf.name))
    
    #create the report in a HTML file

    #print('count gtf items', count_items_from_gtf(gtf.name))
    #nb_features, nb_features_by_ctg, ctg_descr=stat_from_gtf(gtf.name)
    #print('nb_features', nb_features)
    #print('nb_features_by_ctg', nb_features_by_ctg)
    #print('ctg_descr', ctg_descr)
    
    tabInputList=tabInput(fasta.name, gtf.name)
    tab1=tabInputList\
        .to_html()\
        .replace('<table border="0" class="dataframe">','<table class="table table-striped">') # use bootstrap styl

    tabStatList=tabGlobal(fasta.name, gtf.name)
    tab2=tabStatList\
        .to_html()\
        .replace('<table border="0" class="dataframe">','<table class="table table-striped">') 
    
    tabByCtg=tabByContig(gtf.name)
    tab3=tabByCtg\
        .to_html()\
        .replace('<table border="0" class="dataframe">','<table class="table table-striped">') 

    reportFile = '%s/%s' % (output,'simulation_report.html')
    #to append html : a+	
    with open(reportFile,"w+") as f:
       contenu ="""
	   <html>
      <head>
	   <link rel="stylesheet" href="/work/project/sigenae/sarah/intronSeeker/scripts/bootstrap.min.css">
	   <style>body{ margin:0 100; background:whitesmoke; }</style>
	   </head>
	   <body>
	   <h1>SIMULATION REPORT</h1>
	   <h3>Input:  """ + tab1 + """ </h3>
	   <h3>Global statistiques: """ + tab2 + """</h3>
	   <h3>Statistiques by contig modified:</h3>
           <h4>Features types: """ + tab3 + """</h4>
       </body>
       </html>"""
       f.write(contenu)
       f.close() 

    return

if __name__ == '__main__' :
    
    import argparse 
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='fasta')
    parser.add_argument('-g','--gtf', type=argparse.FileType('r'), required=True, dest='gtf') 
    parser.add_argument('-o','--output', type=str, required=True, dest='output')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')

    args = vars(parser.parse_args())
    
    simulationReport(**args)
