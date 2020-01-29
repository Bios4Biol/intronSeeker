#!/usr/bin/env python3

import re
import pickle
import os
import pandas as pd
import gzip
import time
import sys
import numpy as np
from datetime import datetime
from datetime import time as dt_tm
from datetime import date as dt_date
#import plotly as py
#import plotly.plotly as py
#import plotly.tools as plotly_tools
#import plotly.graph_objs as go
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
#import matplotlib.pyplot as plt
#from scipy.stats import gaussian_kde
from IPython.display import HTML
from pprint import pprint
from collections import OrderedDict
import configparser
import concurrent.futures as prl
from itertools import repeat


#step 1  full random simulation : intronSeeker fullRandomSimulation -r -o FRS/ 
#$ conda deactivate
#test : module load system/Python-3.7.4; cd scripts/; 
#python3 simulation2HTML.py -g ../FRS/frs_modifications.gtf -r ../FRS/frs_contigs-modified.fa -f ../FRS/frs_contigs.fa -o HTML/
#python3 simulation2HTML.py -g /work/project/sigenae/sarah/intronSeeker/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -r /work/project/sigenae/sarah/intronSeeker/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /work/project/sigenae/sarah/intronSeeker/FRS/CAS-A/sample1/frs_sample1_contigs.fa -o HTML/
#simulate reads : config/grinder_frs_testA.cfg

def plot_gtf(gtf_file)
    #Tableau 1 - https://plot.ly/python/table/
    df = pd.read_fwf(gft_file, sep = "\t")

    summary_table_1 = go.Figure(data=[go.Table(
    header=dict(values=list(df.columns),
                fill_color='paleturquoise',
                align='left'),
    cells=dict(values=[df.col1, df.col2, df.col3, df.col4, df.col5, df.col6, df.col7, df.col8, df.col9],
               fill_color='lavender',
               align='left'))
    ])
    summary_table_1.show()
    #Mise en forme du tableau
    summary_table_1 = df.describe()
    summary_table_1 = summary_table_1.to_html()
    return summary_table_1 
    
def report_HTML(html_file)
    f = open(html_file'/report.html','w')
    #Report HTML
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
            <h1>2014 technology and CPG stock prices</h1>
            <h3>Summary table: 2014 stock statistics</h3>
                ''' + summary_table_1 + '''   
        </body>
    </html>'''
    f.write(html_string)
    f.close()    
    return f    

def parse_length_introns(gtf_file) :
    table = pd.read_table(gtf_file, sep='\t', header = None, names = ['contig', 'frs', 'intron', 'start', 'end', 'x1', 'x2', 'x3', 'x4'])
    table["length"] = table["end"]-table["start"]
    return table.set_index('length')

#http://www.xavierdupre.fr/app/teachpyx/helpsphinx/i_ex.html
def parse_nb_introns_by_contig(gtf_file) :
    table = pd.read_table(gtf_file, sep='\t', header = None, names = ['contig', 'frs', 'intron', 'start', 'end', 'x1', 'x2', 'x3', 'x4'])
    l=table["contig"]
    d = {}
    for x in l:
        d[x] = d.get(x, 0) + 1
    return d

def simulationReport(fa : str, modifiedfa : str, gtf : str, output : str, prefix : str) :
    output_path = output + "/html";
    if prefix:
        output_path += "_" + prefix;
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.mkdir(output)

    nContigs = sum(1 for _ in open(fa.name))
    nContigsAvecIntrons = sum(1 for _ in open(modifiedfa.name))
    print ('Nombre de total de contigs: ' , nContigsAvecIntrons)
    print ('Nombre de contigs avec introns retenus: ' , nContigsAvecIntrons-nContigs)
    nIntrons = sum(1 for _ in open(gtf.name))
    print ('Nombre total introns retenus: ' , nIntrons)
    
    #Nombre min / max / moyen d'introns retenus par contig
    nbIntronsByContig=parse_nb_introns_by_contig(gtf.name)
    print ('Nb introns par contig: ', nbIntronsByContig)
    
    #Distribution des tailles d introns
    lengthDistrib=parse_length_introns(gtf.name)
    print ('Distribution de la taille des introns: ', lengthDistrib)
    #print uniquement ligne 1 : print ('Distribution de la taille des introns: ', lengthDistrib.head(1))

    #create the report in a HTML file
    #plot_gtf(gtf.name)
    #report_HTML(output.name)
    
    pass


if __name__ == '__main__' :
    
    import argparse 
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='fa')
    parser.add_argument('-r','--reference', type=argparse.FileType('r'), required=True, dest='modifiedfa')
    parser.add_argument('-g','--gtf', type=argparse.FileType('r'), required=True, dest='gtf') 
    parser.add_argument('-o','--output', type=str, required=False, dest='output')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')


    args = vars(parser.parse_args())
    
    simulationReport(**args)
