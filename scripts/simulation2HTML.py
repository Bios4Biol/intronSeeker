#!/usr/bin/env python3

import os
import configparser
import numpy as np
import pandas as pd
import plotly as py
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.subplots as psp
import argparse
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import pickle


import re
import gzip
import time
import sys
import concurrent.futures as prl
import tempfile
from pprint import pprint
from collections import OrderedDict
from itertools import repeat
from Bio import SeqIO

from json import JSONEncoder
import json


#step 1  full random simulation : intronSeeker fullRandomSimulation -r -o FRS/ 
#$ conda deactivate
#module load system/Miniconda3-4.7.10;
#source activate ISeeker_environment;
#cd scripts/; 
#(ISeeker_environment) sigenae@genologin1 /work/project/sigenae/sarah/intronSeeker/scripts $ python3 simulation2HTML.py -m ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -1 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R2.fastq.gz -o HTML -p tests -F  --frs  ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf  -D ../../archives_intronSeeker/TESTS/
#(ISeeker_environment) sigenae@genologin1 /work/project/sigenae/sarah/intronSeeker/scripts $ python3 simulation2HTML.py -m ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o HTML -p tests -F  -1 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R2.fastq.gz
#python3 simulation2HTML.py -m /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/HTML -p TOTO -F  -1 /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R2.fastq.gz -a /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt
#python3 simulation2HTML.py -m /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML -p test1 -F  -1 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R2.fastq.gz -a /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt -r /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_ranks.txt -S /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -H /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.flagstat.txt -bha /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.bam -bhm /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.bam -bsa /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.bam -bsm /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.bam -aia /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_HISAT2/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6

def get_html_header():
    return '''
<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta name="description" content="IntronSeeker simulation report">
    <meta name="author" content="Bardou P., Maman S., Lasguignes E., Oudin F., Cabanettes F., Klopp C.">

    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    <!-- Icons -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/open-iconic/1.1.1/font/css/open-iconic-bootstrap.min.css">
    <!-- Datatable -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.16/css/dataTables.bootstrap4.min.css">

    <!-- Personalised css -->
    <style type="text/css">
      body { font-size: .875rem; }
      .sidebar {
        position: fixed;
        top: 0;
        bottom: 0;
        left: 0;
        z-index: 100; /* Behind the navbar */
        padding: 0;
        box-shadow: inset -1px 0 0 rgba(0, 0, 0, .1);
      }
      .sidebar-sticky {
        position: -webkit-sticky;
        position: sticky;
        top: 48px; /* Height of navbar */
        height: calc(100vh - 48px);
        padding-top: 4rem;
        overflow-x: hidden;
        overflow-y: auto; /* Scrollable contents if viewport is shorter than content. */
      }
      .sidebar .nav-link {
        margin-right: 4px;
        color: #212529;
      }
      .sidebar .nav-link:hover,
      .sidebar .nav-link.active { color: #999; }
      .sidebar-heading {
        font-size: .75rem;
        text-transform: uppercase;
      }
      .navbar-brand {
        padding-top: .75rem;
        padding-bottom: .75rem;
        font-size: 1rem;
        background-color: rgba(0, 0, 0, .7);
        box-shadow: inset -1px 0 0 rgba(0, 0, 0, .7);
				z-index: 200;
      }
      .navbar .form-control {
        padding: .75rem 1rem;
        border-width: 0;
        border-radius: 0;
				z-index: 200;
      }
      .form-control-dark {
        color: #fff;
        background-color: rgba(255, 255, 255, .1);
        border-color: rgba(255, 255, 255, .1);
      }
      .form-control-dark:focus {
        border-color: transparent;
        box-shadow: 0 0 0 3px rgba(255, 255, 255, .25);
      }
      .border-top { border-top: 1px solid #e5e5e5; }
      .border-bottom { border-bottom: 1px solid #e5e5e5; }
	  .valn { vertical-align: middle !important; }
	  .anchor{
  			display: block;
  			height: 83px; /*same height as header*/
  			margin-top: -83px; /*same height as header*/
  			visibility: hidden;
	  }
    </style>

    <!-- Plotly -->
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script> 

    <title>IntronSeeker simulation report</title>
  </head>
'''

def get_html_body1(flagstat=""):
    r = '''
  <body>
    <nav class="navbar navbar-dark sticky-top bg-dark flex-md-nowrap p-0">
      <a class="navbar-brand col-sm-3 col-md-2 mr-0" href="#">IntronSeeker simulation report</a>
    </nav>

    <div class="container-fluid">
      <div class="row">
        <nav class="col-md-2 d-none d-md-block bg-light sidebar">
          <div class="sidebar-sticky">
            <ul class="nav flex-column">
              <li class="nav-item">
                <a class="nav-link" href="#inputs-parameters">
                  <span class="oi oi-file" aria-hidden="true"></span>
                  Input files
                </a>
              </li>
			  <li class="nav-item">
				<a class="nav-link" href="#ref-descr">
				  <span class="oi oi-collapse-down" aria-hidden="true"></span>
					Reference
				</a>
			  </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#gstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Global statistics
			    	</a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#nb_ctg_by_feature">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Nb. of seq. by feature type
			    	</a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#ctg_descr">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Nb. of seq. with same features
			    	</a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#feat_len_dist">
				        <span class="oi oi-graph" aria-hidden="true"></span>
				        Features len. distribution
				    </a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#feat_dist_intron">
				        <span class="oi oi-graph" aria-hidden="true"></span>
				        Introns positions
				    </a>
			    </li>  
			  <li class="nav-item">
				<a class="nav-link" href="#read-descr">
				  <span class="oi oi-collapse-up" aria-hidden="true"></span>
					Reads
				</a>
			  </li>
			  <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#readgstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Global statistics
			    	</a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#abundstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	FRS abundance
			    	</a>
			    </li> '''           
    if flagstat:
        r += '''
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#readastat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Alignment statistics
			    	</a>
			    </li>'''      
    r += '''                    
              <li class="nav-item">
				<a class="nav-link" href="#flag-descr">
				  <span class="oi oi-collapse-up" aria-hidden="true"></span>
					Mapping
				</a>
			  </li>
			  <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#flaghstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Hisat2 statistics
			    	</a>
			    </li>  
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#flagsstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Star statistics
			    	</a>
			    </li> 
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#split">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Split detection
			    	</a>
			    </li>    
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#assemblystat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Intron insertion
			    	</a>
			    </li>   
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#candidatstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Candidats statistics
			    	</a>
			    </li> 
            </ul>
          </div>
          <div style="text-align:center;font-size:smaller;color:darkgrey;margin-top:-25px">
		    Produced by IntronSeeker_v1.0<br>
		    Copyright © 2020, <img style="width:18px;padding-bottom:2px" src="https://www.inrae.fr/themes/custom/inrae_socle/favicon.ico"><!--<img src="http://www.inra.fr/extension/itkinra/design/inra/images/favicon.ico">-->
		    <a style="color:#212529;" href="https://inrae.fr" target="_blank">INRAE</a><br>
		    Designed by the <a style="color:#212529;" href="http://sigenae.org" target="_blank">Sigenae</a> team.
          </div>
        </nav>
		<main role="main" class="col-md-9 ml-sm-auto col-lg-10 pt-3 px-4">
'''
    return r
    

def get_html_inputfiles(fasta:str, mfasta:str, gtf:str, r1:str, r2=""):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pb-2 border-bottom">
			<h1 class="h4">Input files</h1>
			<span class="anchor" id="inputs-parameters"></span>
		</div>
		<div class="d-flex">
			<div class="mt-4 mr-4 pl-0 col-md-6">
				<h5>Files</h5>
				<ul class="list-group">
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">FASTA file <span class="badge badge-success badge-pill ml-4">'''+fasta+'''</span></li>
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Modified FASTA file <span class="badge badge-success badge-pill ml-4">'''+mfasta+'''</span></li>
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">GTF file of modified sequences <span class="badge badge-success badge-pill ml-4">'''+gtf+'''</span></li>
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Read1 FASTQ file <span class="badge badge-success badge-pill ml-4">'''+r1+'''</span></li>'''
    if r2:
        r += '''
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Read2 FASTQ file <span class="badge badge-success badge-pill ml-4">'''+r2+'''</span></li>'''
    r += '''
				</ul>
			</div>
		</div>
'''
    return r


def get_html_seq_descr(global_stat:dict, nb_ctg_by_feature:dict, ctg_descr:dict, gtf:str, pos:dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Sequences</h1>
                <span class="anchor" id="ref-descr"></span>
        </div>
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Global statistics</h5>
                <span class="anchor" id="gstat"></span>
'''+dict_to_table(global_stat,7,True)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Number of sequences by feature type</h5>
                <span class="anchor" id="nb_ctg_by_feature"></span>
'''+dict_to_table(nb_ctg_by_feature,-1,False)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Number of sequences with same feature(s)</h5>
                <span class="anchor" id="nb_ctg_by_feature"></span>
'''+dict_to_table(ctg_descr,-1,False)+'''
            </div>
        </div>
'''

    # Len dist for gtf
    len_by_features, feature_names = len_dist_from_gtf(gtf)
    r += '''
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-6">
                <h5>Features length distribution</h5>
                <span class="anchor" id="feat_len_dist"></span>
'''+plot_dist(len_by_features, feature_names)+'''
            </div>
'''
    r += '''
            <div class="mt-4 mr-0 pl-0 col-md-6">
                <h5>Distribution of introns insertion position along the contigs</h5>
                <span class="anchor" id="feat_dist_intron"></span>
'''+plot_insertion_in_contig(pos)+'''
            </div>
        </div>  
'''
    return r

def get_html_ranks_descr(rank:dict, real_abund_perc:dict):
    r = '''
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <h5>Percentage of abundance of each contig in FRS library</h5>
                <span class="anchor" id="abundstat"></span>
'''+abondance_model(rank, real_abund_perc)+'''
            </div>
        </div>    
'''
    return r


def get_html_mapping_descr(names:dict, mapping_hisat_all:dict, mapping_hisat_mixed:dict, mapping_star_all:dict, mapping_star_mixed:dict, colors:dict, library:dict):
    r = '''
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <h5> Counting table and barplots of mapped covering reads' main characteristics</h5>
                <span class="anchor" id="abundstat"></span>
'''+ plot_covering_reads(*zip(names,[
                    mapping_hisat_all.loc[lambda df : df.covering == True,:],
                    mapping_hisat_mixed.loc[lambda df : df.covering == True,:],
                    mapping_star_all.loc[lambda df : df.covering == True,:],
                    mapping_star_mixed.loc[lambda df : df.covering == True,:]
                    ]),
                colors=colors,
                library=library) +'''
            </div>
        </div>    
'''
    return r


def get_html_reads_descr(global_stat_fastq : dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Reads</h1>
                <span class="anchor" id="read-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
                <span class="anchor" id="readgstat"></span>
'''+dict_to_table(global_stat_fastq,-1,True)+'''
            </div>
        </div>
'''
    return r
    

#def get_html_bam_descr(global_stat_flagstat : dict):
#    r = '''
#        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
#            <h1 class="h4">Bam (remove ?)</h1>
#                <span class="anchor" id="read-descr"></span>
#        </div>
#		<div class="d-flex">
#            <div class="mt-4 mr-0 pl-0 col-md-4">
#                <span class="anchor" id="readastat"></span>
#'''+dict_to_table(global_stat_flagstat,-1,True)+'''
#            </div>
#        </div>
#'''
#    return r    


def get_html_flagstat_descr(global_stat_flagstat_hisat2:dict, global_stat_flagstat_star:dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Mapping</h1>
                <span class="anchor" id="flag-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>HiSAT2</h5>
                <span class="anchor" id="flaghstat"></span>
'''+dict_to_table(global_stat_flagstat_hisat2,-1,True)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>STAR</h5>
                <span class="anchor" id="flagsstat"></span>
'''+dict_to_table(global_stat_flagstat_star,-1,True)+'''
            </div>
        </div>
'''
    return r

#def get_html_flagstat_descr(global_stat_flagstat_hisat2:dict, global_stat_flagstat_star:dict):
#    r = '''
#		<div class="d-flex">
#            <div class="mt-4 mr-0 pl-0 col-md-4">
#                <h5>HiSAT2</h5>
#                <span class="anchor" id="flaghstat"></span>
#'''+dict2_to_table(global_stat_flagstat_hisat2,global_stat_flagstat_star,-1,True)+'''
#            </div>
#        </div>
#'''
#    return r    

def get_html_assemblathon_descr(global_stat_assemblathon:dict):
    r = '''
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
            <h5>Pseudo-assembly comparison with assemblathon statistics</h5>
                <span class="anchor" id="assemblystat"></span>
'''+dict_to_table(global_stat_assemblathon,-1,True)+'''
            </div>
        </div>
'''
    return r

def get_html_table_descr(global_stats_table):
    r = '''
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
            <h5>TO COMPLETE</h5>
                <span class="anchor" id="assemblystat"></span>
'''+dict_to_table(global_stats_table,-1,True)+'''
            </div>
        </div>
'''
    return r


def get_html_split(mapping_hisat_all:dict , mapping_hisat_mixed:dict, mapping_star_all:dict, names:dict, mapping_star_mixed:dict, colors:dict):
    r = '''
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <h5>Barplots of split reads signal detection</h5>
                <span class="anchor" id="split"></span>
'''+plot_class_reads(*zip(names,[
        mapping_hisat_all,
        mapping_hisat_mixed,
        mapping_star_all,
        mapping_star_mixed
        ]),
        colors=colors)+'''
            </div>
        </div>    
'''
    return r


    
def get_html_candidat_descr(global_stat_candidat:dict):
    r = '''
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Candidats</h5>
                <span class="anchor" id="candidatstat"></span>
'''+dict_to_table(global_stat_candidat,-1,True)+'''
            </div>
        </div>  
'''

    return r      
    
def get_html_footer():
    return '''
        </main>
      </div>
    </div>

    <!-- jQuery first, then Popper.js, then Bootstrap JS -->
    <script type="text/javascript" language="javascript" src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
    <script type="text/javascript" language="javascript" src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
    <script type="text/javascript" language="javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>

    <!-- Datatable -->
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/dataTables.bootstrap4.min.js"></script>
  </body>
</html>
'''


# Return HTML table from dict
# Param 1 : dict (key, val)
# Param 2 : int for which line number first col (key) will be right align
# Param 3 : bool to remove first char of the key (used to sort by key)
def dict_to_table(d : dict, i : int, rmfirstchar : bool):
    table = '''
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
'''
    c = 0
    for k, v in sorted(d.items(), key=lambda t: t[0]):
        table += "<tr><td class='valn"
        if(c >= i & i!=-1):
            table += " text-right"
        if(rmfirstchar):
            table += "'>" + k[1:] + "</td><td class='valn text-right'>" + str(v) + "</td></tr>"
        else:
            table += "'>" + k + "</td><td class='valn text-right'>" + str(v) + "</td></tr>"
        c += 1
    table += "</tbody></table>"
    return table

# Return HTML table from dict
# Param 1 : dict (key, val)
# Param 2 : int for which line number first col (key) will be right align
# Param 3 : bool to remove first char of the key (used to sort by key)
def dict2_to_table(d1 : dict, d2:dict, i : int, rmfirstchar : bool):
    table = '''
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
'''
    c = 0
    c2 = 0
    for k, v in sorted(d1.items(), key=lambda t: t[0]):
        table += "<tr><td class='valn"
        if(c >= i & i!=-1):
            table += " text-right"
        if(rmfirstchar):
            table += "'>" + k[1:] + "</td><td class='valn text-right'>" + str(v) + "</td></tr>"
        else:
            table += "'>" + k + "</td><td class='valn text-right'>" + str(v) + "</td></tr>"
        c += 1
    for k2, v2 in sorted(d2.items(), key=lambda t: t[0]):
        table += "<tr><td class='valn"
        if(c2 >= i & i!=-1):
            table += " text-right"
        if(rmfirstchar):
            table += "'>" + k2[1:] + "</td><td class='valn text-right'>" + str(v2) + "</td></tr>"
        else:
            table += "'>" + k2 + "</td><td class='valn text-right'>" + str(v2) + "</td></tr>"
        c2 += 1    
    table += "</tbody></table>"
    return table    


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


# Distribution plot
def plot_dist(len_by_features, feature_names):
    hist_data = len_by_features
    group_labels = feature_names
    colors = ['#333F44', '#37AA9C', '#94F3E4']
    # Create distplot with curve_type set to 'normal'
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=True, colors=colors)
    fig.update_layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')


# Plot introns position on a contigs
def plot_insertion_in_contig(positions) :
    hist = go.Histogram(
            x=positions,
            xbins=dict(
                start=0,
                end=100,
                size=2),
            marker=dict(
                color='purple'
            )
    )
    layout = go.Layout(xaxis=dict(
                           title="% of contig length"),
                       yaxis=dict(
                           title="Count"))
    fig = go.Figure(data=[hist],layout=layout)
    fig.update_layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

# Plot ranks file
def abondance_model(rank:dict, real_abund_perc:dict) :
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x = rank,
            y = real_abund_perc,
            mode = 'lines',
            name = 'Simulated abundance model'
        )
    )

    fig.add_trace(
        go.Scatter(
            x = rank,
            y = real_abund_perc,
            mode = 'lines',
            name ='Waited abundance model'
        )
    )

    fig.update_layout(
        xaxis=dict(title="Contigs"),
        yaxis=dict(title="Relative Abundance percentage",
            range=[-0.25,0.5])
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

# Barplots of split reads signal detection for HiSAT2
def plot_class_reads(*args,**kwargs) :
    series = [[],[]]
    colors = kwargs['colors']
    
    fig = psp.make_subplots(
        rows = 2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.03,
        specs=[[{"type":"table"}],
               [{"type":"bar"}]])

    for name,val in args :
        for idx,val in enumerate([val,val.loc[lambda df : df.contig == df.contig.str.rstrip(".ori"),:]]) :
            s = val.loc[:,'classe'].value_counts()
            s.name=name
            
        
            if idx == 0 :
                visibility = True
            else :
                visibility = False
        
            fig.add_trace(
                go.Bar(
                    visible=visibility,
                    x=s.index[1:],
                    y=s.values[1:]/s.sum()*100,
                    name = s.name,
                    marker_color=colors[s.name]
                ),row=2,col=1)
            
            s['Total'] = s.sum()
            series[idx].append(s)

    for idx,serie in enumerate(series) :
        table = pd.concat(serie,axis=1,sort=False).fillna(0)
        
        
        if idx == 0 :
            visibility = True
        else :
            visibility = False
            
        fig.add_trace(
            go.Table(
                visible=visibility,
                columnwidth=[40,80],
                header = dict(
                    values=["",*[c.join(["<b>","</b>"]) for c in table.columns.to_list()]],
                    fill_color='lightgrey'
                    ),
                cells = dict(
                    values=[[v.join(["<b>","</b>"]) for v in list(table.reset_index().T.values)[0]],
                           *list(table.reset_index().T.values)[1:]],
                    fill=dict(color=['lightgrey','lavender'])
                    )
                ),
            row=1,col=1)
    
    fig.update_layout(
        height=700,
        legend=dict(
            x=-0.4,y=0.40
            ),
        updatemenus=[
            go.layout.Updatemenu(
                buttons=[
                    dict(
                        args=[{"visible":[True,False]}],
                        label='All alignements',
                        method='restyle'
                        ),
                    dict(
                        args=[{"visible":[False,True]}],
                        label='Alignements on contigs with introns',
                        method='restyle'
                        ),
                    ]
                )
            ]
        )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

# Plot : Counting table and barplots of mapped covering reads' main characteristics.
#def plot_covering_reads(*args, **kwargs) :
def plot_covering_reads(*args, library:dict, colors:dict) :
    series=[]
    fig = go.Figure()
    for name, val in args :
        s = pd.Series(name=name)
        s['Covering']=len(val)
        s['Unmapped']=len(val.loc[lambda df : df.mapped == False])
        tmp = val.merge(library,left_on='read',right_index=True,suffixes=("","_lib"))
        s['Mismapped'] = len(tmp.loc[lambda df : (df.contig.str.rstrip('.ori') != df.contig_lib.str.rstrip('.ori'))& (df.mapped == True)])
        s['Unsplit'] = len(val.loc[lambda df : (df.mapped==True)&(df.split==False)])
        s['Missplit'] = len(val.loc[lambda df : (df.split==True)&(df.missplit==True)])
        s['Correct splitting'] = len(val.loc[lambda df : df.classe == 'TP'])
        series.append(s)
        
        to_plot = s[['Unmapped','Unsplit','Correct splitting']]/s['Covering']*100
        fig.add_trace(
            go.Bar(
                    x=to_plot.index,
                    y=to_plot.values,
                    name = s.name,
                    marker_color=colors[s.name]
            ))
    table = pd.concat(series,axis=1,sort=False)
        
    fig.update_layout(
        title='Global mapping results on introns-covering reads',
        xaxis=dict(title='Lectures charecteristics'),
        yaxis=dict(title='Percentage of total covering reads alignements')
        )

    return py.offline.plot(fig, include_plotlyjs=False, output_type='div'), table

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
    
# Return panda which contains gtf features desc (seqref feature start end)  -- WARNING : header=0 car il y a déjà une lettre de titre commentée !!
def parse_candidat(candidat) :
    t = pd.read_table(candidat, usecols=[0,2,3,4,6], names=['ID', 'start', 'end', 'depth', 'filter'],  header=0)
    print(type(t))
    print(t.dtypes)
    return t.set_index('ID')
    

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
        ranks = [ligne.lstrip("# ").split("\t") for ligne in rf.read().rstrip().split("\n")]
    return pd.DataFrame(data = ranks[1:], columns=ranks[0] )

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

# Parse Alignment BAM files
    #bamfile = pysam.AlignmentFile(args_dict["bamfile"], "rb")
def parse_BAM(BamPath:str):    
    bamfile = pysam.AlignmentFile(BamPath, "rb")
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
    return alignments

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


def process_bam(alignments, df_mfasta, df_features, df_library):
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
            df_library.at[record['query_name'],"covering"], # Bool if the read normally covers an intron (i.e. should be split)
            record['cigartuples'] is not None,           # Bool if the read is mapped
            not record['reference_name'] == df_library.at[record['query_name'],"contig"].rstrip(".ori"), # Bool if the read is mapped on right contig
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
                str(df_mfasta.at[record['reference_name'],"sequence"])
            ))
            introns_to_check = df_features.loc[lambda df : df.contig == record['reference_name'],:]
            for limits in zip(introns_to_check["start"],introns_to_check["end"]) :
                row["missplit"] = not (limits[0]-row.start_split in range(-3,4) and limits[1]-row.end_split in range(-3,4))
        else :
            row = begin.append(pd.Series([None,None,None,None,None],
                                         index=["start_split","end_split","split_length","split_flanks","missplit"]))
        rows.append(row)
    
    return pd.DataFrame(rows)

def compute_pos_on_read(cov_lect,intron_start):
    if not cov_lect.complement :
        return (intron_start - cov_lect.start)/(cov_lect.end-cov_lect.start)*100
    else :
        return (cov_lect.end - intron_start)/(cov_lect.end-cov_lect.start)*100

# Return df_cov_lect, a new Dataframe, which contains df_library (lecture, contig, start, end, complement) join with 3 news columns:
# "covering" : one for the intron covering reads (True/False), 
# "intron (name)" : another for the covered intron id (if True) 
# pos_on_read : intron insertion position in read (if True - in term of read length percentage)
def process_intron(intron,lectures) :
    df_cov_lect = pd.DataFrame(lectures.loc[lambda df : 
                         (df.contig+'.modif'== str(intron.contig))
                         & (intron.start > df.start)
                         & (intron.start < df.end)
                          ])              
    df_cov_lect['covering'] = True
    df_cov_lect['intron'] = intron.name
    #print('df_cov_lect[intron]',df_cov_lect['intron'])
    df_cov_lect['pos_on_read'] = df_cov_lect.apply(
            compute_pos_on_read,
            axis=1,
            intron_start=intron.start
            )
    #print('df_cov_lect[pos_on_read]',df_cov_lect['pos_on_read'])        
    '''print('df_cov_lect', df_cov_lect)
    df_cov_lect                contig  start  end  complement  covering                     intron  pos_on_read
    lecture                                                                                        
    203770/1  SEQUENCE685    457  558       False      True  SEQUENCE685.modif|556|897    98.019802'''
    return df_cov_lect

# Return DataFrame Reads
def prlz_process_intron(df_introns,library) :
    df_reads = pd.concat(
        df_introns.apply(
            process_intron,
            axis=1,
            lectures=library
            ).values
        )
    return df_reads   

############
# SUB MAIN #
############
def simulationReport(fasta:str, mfasta:str, gtf:str, r1:str, r2:str, ranksfile:str, allwiintronsassemblathon:str, flagstat:str, flgSTARintrons:str, flgHISAT2introns:str, bamHISATall:str, bamHISATmix: str, bamSTARall: str, bamSTARmix: str, candidat:str, output:str, prefix:str, force:bool, threads:int) :
    output_path = output + "/report"
    if prefix:
        output_path += "_" + prefix

    # Create output dir if not exist
    if not os.path.exists(output) :
        os.mkdir(output)
    
    # Output path filename report html
    output_file = output_path + "_simulation.html"
    if not force:
        try :
            if os.path.exists(output_file):
                raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file already exists.\n')
            exit(1)

	### MEMO 
	# fasta  = sequences used to generate reads
	# mfasta = sequences used to align reads
	
	### Build pandas ###
    # df_fasta    : pandas.DataFrame where each line is a seq description from FASTA file
    # df_mfasta   : pandas.DataFrame where each line is a modified seq description from modified FASTA file
    # df_library  : pandas.DataFrame where each line is a read description from R1 (& R2) fastq file(s)
    # df_features : pandas.DataFrame where each line is a simulated features description 
    # df_candidat : pandas.DataFrame where each line is a candidat description
    df_fasta  = parse_fasta(fasta.name, False)
    df_mfasta = parse_fasta(mfasta.name, True)

    if r2 :
        df_library = parse_library(r1.name, r2.name)
    else :
        df_library = parse_library(r1.name)

    df_features = parse_gtf(gtf.name)
    
    # Add a column to df_fasta with the "fasta" length (without any simulated features)
    df_mfasta["short_length"] = df_mfasta.apply(
        compute_tr_length,
        axis = 1,
        df_features=df_features
    )
    
    # Add two columns to df_features:
    #  1- the true insertion position of the simulated feature (in term of mfasta length percentage)
    #  2- the borders of the simulated features (in term of nucleotides)
    df_features = df_features.join(
        other = df_features.apply(
            compute_pos_on_mfasta,
            axis=1,
            df_mfasta=df_mfasta
        )
    )
    
    print('fasta head :' , df_fasta.head(5))
    print('mfasta head :' , df_mfasta.head(5))
    print('library head :' , df_library.head(5))
    print('features head :', df_features.head(5))

    # HEADER
    html = get_html_header()
    if flagstat:
        html += get_html_body1(flagstat.name)
    else:
        html += get_html_body1()

    '''elif candidat:
        html += get_html_body1(candidat.name)
    elif ranksfile:
        html += get_html_body1(ranksfile.name)    
    elif candidat and flagstat and ranksfile:
        html += get_html_body1(flagstat.name, candidat.name, ranksfile.name)'''    

    # INPUT FILES
    if r2:
        html += get_html_inputfiles(os.path.basename(fasta.name), os.path.basename(mfasta.name), os.path.basename(gtf.name), os.path.basename(r1.name), os.path.basename(r2.name))
    else:
        html += get_html_inputfiles(os.path.basename(fasta.name), os.path.basename(mfasta.name), os.path.basename(gtf.name), os.path.basename(r1.name))

    # SEQUENCE STAT
    # Global stat
    nb_distinct_features, nb_ctg_by_feature, ctg_descr = stat_from_gtf(gtf.name)
    global_stat = dict()
    global_stat["0FASTA file - Number of seq."]            = df_fasta.shape[0]
    global_stat["1FASTA file - Mean seq. length"]          = int(df_fasta['length'].mean())
    global_stat["2Modified FASTA file - Number of seq."]   = df_mfasta.shape[0]
    global_stat["3Modified FASTA file - Mean seq. length"] = int(df_mfasta['length'].mean())
    global_stat["4Number of modified sequences"]           = df_mfasta.loc[df_mfasta.index.str.contains(".modif")].shape[0]
    global_stat["5Number of distinct features in GTF"]     = df_features.feature.value_counts().shape[0]
    global_stat["6Number of features in GTF"]              = df_features.shape[0]
    c = 7
    for k, v in (df_features.feature.value_counts()).items() :
        global_stat[str(c)+k] = v
        c+=1
    
    html += get_html_seq_descr(global_stat, nb_ctg_by_feature, ctg_descr, gtf.name, df_features['pos_on_contig'])
    
    ranks = parse_rank_file(ranksfile.name)
    real = pd.DataFrame((df_library.groupby('contig').size()/len(df_library))*100,columns = ['real_abund_perc']).reset_index()
    df_abund = ranks.merge(real,right_on = 'contig',left_on='seq_id',suffixes = ('_grinder','_real'))
    print('abund',df_abund)
    html += get_html_ranks_descr(df_abund['rank'], df_abund['real_abund_perc'])

    #https://stackoverflow.com/questions/45759966/counting-unique-values-in-a-column-in-pandas-dataframe-like-in-qlik
    print('6Number of features in GTF:', df_features['contig'].nunique())
    #https://stackoverflow.com/questions/45759966/counting-unique-values-in-a-column-in-pandas-dataframe-like-in-qlik
    print('5Number of distinct features in GTF:', df_features['feature'].nunique())
    #ou
    df2=df_features.groupby('feature')['contig'].nunique()
    print('**Number of distinct features in GTF:', df2.count())
    #Replace nb_ctg_by_feature ?   # Number of ctg by feature from all GTF lines (Ex: "Exon" see in X ctg, "Intron" see in Y ctg, ...)
    #https://stackoverflow.com/questions/38309729/count-unique-values-with-pandas-per-groups/38309823
    #https://www.shanelynn.ie/summarising-aggregation-and-grouping-data-in-python-pandas/
    #print('**TEST**Number of sequences by feature type (1):', df_features.groupby('feature')['contig'].count())  #series
    print('**TEST**Number of sequences by feature type (2):', df_features.groupby('feature')[['contig']].nunique())  #panda dataframe
    print('**TEST**Number of sequences by feature type (2bis):', df_features.groupby('contig')[['feature']].nunique())  
    #print('**TEST**Number of sequences by feature type (4):', df_features.groupby('feature', as_index=False).agg({"contig": "count"})) 
    
    #df_features.insert(1, 'nbfeature', 1)
    #print ('new df_features', df_features)
    #pddf=df_features.groupby(['feature','contig'])['nbfeature'].nunique()
    #print('**TEST**Number of sequences by feature type (5):', pddf)
    #df_features = df_features.drop(columns='nbfeature')
    print('**Number of sequences by feature type: (6)', (df_features.groupby('feature')['contig'].nunique()).values) # type ; pandas.core.series.Series
    print('**Number of sequences by feature type (7):', (df_features.groupby('feature')['contig'].nunique()).index) 
    '''
    #Methode troooop longue 
    i=0
    k=0
    c=""
    f=""
    while i < df_features.shape[0]:
        if c != df_features['contig'][i]:
            k += 1
            #print('Same contig : ', str(k) + " " + df_features['contig'][i]+ ":" + str(k) + " " + df_features['feature'][i] +"(s)")
            k=0
        else:
            i+=1
            k+=1
            if f == df_features['feature'][i]:
                k+=1
                #print('same contig, same feature : ', str(k) + " " + df_features['contig'][i]+ ":" + str(k) + " " + f +"(s)") 
            else:
                k+=1
                #print('same contig, feature diff : ',str(k) + " " + df_features['contig'][i]+ ":" + str(k) + " " + f +"(s)")
        c=df_features['contig'][i] 
        f=df_features['feature'][i]  
        i += 1
        k=0
        
    #ctg_descr            = dict()   # Number of features profiles by ctg (Ex: "1 Exon & 2 Intron" see in X ctg, "3 Introns" see in Y ctg, ...)
'''
    

    # READS STAT
    # Global stat
    global_stat_fastq = dict()
    global_stat_fastq["0Number of fragments"] = df_library['contig'].count()
    global_stat_fastq["1Mean coverage"] = 0
    for i,row in df_library.iterrows():
        global_stat_fastq["1Mean coverage"] += row['end'] - row['start'] + 1
    global_stat_fastq["1Mean coverage"] /= (global_stat["1FASTA file - Mean seq. length"] * global_stat["0FASTA file - Number of seq."])
    global_stat_fastq["1Mean coverage"] = round(global_stat_fastq["1Mean coverage"], 2)
    #to remove ?
    global_stat_fastq["2Min fragments length"] = df_features['length'].min()
    global_stat_fastq["3Max fragments length"] = df_features['length'].max()
    global_stat_fastq["4Mean fragments length"] = round(df_features['length'].mean())
        
    html += get_html_reads_descr(global_stat_fastq)
    
    ## ALIGNMENT STAT
    '''if flagstat:
        nbreads, mapped, paired, proper = parse_flagstat2(flagstat.name)
        global2_stat_flagstat = dict()
        global2_stat_flagstat["0Total number of reads passing quality controls"] = nbreads
        global2_stat_flagstat["1Mapped"]                                         = mapped
        global2_stat_flagstat["2Paired in sequencing"]                           = paired
        global2_stat_flagstat["3Properly paired"]                                = proper
        html += get_html_bam_descr(global2_stat_flagstat)'''

    ## ALIGNMENT STATS
    df_flag_all_hisat = parse_flagstat(flgHISAT2introns.name, len(df_library),"All with introns - Hisat2")
    #df_flag_half_hisat = parse_flagstat(flgSTAR, len(df_library),"Mix-states contigs - Hisat2")
    df_flag_all_star = parse_flagstat(flgSTARintrons.name, len(df_library),"All with introns - STAR")
    #df_flag_half_star = parse_flagstat(flgHISAT2, len(df_library),"Mix-states contigs - STAR")    

    print(pd.concat([df_flag_all_hisat,df_flag_all_star],axis=1,sort=False).fillna(0))
    df_flag_all=pd.concat([df_flag_all_hisat,df_flag_all_star],axis=1,sort=False).fillna(0)
    global_stat_flagstat_hisat2 = dict()
    print(df_flag_all.shape[0])#lignes 
    print(df_flag_all.shape[1])#colonnes 
    print(df_flag_all.iloc[0,1])
    global_stat_flagstat_hisat2["0Total counts of reads to map"] = str(len(df_library))
    global_stat_flagstat_hisat2["1Total count"]     = df_flag_all.iloc[0,0]
    global_stat_flagstat_hisat2["2Secondary"]       = df_flag_all.iloc[1,0]
    global_stat_flagstat_hisat2["3Mapped"]          = df_flag_all.iloc[2,0]
    global_stat_flagstat_hisat2["4Properly paired"] = df_flag_all.iloc[3,0]
    global_stat_flagstat_hisat2["5Singletons"]      = df_flag_all.iloc[4,0]
    global_stat_flagstat_star= dict()
    global_stat_flagstat_star["0Total counts of reads to map"] = str(len(df_library))
    global_stat_flagstat_star["1Total count"]     = df_flag_all.iloc[0,1]
    global_stat_flagstat_star["2Secondary"]       = df_flag_all.iloc[1,1]
    global_stat_flagstat_star["3Mapped"]          = df_flag_all.iloc[2,1]
    global_stat_flagstat_star["4Properly paired"] = df_flag_all.iloc[3,1]
    global_stat_flagstat_star["5Singletons"]      = df_flag_all.iloc[4,1]
    html += get_html_flagstat_descr(global_stat_flagstat_hisat2, global_stat_flagstat_star)
    

    #MAPPING
    #Compare assemblathon files
    #df_assemblathon_wo = parse_assemblathon(wointronsassemblathon.name,"FRS without introns")
    df_assemblathon_all = parse_assemblathon(allwiintronsassemblathon.name, "FRS all modified")
    #df_assemblathon_half = parse_assemblathon(halfwiintronsassemblathon.name,"FRS Mixed states")
    #df_assemblathon=pd.concat([df_assemblathon_wo,df_assemblathon_all.name,df_assemblathon_half],axis=1)
    print('df_assemblathon_all', df_assemblathon_all)
    global_stat_assemblathon = dict()
    print(df_assemblathon_all.shape[0])#lignes 
    print(df_assemblathon_all.shape[1])#colonnes 
    global_stat_assemblathon["0Number of contigs"]         = df_assemblathon_all.iloc[0,0]
    global_stat_assemblathon["1Total size of contigs"]     = df_assemblathon_all.iloc[1,0]
    global_stat_assemblathon["2Longest contig"]            = df_assemblathon_all.iloc[2,0]
    global_stat_assemblathon["3Shortest contiged"]         = df_assemblathon_all.iloc[3,0]
    global_stat_assemblathon["4Number of contigs > 1K nt"] = df_assemblathon_all.iloc[4,0]
    global_stat_assemblathon["5N50 contig length"]         = df_assemblathon_all.iloc[5,0]
    global_stat_assemblathon["6L50 contig count"]          = df_assemblathon_all.iloc[6,0]
    html += get_html_assemblathon_descr(global_stat_assemblathon)
    '''
    #Analysis of alignments by seaching split reads and comparing with simulated introns   
    #Split read signal analysis
    #Effectives table and barplots of split reads signal detection for STAR and HiSAT2 for the two types of reference.
    # # Add three columns on df_library : one for the intron covering reads (True/False), another for the covered intron id (if True) and the last
    # for the intron insertion position in read (if True - in term of read length percentage)
    # (precision : the function is called on df_library DataFrame but it returns a library-like DataFrame)
    with prl.ProcessPoolExecutor(max_workers=threads) as ex :
        introns_split = np.array_split(df_features,ex._max_workers)
        library_cov = pd.concat(ex.map(prlz_process_intron,introns_split,repeat(df_library,ex._max_workers)))
    
    df_library = df_library.join(library_cov,lsuffix='',rsuffix='_cov').loc[:,library_cov.columns]
    df_library.loc[lambda df : df.covering != True, "covering"] = False
    #process_bam(alignments, contigs, introns, library)
    mapping_hisat_all   = pd.read_pickle(process_bam(parse_BAM(bamHISATall.name), df_mfasta, df_features, df_library))  #read_pickle takes as input as compressed file or a dataframe. Here it is a dataframe.
    mapping_hisat_mixed = pd.read_pickle(process_bam(parse_BAM(bamHISATmix.name), df_mfasta, df_features, df_library))
    mapping_star_all    = pd.read_pickle(process_bam(parse_BAM(bamSTARall.name), df_mfasta, df_features, df_library))
    mapping_star_mixed  = pd.read_pickle(process_bam(parse_BAM(bamSTARmix.name), df_mfasta, df_features, df_library))
    names = ['All with introns - Hisat2','Mix-states contigs - Hisat2','All with introns - STAR','Mix-states contigs - STAR'] 
    colors = {'All with introns - Hisat2':"limegreen",
          'Mix-states contigs - Hisat2':"forestgreen",
          'All with introns - STAR':"darkorange",
          'Mix-states contigs - STAR':"chocolate"}
    html += get_html_split(mapping_hisat_all, mapping_hisat_mixed, mapping_star_all, mapping_star_mixed, names, colors)
    '''

    #Counting table and barplots of mapped covering reads' main characteristics
    '''
    library = pd.read_pickle(df_library)
    html += get_html_mapping_descr(names, mapping_hisat_all, mapping_hisat_mixed, mapping_star_all, mapping_star_mixed, colors, library)
    fig , table = plot_covering_reads(*zip(names,[
                    mapping_hisat_all.loc[lambda df : df.covering == True,:],
                    mapping_hisat_mixed.loc[lambda df : df.covering == True,:],
                    mapping_star_all.loc[lambda df : df.covering == True,:],
                    mapping_star_mixed.loc[lambda df : df.covering == True,:]
                    ]),
                colors=colors,
                library=library)
    global_stats_table=[]
    global_stats_table['0titre']=table[0,1]
    html += get_html_table_descr(global_stats_table)
    '''

    ## SPLITREADSEARCH STAT
    if candidat:
        df_candidat = parse_candidat(candidat.name)
        print('candidat head :' , df_candidat.head(5))
        global_stat_candidat = dict()
        global_stat_candidat["0Number of candidats"] = df_candidat.shape[0]
        global_stat_candidat["1Mean length"] = 0
        nbPASS = 0
        for i,row in df_candidat.iterrows():
            global_stat_candidat["1Mean length"] += row['end'] - row['start'] + 1
            if "PASS" in row['filter']:
                nbPASS += 1
        global_stat_candidat["1Mean length"] /= (global_stat_candidat["0Number of candidats"])
        global_stat_candidat["1Mean length"] = round(global_stat_candidat["1Mean length"], 2)
        global_stat_candidat["2Mean depth"]= round(df_candidat['depth'].mean(), 2)
        global_stat_candidat["3Number of canonic junction"] = nbPASS
        html += get_html_candidat_descr(global_stat_candidat)
            
    # FOOTER
    html += "<br/>"
    html += get_html_footer()

    with open(output_file, "w") as f:
        f.write(html)

if __name__ == '__main__' :
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='fasta')
    parser.add_argument('-m','--modifiedfasta', type=argparse.FileType('r'), required=True, dest='mfasta')
    parser.add_argument('-g','--gtf', type=argparse.FileType('r'), required=True, dest='gtf')
    parser.add_argument('-1','--R1', type=argparse.FileType('r'), required=True, dest='r1')
    parser.add_argument('-2','--R2', type=argparse.FileType('r'), required=False, dest='r2')
    parser.add_argument('-a','--flagstat', type=argparse.FileType('r'), required=False, dest='flagstat')
    #flgSTARintrons:str, flgSTAR:str, flgHISAT2introns:str, flgHISAT2:str
    parser.add_argument('-S','--flgSTARintrons', type=argparse.FileType('r'), required=True, dest='flgSTARintrons')
    #parser.add_argument('-s','--flgSTAR', type=argparse.FileType('r'), required=True, dest='flgSTAR')
    parser.add_argument('-H','--flgHISAT2introns', type=argparse.FileType('r'), required=True, dest='flgHISAT2introns')
    #parser.add_argument('-h','--flgHISAT2', type=argparse.FileType('r'), required=True, dest='flgHISAT2')
    parser.add_argument('-bha','--mappingHISATall', type=argparse.FileType('r'), required=True, dest='bamHISATall')
    parser.add_argument('-bhm','--mappingHISATmix', type=argparse.FileType('r'), required=True, dest='bamHISATmix')
    parser.add_argument('-bsa','--mappingSTARall', type=argparse.FileType('r'), required=True, dest='bamSTARall')
    parser.add_argument('-bsm','--mappingSTARmix', type=argparse.FileType('r'), required=True, dest='bamSTARmix')
    #parser.add_argument('-awi','--FRS_without_introns_assemblathon', type=argparse.FileType('r'), required=False, dest='wointronsassemblathon')
    parser.add_argument('-aia','--FRS_all_modified_assemblathon', type=argparse.FileType('r'), required=True, dest='allwiintronsassemblathon') 
    #parser.add_argument('-aim','--FRS_Mixed_states_assemblathon', type=argparse.FileType('r'), required=False, dest='halfwiintronsassemblathon')
    parser.add_argument('-r','--ranksfile', type=argparse.FileType('r'), required=True, dest='ranksfile')
    parser.add_argument('-c','--candidat', type=argparse.FileType('r'), required=True, dest='candidat')
    parser.add_argument('-o','--output', type=str, required=True, dest='output')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser.add_argument('-t','--threads', type=int, default=1, required=False, dest='threads')
    parser.add_argument('-F', '--force', action='store_true', default=False, dest='force')

    args = vars(parser.parse_args())
    
    simulationReport(**args)
