#!/usr/bin/env python3

#<IntronSeeker searches introns by splice-realigning reads on contigs.>
#Copyright (C) <2019> INRAE
#<Sarah Maman, Philippe Bardou, Emilien Lasguignes, Faustine Oudin, Floréal Cabanettes, Christophe Klopp>
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.


import os
import argparse 
from argparse import ArgumentParser
import configparser # To parse parameters file
import numpy as np  # For Split read signal analysis
import pandas as pd
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import pickle
import glob
import json
import csv
import subprocess as sp # To run subprocess
import concurrent.futures as prl # For Split read signal analysis
from itertools import repeat     # For Split read signal analysis
# Import all functions from internal modules
from simulation2HTMLparse import *
from simulation2HTMLtags import *
from simulation2HTMLplots import *


# source activate ISeeker_environment;
# cd scripts/; 
# python3 simulation2HTML.py -F --config_file ../config/simulation2HTML_example.cfg


############
# SUB MAIN #
############
def simulationReport(   config_file: str,fasta:str, mfasta:str, gtf:str, r1:str, r2:str, ranks:str,
                        flagstat:str, candidat:str, split:str,
                        output:str, prefix:str, force:bool, threads:int ) :

    output_path = output + "/report"
    if prefix:
        output_path += "_" + prefix

    # Create output dir if not exist
    if not os.path.exists(output) :
        os.makedirs(output)
    
    # Output path filename report html
    output_file = output_path + ".html"
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
    
    # if simulation ?
    if mfasta:
        df_mfasta = parse_fasta(mfasta.name,True)
        # # Drop rows without not modified contigs (only sequence.modif)
    if r2 :
        df_library = parse_library(r1.name, r2.name, mfasta)
    else :
        df_library = parse_library(r1.name, mfasta)
		
    if gtf:
        df_features = parse_gtf(gtf.name)
    else:
        d = {'features': ['featureID'], 'contig': ['contig'],'feature': ['feature'],'start': [0],'end': [0],'length': [0],'flanks': ['CT_AC'],'pos_on_contig': [0]}
        df_features = pd.DataFrame(data=d)
        

    # Add to df_library the read positions on the modified contig
    # if simulation ?
    if mfasta:
        df_library['mstart'] = df_library['start']
        df_library['mend']   = df_library['end']
        for i, feature in df_features.iterrows():
            c = feature['contig'].replace('.modif','')
            s = feature['start']
            e = feature['end']
            l = feature['length']
            
            for v in df_library.loc[(df_library['contig'] == c) & (df_library['mstart'] > s)].index.values:
                df_library.at[v, 'mstart'] += l
            for v in df_library.loc[(df_library['contig'] == c) & (df_library['mend'] > s)].index.values:
                df_library.at[v, 'mend'] += l
    
    # Add a column to df_fasta with the "fasta" length (without any simulated features)
    # if simulation ?
    if mfasta and gtf:
        df_mfasta["short_length"] = df_mfasta.apply(
            compute_tr_length,
            axis = 1,
            df_features=df_features
        )
    
    # Add two columns to df_features:
    #  1- the true insertion position of the simulated feature (in term of mfasta length percentage)
    #  2- the borders of the simulated features (in term of nucleotides)

    # if simulation ?
    if mfasta and gtf:
        df_features = df_features.join(
            other = df_features.apply(
                compute_pos_on_mfasta,
                axis=1,
                df_mfasta=df_mfasta
            )
        )    

    # HEADER
    html = get_html_header()
    
    # if simulation ?
    if mfasta and gtf:
        inputfiles = [
            "Contig FASTA#" + os.path.basename(fasta.name),
            "Contig FASTA with feature(s)#" + os.path.basename(mfasta.name),
            "GTF of contig with feature(s)#" + os.path.basename(gtf.name),
             "Read1 FASTQ#" + os.path.basename(r1.name)
        ]
    else:
        inputfiles = [
            "Contig FASTA#" + os.path.basename(fasta.name),
             "Read1 FASTQ#" + os.path.basename(r1.name)
        ]   
    if r2:
        inputfiles.append("Read2 FASTQ#" + os.path.basename(r2.name))
    ranks_file=""    
    if ranks:
        ranks_file=ranks.name
        inputfiles.append("Ranks#" + os.path.basename(ranks.name))
    flagstat_file = ""
    if flagstat :
        flagstat_file = flagstat.name
        inputfiles.append("Flagstat#" + os.path.basename(flagstat_file))
    split_file=""    
    if split:
        split_file=split.name
        inputfiles.append("Split#" + os.path.basename(split.name))
    candidat_file=""    
    if candidat:
        candidat_file=candidat.name
        inputfiles.append("Candidat#" + os.path.basename(candidat.name))
    
    if mfasta:
        html += get_html_body1(flagstat_file, split_file, candidat_file, mfasta)
    else:
        html += get_html_body1(flagstat_file, split_file, candidat_file)		
   
    # INPUT FILES
    html += get_html_inputfiles(inputfiles)

    # SEQUENCE STAT
    # Global stat
    
    if gtf:
        nb_distinct_features, nb_ctg_by_feature, ctg_descr = stat_from_gtf(gtf.name)
    else:
        nb_distinct_features = 0
        nb_ctg_by_feature    = 0
        ctg_descr            = 0
    global_stat = dict()
    global_stat["001Contig FASTA - Number of sequences"]  = len(df_fasta.index)
    global_stat["002Contig FASTA - Mean sequence length"] = int(df_fasta['length'].mean())
    # if simulation ?
    if mfasta:
        global_stat["003Contig FASTA with feature(s) - Number of sequences"]  = len(df_mfasta.index)
        global_stat["004Contig FASTA with feature(s) - Mean sequence length"] = int(df_mfasta['length'].mean())
        global_stat["005Number of modified sequences"]       = df_mfasta.loc[df_mfasta.index.str.contains(".modif")].shape[0]
    if gtf:
        global_stat["006Number of distinct features in GTF"] = df_features.feature.value_counts().shape[0]
        global_stat["007Number of features in GTF"]          = len(df_features.index)
    if mfasta:
        c = 7
    else:
        c = 4	
        
    if gtf:
        for k, v in (df_features.feature.value_counts()).items() :
            global_stat["0"+str(c)+k+' (s) '] = v
            c+=1

    # ASSEMBLATHON on fasta files        
    nbContigs, totContigSize, longestContig, shortestContig, nbContigsSup1K, n50, l50, meanContigSize = run_assemblathon(fasta.name)
    global_stat_assemblathon_fasta = dict()
    global_stat_assemblathon_fasta["01Number of contigs"]         = nbContigs
    global_stat_assemblathon_fasta["02Mean contigs length"]       = round(meanContigSize, 0)
    global_stat_assemblathon_fasta["03Total size of contigs"]     = totContigSize
    global_stat_assemblathon_fasta["04Longest contig"]            = longestContig
    global_stat_assemblathon_fasta["05Shortest contig"]           = shortestContig
    global_stat_assemblathon_fasta["06Number of contigs > 1K nt"] = nbContigsSup1K
    global_stat_assemblathon_fasta["07N50 contig length"]         = n50
    global_stat_assemblathon_fasta["08L50 contig count"]          = l50
    
    # if simulation ?
    if mfasta:
        nbContigs, totContigSize, longestContig, shortestContig, nbContigsSup1K, n50, l50, meanContigSize = run_assemblathon(mfasta.name)
        global_stat_assemblathon_mfasta = dict()
        global_stat_assemblathon_mfasta["01Number of contigs"]         = nbContigs
        global_stat_assemblathon_mfasta["02Mean contigs length"]       = round(meanContigSize, 0)
        global_stat_assemblathon_mfasta["03Total size of contigs"]     = totContigSize
        global_stat_assemblathon_mfasta["04Longest contig"]            = longestContig
        global_stat_assemblathon_mfasta["05Shortest contig"]           = shortestContig
        global_stat_assemblathon_mfasta["06Number of contigs > 1K nt"] = nbContigsSup1K
        global_stat_assemblathon_mfasta["07N50 contig length"]         = n50
        global_stat_assemblathon_mfasta["08L50 contig count"]          = l50
    # if simulation ?
    if mfasta and gtf:
        html += get_html_seq_descr_simulation(global_stat, nb_ctg_by_feature, ctg_descr, gtf.name, df_features['pos_on_contig'], df_fasta, df_mfasta, global_stat_assemblathon_fasta, global_stat_assemblathon_mfasta)
    else:
        html += get_html_seq_descr_real(global_stat, df_fasta)
		
    # READS STAT
    # Global stat
    global_stat_fastq = dict()
    global_stat_fastq["01Number of reads"]  = len(df_library.index)
    global_stat_fastq["02Mean coverage"]    = df_library['length'].sum()
    deno = global_stat["002Contig FASTA - Mean sequence length"] * global_stat["001Contig FASTA - Number of sequences"]
    if deno:
        global_stat_fastq["02Mean coverage"] /= deno
    else :
        global_stat_fastq["02Mean coverage"] = "NaN"
    global_stat_fastq["02Mean coverage"]    = round(global_stat_fastq["02Mean coverage"], 2)
    global_stat_fastq["03Min reads length"] = df_library['length'].min()
    global_stat_fastq["04Max reads length"] = df_library['length'].max()
    global_stat_fastq["05Mean reads length"]= round(df_library['length'].mean(), 2)
    html += get_html_reads_descr(global_stat_fastq, df_library) 

    # # ABUNDANCE number of reads by contig
    # # Build a dataframe with:   
    # #   ctg
    # #   abund_perc => (number of read on this contig / number of reads) * 100
    # #   requested  => if grinder ranks output file is given ...
    # #   norm       => (((number of read on this contig / contig len) * mean len of all contigs) / number of reads) * 100
    # df_tmp = pd.DataFrame((df_library.groupby('contig').size()/len(df_library))*100, columns = ['abund_perc'])
    # df_fasta = df_fasta.assign(real=df_tmp.abund_perc.values)
    # if ranks:
    #     df_tmp = parse_rank_file(ranks_file)
    #     df_fasta = df_fasta.assign(rank=df_tmp['rank'].values)
    #     df_fasta = df_fasta.assign(waiting=df_tmp.rel_abund_perc.values)
    # df_tmp = pd.DataFrame(
    #     (((df_library.groupby('contig').size()/df_fasta['length'])*(df_fasta['length'].mean()))/df_library.shape[0])*100,
    #     columns = ['norm'])
    # df_fasta = df_fasta.assign(norm=df_tmp['norm'].values)
    # del df_tmp
    # if len(df_fasta.index) <= 100:
    #     html += get_html_abundance(df_fasta.sample, "Pourcentage of reads abundance for each contigs")
    # elif len(df_fasta.index) * 0.2 > 100:
    #     html += get_html_abundance(df_fasta.sample(100), "Pourcentage of reads abundance for 100 random contigs")
    # else:
    #     html += get_html_abundance(df_fasta.sample(frac=.2), "Pourcentage of reads abundance for 20\% \random contigs")
    
    ## ALIGNMENT STATS
    if flagstat:
        nbreads, mapped, mappercent, paired, proper, properpercent, secondary, singletons=parse_flagstat(flagstat_file)
        global_stat_flagstat = dict()
        global_stat_flagstat["01Number of reads"] = paired
        global_stat_flagstat["02Number of lines in the BAM"] = nbreads
        global_stat_flagstat["03Number of mapped reads"] = mapped
        global_stat_flagstat["04Percentage of mapped reads"] = mappercent
        global_stat_flagstat["05Number of properly paired reads"] = proper
        global_stat_flagstat["06Percentage of properly paired reads"] = properpercent
        global_stat_flagstat["07Secondary alignments"] = secondary
        global_stat_flagstat["08Singletons alignements"] = singletons
        html += get_html_flagstat_descr(global_stat_flagstat)
   
    html += get_html_results()
   
    ## SPLITREADSEARCH STAT
    if split:
        df_split=parse_split(split.name)   
        global_stat_split = dict()
        global_stat_split["01Number of reads overlapping introns"] = len(df_split.index)
        global_stat_split["02Mean length of introns"] = round(df_split['split_length'].mean(), 2)

        nbCanonic = 0
        nbOtherJunctions = 0
        df_split.sort_values(by=['split_borders'])
        c = 3
        n = 0
        for k, v in (df_split['split_borders'].value_counts()).items() :
            if k == "GT_AG" or k == "CT_AC":
                nbCanonic += v
            elif n < 5:
                global_stat_split["0"+str(c)+"Junction "+k] = v
                n += 1
                c += 1
            else:
                nbOtherJunctions += v
        global_stat_split["0"+str(c)+"Other junctions"] = nbOtherJunctions
        global_stat_split["0"+str(c+1)+"Canonical junction (GT_AG or CT_AC)"] = nbCanonic

        html += get_html_split_descr(global_stat_split)   

    ## CANDIDATS statistics - detected introns
    if candidat:
        df_candidat, mindepth, maxlen = parse_candidat(candidat.name)

        # Definition dict
        definitions = dict()
        definitions['DP']   = "Filtered because of depth (<= "+ str(mindepth)+ ")"
        definitions['LEN']  = "Filtered because of length (>= "+ str(maxlen)+ "%)"    
        definitions['SS']   = "Filtered because of non canonical junction"
        definitions['OI']   = "Filtered because of overlapping introns"
        definitions['PASS'] = "Number"
        
        # Detected introns (stat table and depth distribution)
        global_stat_detected_introns = dict()
        global_stat_detected_introns["01Number"] = len(df_candidat.index)
        global_stat_detected_introns["02Min length"]  = df_candidat.iloc[0]['end'] - df_candidat.iloc[0]['start'] + 1 
        global_stat_detected_introns["03Max length"]  = 0
        global_stat_detected_introns["04Mean length"] = 0
        for i,row in df_candidat.iterrows():
            l = row['end'] - row['start'] + 1
            global_stat_detected_introns["02Min length"] = min(l, global_stat_detected_introns["02Min length"])
            global_stat_detected_introns["03Max length"] = max(l, global_stat_detected_introns["03Max length"])
            global_stat_detected_introns["04Mean length"] += row['end'] - row['start'] + 1
        deno = global_stat_detected_introns["01Number"]
        if deno :
            global_stat_detected_introns["04Mean length"] /= deno
        else :
            global_stat_detected_introns["04Mean length"] = "NaN"
        global_stat_detected_introns["04Mean length"] = round(global_stat_detected_introns["04Mean length"], 2)
        global_stat_detected_introns["05Min depth"]  = df_candidat['depth'].min()
        global_stat_detected_introns["06Max depth"]  = df_candidat['depth'].max()
        global_stat_detected_introns["07Mean depth"]  = round(df_candidat['depth'].mean(), 2)
        html += get_html_detected(global_stat_detected_introns, df_candidat)

        # Filtered detected introns
        global_stat_f_detected_introns = dict()
        global_stat_f_detected_introns["01" + definitions['PASS']] = 0
        global_stat_f_detected_introns["02" + definitions['DP']]   = 0
        global_stat_f_detected_introns["03" + definitions['LEN']]  = 0
        global_stat_f_detected_introns["04" + definitions['SS']]   = 0
        global_stat_f_detected_introns["05" + definitions['OI']]   = 0
        global_stat_f_detected_introns["06Min length"]  = 99999999
        global_stat_f_detected_introns["07Max length"]  = 0 
        global_stat_f_detected_introns["08Mean length"] = 0
        global_stat_f_detected_introns["09Min depth"]   = 99999999
        global_stat_f_detected_introns["10Max depth"]   = 0 
        global_stat_f_detected_introns["11Mean depth"]  = 0
        for i, v in (df_candidat['filter'].items()) :
            if "PASS" in v:
                global_stat_f_detected_introns["01" + definitions['PASS']] += 1
                l = df_candidat.loc[i]['end'] - df_candidat.loc[i]['start'] + 1
                global_stat_f_detected_introns["06Min length"] = min(l, global_stat_f_detected_introns["06Min length"])
                global_stat_f_detected_introns["07Max length"] = max(l, global_stat_f_detected_introns["07Max length"])
                global_stat_f_detected_introns["08Mean length"]+= l
                l = df_candidat.loc[i]['depth']
                global_stat_f_detected_introns["09Min depth"]  = min(l, global_stat_f_detected_introns["09Min depth"])
                global_stat_f_detected_introns["10Max depth"]  = max(l, global_stat_f_detected_introns["10Max depth"])
                global_stat_f_detected_introns["11Mean depth"] += l
            else :
                if "DP" in v:
                    global_stat_f_detected_introns["02" + definitions['DP']] += 1
                if "LEN" in v:
                        global_stat_f_detected_introns["03" + definitions['LEN']] += 1
                if "SS" in v:
                    global_stat_f_detected_introns["04" + definitions['SS']] += 1
                if "OI" in v:
                    global_stat_f_detected_introns["05" + definitions['OI']] += 1

        deno = global_stat_f_detected_introns["01" + definitions['PASS']]
        if deno :
            global_stat_f_detected_introns["08Mean length"] /= deno
            global_stat_f_detected_introns["08Mean length"]  = round(global_stat_f_detected_introns["08Mean length"], 2)
            global_stat_f_detected_introns["11Mean depth"]  /= deno
            global_stat_f_detected_introns["11Mean depth"]   = round(global_stat_f_detected_introns["11Mean depth"], 2)
        else :
            global_stat_f_detected_introns["08Mean length"] = "NaN"
            global_stat_f_detected_introns["11Mean depth"]  = "NaN"
                          
        if mfasta and gtf:
            # Detectable features (filter features because of threshold: mindepth and maxlength)
            # Add a column to df_features for the DP (using df_library)
            df_features["depth"] = df_features.apply(
                compute_dp,
                axis = 1,
                df_library=df_library
            )
            # Add a column to df_features for the length of the corresponding contig (using df_mfasta)
            df_features["ctg_length"] = df_features.apply(
                compute_len,
                axis = 1,
                df_mfasta=df_mfasta
            )
            global_stat_detectable_features = dict()
            global_stat_detectable_features["01" + definitions['PASS']] = 0
            global_stat_detectable_features["02" + definitions['DP']]   = 0
            global_stat_detectable_features["03" + definitions['LEN']]  = 0
            global_stat_detectable_features["04" + definitions['SS']]   = 0
            global_stat_detectable_features["05" + definitions['OI']]   = "-"
            global_stat_detectable_features["06Min length"]  = 99999999
            global_stat_detectable_features["07Max length"]  = 0 
            global_stat_detectable_features["08Mean length"] = 0
            global_stat_detectable_features["09Min depth"]   = 99999999
            global_stat_detectable_features["10Max depth"]   = 0 
            global_stat_detectable_features["11Mean depth"]  = 0
            for index, row in df_features.iterrows():            
                PASS = True
                if (row['flanks'] != 'GT_AG' and row['flanks'] != 'CT_AC'): 
                    global_stat_detectable_features["04" + definitions['SS']] += 1    
                    PASS = False 
                if ((((row['end'] - row['start'] + 1) / row['ctg_length']) * 100) >= int(maxlen)):
                    global_stat_detectable_features["03" + definitions['LEN']] += 1
                    PASS = False
                if (row['depth'] <= int(mindepth)):
                    global_stat_detectable_features["02" + definitions['DP']] += 1
                    PASS = False
                if PASS:
                    global_stat_detectable_features["01" + definitions['PASS']] += 1
                    l = row['end'] - row['start'] + 1
                    global_stat_detectable_features["06Min length"] = min(l, global_stat_detectable_features["06Min length"])
                    global_stat_detectable_features["07Max length"] = max(l, global_stat_detectable_features["07Max length"])
                    global_stat_detectable_features["08Mean length"]+= l
                    l = row['depth']
                    global_stat_detectable_features["09Min depth"]  = min(l, global_stat_detectable_features["09Min depth"])
                    global_stat_detectable_features["10Max depth"]  = max(l, global_stat_detectable_features["10Max depth"])
                    global_stat_detectable_features["11Mean depth"] += l
            deno = global_stat_detectable_features["01" + definitions['PASS']]
            if deno :
                global_stat_detectable_features["08Mean length"] /= deno
                global_stat_detectable_features["08Mean length"] =  round(global_stat_detectable_features["08Mean length"], 2)
                global_stat_detectable_features["11Mean depth"]  /= deno
                global_stat_detectable_features["11Mean depth"]  =  round(global_stat_detectable_features["11Mean depth"], 2)
            else :
                global_stat_detectable_features["08Mean length"] = "NaN"
                global_stat_detectable_features["11Mean depth"]  = "NaN"
            html += get_html_candidat(global_stat_f_detected_introns, global_stat_detectable_features)
        else :
            html += get_html_candidat(global_stat_f_detected_introns)
 

        # Detected : contig name / filtered_detected_too_complex
        global_stat_too_complex_detected = dict()
        
        df_too_complex_detected = df_candidat.loc[df_candidat['filter'] == "PASS"][['reference']].groupby(['reference']) \
                             .size() \
                             .nlargest(10) \
                             .reset_index(name='top10')  
        
        # Add filter information for each complex intron in top 10 of contigs with the highest number of detected introns
        cmp = 0
        for k, v in df_too_complex_detected['reference'].items() :
            nbPASS =len(df_candidat[(df_candidat['reference'] == v) & (df_candidat['filter'] == "PASS" )] )
            nbDP   =len(df_candidat[(df_candidat['reference'] == v) & (df_candidat['filter'] == "DP" )] ) 
            nbOI   =len(df_candidat[(df_candidat['reference'] == v) & (df_candidat['filter'] == "OI" )] ) 
            nbSS   =len(df_candidat[(df_candidat['reference'] == v) & (df_candidat['filter'] == "SS" )] ) 
            nbLEN  =len(df_candidat[(df_candidat['reference'] == v) & (df_candidat['filter'] == "LEN" )] ) 
            global_stat_too_complex_detected["0"+str(cmp)+str(v)] = str(nbPASS)+" PASS; "+str(nbDP)+" DP; "+str(nbOI)+" IO; "+str(nbSS)+" SS; "+str(nbLEN)+" LEN"
            cmp += 1    
         
        html += get_html_too_complex(global_stat_too_complex_detected)

        # if simulation ?
        if mfasta and  gtf:
            eval_def = dict()
            eval_def["TP"] = "TP = Detected introns &#8745; Features"
            eval_def["FP"] = "FP = Detected introns &#8713; Features"
            eval_def["FN"] = "FN = Features &#8713; Detected introns"
            eval_def["Se"] = "Sensibility (Se) = TP / (TP+FN)"
            eval_def["Sp"] = "Specificity (Sp) = TP / (TP+FP)"
            eval_def["F1"] = "F1 score = (2*Se*Sp) / (Se+Sp)"

            eval_stat = dict()
            eval_stat["01Number of detected introns"] = global_stat_detected_introns["01Number"] 
            eval_stat["02Number of features"] = global_stat["007Number of features in GTF"] 
            eval_stat["03"+eval_def["TP"]] = 0
            eval_stat["04"+eval_def["FP"]] = 0
            eval_stat["05"+eval_def["FN"]] = 0
            eval_stat["06"+eval_def["Se"]] = 0
            eval_stat["07"+eval_def["Sp"]] = 0
            eval_stat["08"+eval_def["F1"]] = 0
            for index, row in df_candidat.iterrows():
                ok = len(df_features.loc[lambda df :
                            (df['contig'] == row['reference']) &
                            (df['start']  == row['start']) &
                            (df['end']    == row['end'])])  
                if ok == 1:
                    eval_stat["03"+eval_def["TP"]] += 1         
                else:
                    eval_stat["04"+eval_def["FP"]] += 1
                    
            eval_stat["05"+eval_def["FN"]] = eval_stat["02Number of features"] - eval_stat["03"+eval_def["TP"]]
            
            
            deno = eval_stat["03"+eval_def["TP"]] + eval_stat["05"+eval_def["FN"]]
            if deno :
                eval_stat["06"+eval_def["Se"]] = eval_stat["03"+eval_def["TP"]] / deno * 100
            else :
                eval_stat["06"+eval_def["Se"]] = "NaN"
            deno = (eval_stat["03"+eval_def["TP"]] + eval_stat["04"+eval_def["FP"]])
            if deno :
                eval_stat["07"+eval_def["Sp"]] = eval_stat["03"+eval_def["TP"]] / deno * 100
            else :
                eval_stat["07"+eval_def["Sp"]] = "NaN"
            if eval_stat["06"+eval_def["Se"]] != "NaN" and eval_stat["07"+eval_def["Sp"]] != "NaN" :
                deno = eval_stat["06"+eval_def["Se"]] + eval_stat["07"+eval_def["Sp"]]
                if deno :
                    eval_stat["08"+eval_def["F1"]] = round(2*(eval_stat["06"+eval_def["Se"]]*eval_stat["07"+eval_def["Sp"]]) / deno, 2) 
                else :
                    eval_stat["08"+eval_def["F1"]] = "NaN"
            else :
                eval_stat["08"+eval_def["F1"]] = "NaN"
            if isinstance(eval_stat["06"+eval_def["Se"]], float) :
                eval_stat["06"+eval_def["Se"]] = round(eval_stat["06"+eval_def["Se"]], 2)
            if isinstance(eval_stat["07"+eval_def["Sp"]], float) :
                eval_stat["07"+eval_def["Sp"]] = round(eval_stat["07"+eval_def["Sp"]], 2)


            eval_f_stat = dict()
            eval_f_stat["01Number of detected introns"] = global_stat_f_detected_introns["01" + definitions['PASS']]
            eval_f_stat["02Number of features"] = global_stat_detectable_features["01" + definitions['PASS']]
            eval_f_stat["03"+eval_def["TP"]] = 0
            eval_f_stat["04"+eval_def["FP"]] = 0
            eval_f_stat["05"+eval_def["FN"]] = 0
            eval_f_stat["06"+eval_def["Se"]] = 0
            eval_f_stat["07"+eval_def["Sp"]] = 0
            eval_f_stat["08"+eval_def["F1"]] = 0
            
            if gtf:
                print('features head :', df_features.head(5), '\n\n')
            print('candidat head :', df_candidat.head(5), '\n\n')
            print('library head :', df_library, '\n\n')
            for index, row in df_candidat.iterrows():
                if "PASS" in row['filter']:
                    ok = len(df_features.loc[lambda df :
                             (df['contig'] == row['reference']) &
                             (df['start']  == row['start']) &
                             (df['end']    == row['end']) &
                             (df['depth'] > int(mindepth)) &
                             ((((df['end'] - df['start'] + 1) / df['ctg_length']) * 100) < int(maxlen)) &
                             (df['flanks'].str.contains('GT_AG|CT_AC'))])
                    if ok == 1:
                        eval_f_stat["03"+eval_def["TP"]] += 1
                    else:
                        eval_f_stat["04"+eval_def["FP"]] += 1
                       
            eval_f_stat["05"+eval_def["FN"]] = eval_f_stat["02Number of features"] - eval_f_stat["03"+eval_def["TP"]]
            

            deno = eval_f_stat["03"+eval_def["TP"]] + eval_f_stat["05"+eval_def["FN"]]
            if deno :
                eval_f_stat["06"+eval_def["Se"]] = eval_f_stat["03"+eval_def["TP"]] / deno * 100
            else :
                eval_f_stat["06"+eval_def["Se"]] = "NaN"
            deno = eval_f_stat["03"+eval_def["TP"]] + eval_f_stat["04"+eval_def["FP"]]
            if deno :
                eval_f_stat["07"+eval_def["Sp"]] = eval_f_stat["03"+eval_def["TP"]] / deno * 100
            else :
                eval_f_stat["07"+eval_def["Sp"]] = "NaN"
            if eval_f_stat["06"+eval_def["Se"]] != "NaN" and eval_f_stat["07"+eval_def["Sp"]] != "NaN" :
                deno = eval_f_stat["06"+eval_def["Se"]] + eval_f_stat["07"+eval_def["Sp"]]
                if deno :
                    eval_f_stat["08"+eval_def["F1"]] = round(2*(eval_f_stat["06"+eval_def["Se"]]*eval_f_stat["07"+eval_def["Sp"]]) / deno, 2) 
                else :
                    eval_f_stat["08"+eval_def["F1"]] = "NaN"
            else :
                eval_f_stat["08"+eval_def["F1"]] = "NaN"
            if isinstance(eval_f_stat["06"+eval_def["Se"]], float) :
                eval_f_stat["06"+eval_def["Se"]] = round(eval_f_stat["06"+eval_def["Se"]], 2)
            if isinstance(eval_f_stat["07"+eval_def["Sp"]], float) :
                eval_f_stat["07"+eval_def["Sp"]] = round(eval_f_stat["07"+eval_def["Sp"]], 2)

            html += get_html_eval(eval_stat, eval_f_stat)
    
    # GLOSSARY
    html += get_html_glossary()

    # FOOTER
    html += "<br/>"
    html += get_html_footer()

    with open(output_file, "w") as f:
        f.write(html)

    # print('DATAFRAMES:\n\n')
    # print('fasta head :' , df_fasta.head(5), '\n\n')
    # if mfasta:
    #     print('mfasta head :' , df_mfasta.head(5), '\n\n')    
    # # print('mfasta head :' , df_mfasta.head(5), '\n\n')    
    # print('library head :' , df_library.head(5), '\n\n')
    # print('features head :', df_features.head(5), '\n\n')
    # if split:
    #     print('df_split', df_split, '\n\n')
    # if candidat:
    #     print('df_candidat', df_candidat, '\n\n')


if __name__ == '__main__' :
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--config_file')
    args, left_argv = parser.parse_known_args()
    if args.config_file:
        with open(args.config_file, 'r') as f:
            config = configparser.ConfigParser()
            config.read([args.config_file])
    
    parser = ArgumentParser()
    parser.add_argument('--config-file', type=argparse.FileType('r'), required=False, help="Provide a config file")
    parser.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='fasta', help="Path to the reference FASTA file.")
    parser.add_argument('-m','--modified-fasta', type=argparse.FileType('r'), required=False, dest='mfasta', help="Path to the modified FASTA file.")
    parser.add_argument('-g','--gtf', type=argparse.FileType('r'), required=False, dest='gtf', help="GTF filename which contains the genome annotation.")
    parser.add_argument('-1','--R1', type=argparse.FileType('r'), required=True, dest='r1', help="Name of the  FASTQ  file  which  contains  the  single-end   reads library. If paired-end, filename of #1 reads mates")
    parser.add_argument('-2','--R2', type=argparse.FileType('r'), required=False, dest='r2', help="Only for a paired-end library, filename of #2 reads mates.")
    parser.add_argument('--flagstat', type=argparse.FileType('r'), required=False, dest='flagstat', help="Path to flagstat file.")
    parser.add_argument('-r','--ranks', type=argparse.FileType('r'), required=False, dest='ranks', help="Path to ranks file.")
    parser.add_argument('-c','--candidat', type=argparse.FileType('r'), required=False, dest='candidat', help="Path to candidat file.")
    parser.add_argument('-s','--split', type=argparse.FileType('r'), required=False, dest='split', help="Path to split file.")
    parser.add_argument('-o','--output', type=str, required=True, dest='output', help="Output dir name.")
    parser.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix', help="Prefix for output files name.")
    parser.add_argument('-t','--threads', type=int, default=1, required=False, dest='threads', help="Number of threads [1]")
    parser.add_argument('-F', '--force', action='store_true', default=False, dest='force', help="Force to overwrite output files.")

    try:
        config
    except NameError:
        pass
    else:
        for k, v in config.items("Defaults"):
            config_args={str(k): str(v)}
            # Use values from configuration file by default
            parser.set_defaults(**config_args)
            # Reset `required` attribute when provided from config file
            for action in parser._actions:
                if action.dest in config_args:
                    action.required = False

    # Override with command line arguments when provided
    args = vars(parser.parse_args(left_argv, args))

    simulationReport(**args)
