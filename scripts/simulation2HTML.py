#!/usr/bin/env python3

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
# scp  /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML/*FRS_CASA_sample1_n1000_r_STAR*.html smaman@genologin.toulouse.inra.fr:/save/smaman/public_html/intronSeeker/.
# See result : http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/report_FRS_CASA_sample1_n1000_r_STAR_simulation.html

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
    # SARAH / CAS REAL
    if mfasta:
        df_mfasta = parse_fasta(mfasta.name,True)
        # # Drop rows without not modified contigs (only sequence.modif)
        # df_mfasta = df_mfasta[df_mfasta.index.str.contains("modif")]
    # df_mfasta = parse_fasta(mfasta.name, True)    
    # SARAH / CAS REAL

    if r2 :
        df_library = parse_library(r1.name, r2.name)
    else :
        df_library = parse_library(r1.name)

    df_features = parse_gtf(gtf.name)

    # Add a column to df_fasta with the "fasta" length (without any simulated features)
    # SARAH / CAS REAL
    if mfasta:
        df_mfasta["short_length"] = df_mfasta.apply(
            compute_tr_length,
            axis = 1,
            df_features=df_features
        )

    # df_mfasta["short_length"] = df_mfasta.apply(
    #     compute_tr_length,
    #     axis = 1,
    #     df_features=df_features
    # )
    # SARAH / CAS REAL
    
    # Add two columns to df_features:
    #  1- the true insertion position of the simulated feature (in term of mfasta length percentage)
    #  2- the borders of the simulated features (in term of nucleotides)

    # SARAH / CAS REAL
    if mfasta:
        df_features = df_features.join(
            other = df_features.apply(
                compute_pos_on_mfasta,
                axis=1,
                df_mfasta=df_mfasta
            )
        )    
    # df_features = df_features.join(
    #     other = df_features.apply(
    #         compute_pos_on_mfasta,
    #         axis=1,
    #         df_mfasta=df_mfasta
    #     )
    # )
    # SARAH / CAS REAL

    # HEADER
    html = get_html_header()
    
    inputfiles = [
        "Contig FASTA#" + os.path.basename(fasta.name),
        "Contig FASTA with feature(s)#" + os.path.basename(mfasta.name),
        "GTF of contig with feature(s)#" + os.path.basename(gtf.name),
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

    html += get_html_body1(flagstat_file, split_file, candidat_file)
   
    # INPUT FILES
    html += get_html_inputfiles(inputfiles)

    # SEQUENCE STAT
    # Global stat
    nb_distinct_features, nb_ctg_by_feature, ctg_descr = stat_from_gtf(gtf.name)
    global_stat = dict()
    global_stat["0Contig FASTA - Number of sequences"]  = len(df_fasta.index)
    global_stat["1Contig FASTA - Mean sequence length"] = int(df_fasta['length'].mean())
    global_stat["2Contig FASTA with feature(s) - Number of sequences"]  = len(df_mfasta.index)
    global_stat["3Contig FASTA with feature(s) - Mean sequence length"] = int(df_mfasta['length'].mean())
    global_stat["4Number of modified sequences"]       = df_mfasta.loc[df_mfasta.index.str.contains(".modif")].shape[0]
    global_stat["5Number of distinct features in GTF"] = df_features.feature.value_counts().shape[0]
    global_stat["6Number of features in GTF"]          = len(df_features.index)
    c = 7
    for k, v in (df_features.feature.value_counts()).items() :
        global_stat[str(c)+k] = v
        c+=1

    # ASSEMBLATHON on fasta files        
    nbContigs, totContigSize, longestContig, shortestContig, nbContigsSup1K, n50, l50, meanContigSize = run_assemblathon(fasta.name)
    global_stat_assemblathon_fasta = dict()
    global_stat_assemblathon_fasta["0Number of contigs"]         = nbContigs
    global_stat_assemblathon_fasta["1Mean contigs length"]       = round(meanContigSize, 0)
    global_stat_assemblathon_fasta["2Total size of contigs"]     = totContigSize
    global_stat_assemblathon_fasta["3Longest contig"]            = longestContig
    global_stat_assemblathon_fasta["4Shortest contig"]           = shortestContig
    global_stat_assemblathon_fasta["5Number of contigs > 1K nt"] = nbContigsSup1K
    global_stat_assemblathon_fasta["6N50 contig length"]         = n50
    global_stat_assemblathon_fasta["7L50 contig count"]          = l50

    if mfasta:
        nbContigs, totContigSize, longestContig, shortestContig, nbContigsSup1K, n50, l50, meanContigSize = run_assemblathon(mfasta.name)
        global_stat_assemblathon_mfasta = dict()
        global_stat_assemblathon_mfasta["0Number of contigs"]         = nbContigs
        global_stat_assemblathon_mfasta["1Mean contigs length"]       = round(meanContigSize, 0)
        global_stat_assemblathon_mfasta["2Total size of contigs"]     = totContigSize
        global_stat_assemblathon_mfasta["3Longest contig"]            = longestContig
        global_stat_assemblathon_mfasta["4Shortest contig"]           = shortestContig
        global_stat_assemblathon_mfasta["5Number of contigs > 1K nt"] = nbContigsSup1K
        global_stat_assemblathon_mfasta["6N50 contig length"]         = n50
        global_stat_assemblathon_mfasta["7L50 contig count"]          = l50
    html += get_html_seq_descr(global_stat, nb_ctg_by_feature, ctg_descr, gtf.name, df_features['pos_on_contig'], df_fasta, df_mfasta, global_stat_assemblathon_fasta, global_stat_assemblathon_mfasta)

    # READS STAT
    # Global stat
    global_stat_fastq = dict()
    global_stat_fastq["0Number of reads"]  = len(df_library.index)
    global_stat_fastq["1Mean coverage"]    = df_library['length'].sum()
    deno = global_stat["1Contig FASTA - Mean sequence length"] * global_stat["0Contig FASTA - Number of sequences"]
    if deno:
        global_stat_fastq["1Mean coverage"] /= deno
    else :
        global_stat_fastq["1Mean coverage"] = "NaN"
    global_stat_fastq["1Mean coverage"]    = round(global_stat_fastq["1Mean coverage"], 2)
    global_stat_fastq["2Min reads length"] = df_library['length'].min()
    global_stat_fastq["3Max reads length"] = df_library['length'].max()
    global_stat_fastq["4Mean reads length"]= round(df_library['length'].mean(), 2)
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
        global_stat_flagstat["0Number of reads"] = paired
        global_stat_flagstat["1Number of lines in the BAM"] = nbreads
        global_stat_flagstat["2Number of mapped reads"] = mapped
        global_stat_flagstat["3Percentage of mapped reads"] = mappercent
        global_stat_flagstat["4Number of properly paired reads"] = proper
        global_stat_flagstat["5Percentage of properly paired reads"] = properpercent
        global_stat_flagstat["6Secondary alignments"] = secondary
        global_stat_flagstat["7Singletons alignements"] = singletons
        html += get_html_flagstat_descr(global_stat_flagstat)
   
    html += get_html_results()
   
    ## SPLITREADSEARCH STAT
    if split:
        df_split=parse_split(split.name)   
        global_stat_split = dict()
        global_stat_split["0Number of reads overlapping introns"] = len(df_split.index)
        global_stat_split["1Mean length of introns"] = round(df_split['split_length'].mean(), 2)

        nbCanonic = 0
        nbOtherJunctions = 0
        df_split.sort_values(by=['split_borders'])
        c = 2
        n = 0
        for k, v in (df_split['split_borders'].value_counts()).items() :
            if k == "GT_AG" or k == "CT_AC":
                nbCanonic += v
            elif n < 6:
                global_stat_split[str(c)+"Junction "+k] = v
                n += 1
                c += 1
            else:
                nbOtherJunctions += v
        global_stat_split[str(c)+"Other junctions"] = nbOtherJunctions
        global_stat_split[str(c+1)+"Canonical junction (GT_AG or CT_AC)"] = nbCanonic
        
        html += get_html_split_descr(global_stat_split)   

    ## CANDIDATS statistics - detected introns
    if candidat:
        df_candidat, mindepth, maxlen = parse_candidat(candidat.name)

        # Definition dict
        definitions = dict()
        definitions['DP']   = "Filtered because of depth (<= "+ str(mindepth)+ ")"
        definitions['LEN']  = "Filtered because of length (>= "+ str(maxlen)+ "%)"    
        definitions['SS']   = "Filtered because of non canonical junction"
        definitions['PASS'] = "Number"
        
        # Detected introns (stat table and depth distribution)
        global_stat_detected_introns = dict()
        global_stat_detected_introns["0Number"] = len(df_candidat.index)
        global_stat_detected_introns["1Min length"]  = df_candidat.iloc[0]['end'] - df_candidat.iloc[0]['start'] + 1 
        global_stat_detected_introns["2Max length"]  = 0
        global_stat_detected_introns["3Mean length"] = 0
        for i,row in df_candidat.iterrows():
            l = row['end'] - row['start'] + 1
            global_stat_detected_introns["1Min length"] = min(l, global_stat_detected_introns["1Min length"])
            global_stat_detected_introns["2Max length"] = max(l, global_stat_detected_introns["2Max length"])
            global_stat_detected_introns["3Mean length"] += row['end'] - row['start'] + 1
        deno = global_stat_detected_introns["0Number"]
        if deno :
            global_stat_detected_introns["3Mean length"] /= deno
        else :
            global_stat_detected_introns["3Mean length"] = "NaN"
        global_stat_detected_introns["3Mean length"] = round(global_stat_detected_introns["3Mean length"], 2)
        global_stat_detected_introns["4Min depth"]  = df_candidat['depth'].min()
        global_stat_detected_introns["5Max depth"]  = df_candidat['depth'].max()
        global_stat_detected_introns["6Mean depth"]  = round(df_candidat['depth'].mean(), 2)
        html += get_html_detected(global_stat_detected_introns, df_candidat)

        # Filtered detected introns
        global_stat_f_detected_introns = dict()
        global_stat_f_detected_introns["0" + definitions['PASS']] = 0
        global_stat_f_detected_introns["1" + definitions['DP']]   = 0
        global_stat_f_detected_introns["2" + definitions['LEN']]  = 0
        global_stat_f_detected_introns["3" + definitions['SS']]   = 0
        global_stat_f_detected_introns["4Min length"]  = 99999999
        global_stat_f_detected_introns["5Max length"]  = 0 
        global_stat_f_detected_introns["6Mean length"] = 0
        global_stat_f_detected_introns["7Min depth"]   = 99999999
        global_stat_f_detected_introns["8Max depth"]   = 0 
        global_stat_f_detected_introns["9Mean depth"]  = 0
        for i, v in (df_candidat['filter'].items()) :
            if "PASS" in v:
                global_stat_f_detected_introns["0" + definitions['PASS']] += 1
                l = df_candidat.loc[i]['end'] - df_candidat.loc[i]['start'] + 1
                global_stat_f_detected_introns["4Min length"] = min(l, global_stat_f_detected_introns["4Min length"])
                global_stat_f_detected_introns["5Max length"] = max(l, global_stat_f_detected_introns["5Max length"])
                global_stat_f_detected_introns["6Mean length"]+= l
                l = df_candidat.loc[i]['depth']
                global_stat_f_detected_introns["7Min depth"]  = min(l, global_stat_f_detected_introns["7Min depth"])
                global_stat_f_detected_introns["8Max depth"]  = max(l, global_stat_f_detected_introns["8Max depth"])
                global_stat_f_detected_introns["9Mean depth"] += l
            else :
                if "DP" in v:
                    global_stat_f_detected_introns["1" + definitions['DP']] += 1
                if "LEN" in v:
                        global_stat_f_detected_introns["2" + definitions['LEN']] += 1
                if "SS" in v:
                    global_stat_f_detected_introns["3" + definitions['SS']] += 1
        deno = global_stat_f_detected_introns["0" + definitions['PASS']]
        if deno :
            global_stat_f_detected_introns["6Mean length"] /= deno
            global_stat_f_detected_introns["6Mean length"]  = round(global_stat_f_detected_introns["6Mean length"], 2)
            global_stat_f_detected_introns["9Mean depth"]  /= deno
            global_stat_f_detected_introns["9Mean depth"]   = round(global_stat_f_detected_introns["9Mean depth"], 2)
        else :
            global_stat_f_detected_introns["6Mean length"] = "NaN"
            global_stat_f_detected_introns["9Mean depth"]  = "NaN"

          
       
        # if simulation ?                   
        if mfasta :
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
            global_stat_detectable_features["0" + definitions['PASS']] = 0
            global_stat_detectable_features["1" + definitions['DP']]   = 0
            global_stat_detectable_features["2" + definitions['LEN']]  = 0
            global_stat_detectable_features["3" + definitions['SS']]   = 0
            global_stat_detectable_features["4Min length"]  = 99999999
            global_stat_detectable_features["5Max length"]  = 0 
            global_stat_detectable_features["6Mean length"] = 0
            global_stat_detectable_features["7Min depth"]   = 99999999
            global_stat_detectable_features["8Max depth"]   = 0 
            global_stat_detectable_features["9Mean depth"]  = 0
            for index, row in df_features.iterrows():            
                PASS = True
                if (row['flanks'] != 'GT_AG' and row['flanks'] != 'CT_AC'): 
                    global_stat_detectable_features["3" + definitions['SS']] += 1    
                    PASS = False 
                if ((((row['end'] - row['start'] + 1) / row['ctg_length']) * 100) >= int(maxlen)):
                    global_stat_detectable_features["2" + definitions['LEN']] += 1
                    PASS = False
                if (row['depth'] <= int(mindepth)):
                    global_stat_detectable_features["1" + definitions['DP']] += 1
                    PASS = False
                if PASS:
                    global_stat_detectable_features["0" + definitions['PASS']] += 1
                    l = row['end'] - row['start'] + 1
                    global_stat_detectable_features["4Min length"] = min(l, global_stat_detectable_features["4Min length"])
                    global_stat_detectable_features["5Max length"] = max(l, global_stat_detectable_features["5Max length"])
                    global_stat_detectable_features["6Mean length"]+= l
                    l = row['depth']
                    global_stat_detectable_features["7Min depth"]  = min(l, global_stat_detectable_features["7Min depth"])
                    global_stat_detectable_features["8Max depth"]  = max(l, global_stat_detectable_features["8Max depth"])
                    global_stat_detectable_features["9Mean depth"] += l
            deno = global_stat_detectable_features["0" + definitions['PASS']]
            if deno :
                global_stat_detectable_features["6Mean length"] /= deno
                global_stat_detectable_features["6Mean length"] =  round(global_stat_detectable_features["6Mean length"], 2)
                global_stat_detectable_features["9Mean depth"]  /= deno
                global_stat_detectable_features["9Mean depth"]  =  round(global_stat_detectable_features["9Mean depth"], 2)
            else :
                global_stat_detectable_features["6Mean length"] = "NaN"
                global_stat_detectable_features["9Mean depth"]  = "NaN"
            html += get_html_candidat(global_stat_f_detected_introns, global_stat_detectable_features)
        else :
            html += get_html_candidat(global_stat_f_detected_introns)
 

        # Test Sarah
        global_stat_too_complex_detected = dict()
        # https://stackoverflow.com/questions/40454030/count-and-sort-with-pandas
        # DETECTED : Nom de contig / filtered_detected_too_complex
        df_too_complex_detected = df_candidat[['reference']].groupby(['reference']) \
                             .size() \
                             .nlargest(10) \
                             .reset_index(name='top10')        
        print('df_too_complex_detected ', df_too_complex_detected )  
        
        cmp = 0
        for k, v in df_too_complex_detected['reference'].items() :
            global_stat_too_complex_detected[str(cmp)+str(v)] = df_too_complex_detected.loc[k]['top10']
            cmp += 1
        
        df_too_complex_detected_filtered = df_candidat[['reference']].loc[df_candidat['filter'].str.contains('PASS')].groupby(['reference']) \
                             .size() \
                             .nlargest(10) \
                             .reset_index(name='top10')        
        print('df_too_complex_detected_filtered', df_too_complex_detected_filtered)

        # DETECTABLE :  nom contig detectable / nb too complex introns
        global_stat_too_complex_detectable = dict()
        if mfasta:
            df_too_complex_detectable = df_features[['contig']].groupby(['contig']) \
                             .size() \
                             .nlargest(10) \
                             .reset_index(name='top10')        
            cmp = 0
            for k, v in df_too_complex_detectable['contig'].items() :
                global_stat_too_complex_detectable[str(cmp)+str(v)] = df_too_complex_detectable.loc[k]['top10']
                cmp += 1
                     
            #Add table "too complex" in html report    
            html += get_html_too_complex(global_stat_too_complex_detected, global_stat_too_complex_detectable)
        else:
            html += get_html_too_complex(global_stat_too_complex_detected)


        # if simulation ?
        if mfasta:
            eval_def = dict()
            eval_def["TP"] = "TP = Detected introns &#8745; Features"
            eval_def["FP"] = "FP = Detected introns &#8713; Features"
            eval_def["FN"] = "FN = Features &#8713; Detected introns"
            eval_def["Se"] = "Sensibility (Se) = TP / (TP+FN)"
            eval_def["Sp"] = "Specificity (Sp) = TP / (TP+FP)"
            eval_def["F1"] = "F1 score = (2*Se*Sp) / (Se+Sp)"
            
            eval_stat = dict()
            eval_stat["0Number of detected introns"] = global_stat_detected_introns["0Number"] 
            eval_stat["1Number of features"] = global_stat["6Number of features in GTF"] 
            eval_stat["2"+eval_def["TP"]] = 0
            eval_stat["3"+eval_def["FP"]] = 0
            eval_stat["4"+eval_def["FN"]] = 0
            eval_stat["5"+eval_def["Se"]] = 0
            eval_stat["6"+eval_def["Sp"]] = 0
            eval_stat["7"+eval_def["F1"]] = 0
            for index, row in df_candidat.iterrows():
                ok = len(df_features.loc[lambda df :
                            (df['contig'] == row['reference']) &
                            (df['start']  == row['start']) &
                            (df['end']    == row['end'])])  
                if ok == 1:
                    eval_stat["2"+eval_def["TP"]] += 1         
                else:
                    eval_stat["3"+eval_def["FP"]] += 1
            eval_stat["4"+eval_def["FN"]] = eval_stat["1Number of features"] - eval_stat["2"+eval_def["TP"]]
            deno = eval_stat["2"+eval_def["TP"]] + eval_stat["4"+eval_def["FN"]]
            if deno :
                eval_stat["5"+eval_def["Se"]] = eval_stat["2"+eval_def["TP"]] / deno * 100
            else :
                eval_stat["5"+eval_def["Se"]] = "NaN"
            deno = (eval_stat["2"+eval_def["TP"]] + eval_stat["3"+eval_def["FP"]])
            if deno :
                eval_stat["6"+eval_def["Sp"]] = eval_stat["2"+eval_def["TP"]] / deno * 100
            else :
                eval_stat["6"+eval_def["Sp"]] = "NaN"
            if eval_stat["5"+eval_def["Se"]] != "NaN" and eval_stat["6"+eval_def["Sp"]] != "NaN" :
                deno = eval_stat["5"+eval_def["Se"]] + eval_stat["6"+eval_def["Sp"]]
                if deno :
                    eval_stat["7"+eval_def["F1"]] = round(2*(eval_stat["5"+eval_def["Se"]]*eval_stat["6"+eval_def["Sp"]]) / deno, 2) 
                else :
                    eval_stat["7"+eval_def["F1"]] = "NaN"
            else :
                eval_stat["7"+eval_def["F1"]] = "NaN"
            if isinstance(eval_stat["5"+eval_def["Se"]], float) :
                eval_stat["5"+eval_def["Se"]] = round(eval_stat["5"+eval_def["Se"]], 2)
            if isinstance(eval_stat["6"+eval_def["Sp"]], float) :
                eval_stat["6"+eval_def["Sp"]] = round(eval_stat["6"+eval_def["Sp"]], 2)

            eval_f_stat = dict()
            eval_f_stat["0Number of detected introns"] = global_stat_f_detected_introns["0" + definitions['PASS']]
            eval_f_stat["1Number of features"] = global_stat_detectable_features["0" + definitions['PASS']]
            eval_f_stat["2"+eval_def["TP"]] = 0
            eval_f_stat["3"+eval_def["FP"]] = 0
            eval_f_stat["4"+eval_def["FN"]] = 0
            eval_f_stat["5"+eval_def["Se"]] = 0
            eval_f_stat["6"+eval_def["Sp"]] = 0
            eval_f_stat["7"+eval_def["F1"]] = 0
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
                        eval_f_stat["2"+eval_def["TP"]] += 1         
                    else:
                        eval_f_stat["3"+eval_def["FP"]] += 1
            eval_f_stat["4"+eval_def["FN"]] = eval_f_stat["1Number of features"] - eval_f_stat["2"+eval_def["TP"]]
            deno = eval_f_stat["2"+eval_def["TP"]] + eval_f_stat["4"+eval_def["FN"]]
            if deno :
                eval_f_stat["5"+eval_def["Se"]] = eval_f_stat["2"+eval_def["TP"]] / deno * 100
            else :
                eval_f_stat["5"+eval_def["Se"]] = "NaN"
            deno = eval_f_stat["2"+eval_def["TP"]] + eval_f_stat["3"+eval_def["FP"]]
            if deno :
                eval_f_stat["6"+eval_def["Sp"]] = eval_f_stat["2"+eval_def["TP"]] / deno * 100
            else :
                eval_f_stat["6"+eval_def["Sp"]] = "NaN"
            if eval_f_stat["5"+eval_def["Se"]] != "NaN" and eval_f_stat["6"+eval_def["Sp"]] != "NaN" :
                deno = eval_f_stat["5"+eval_def["Se"]] + eval_f_stat["6"+eval_def["Sp"]]
                if deno :
                    eval_f_stat["7"+eval_def["F1"]] = round(2*(eval_f_stat["5"+eval_def["Se"]]*eval_f_stat["6"+eval_def["Sp"]]) / deno, 2) 
                else :
                    eval_f_stat["7"+eval_def["F1"]] = "NaN"
            else :
                eval_f_stat["7"+eval_def["F1"]] = "NaN"
            if isinstance(eval_f_stat["5"+eval_def["Se"]], float) :
                eval_f_stat["5"+eval_def["Se"]] = round(eval_f_stat["5"+eval_def["Se"]], 2)
            if isinstance(eval_f_stat["6"+eval_def["Sp"]], float) :
                eval_f_stat["6"+eval_def["Sp"]] = round(eval_f_stat["6"+eval_def["Sp"]], 2)

            html += get_html_eval(eval_stat, eval_f_stat)

    # GLOSSARY
    html += get_html_glossary()

    # FOOTER
    html += "<br/>"
    html += get_html_footer()

    with open(output_file, "w") as f:
        f.write(html)

    print('DATAFRAMES:\n\n')
    print('fasta head :' , df_fasta.head(5), '\n\n')
    # SARAH / CAS REAL
    if mfasta:
        print('mfasta head :' , df_mfasta.head(5), '\n\n')    
    # print('mfasta head :' , df_mfasta.head(5), '\n\n')    
    # SARAH / CAS REAL
    print('library head :' , df_library.head(5), '\n\n')
    print('features head :', df_features.head(5), '\n\n')
    if split:
        print('df_split', df_split, '\n\n')
    if candidat:
        print('df_candidat', df_candidat, '\n\n')

    # SARAH : main stats in a json file
    # https://stackabuse.com/reading-and-writing-json-to-a-file-in-python/  
    genowebReportPath='http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA/report_'+prefix+'_simulation.html'  
    data = {}
    data[output_path] = []
    data[output_path].append({
        'Prefix'                 : prefix,
        'Output path'            : output,
        'Nb seq'                 : global_stat["0Contig FASTA - Number of sequences"],
        'Mean seq. length'       : global_stat["1Contig FASTA - Mean sequence length"],
        'Nb contigs'             : global_stat["2Contig FASTA with feature(s) - Number of sequences"],
        'Mean contigs length'    : global_stat["3Contig FASTA with feature(s) - Mean sequence length"],
        'Nb reads'               : global_stat_fastq["0Number of reads"],
        'Nb prop reads'          : global_stat_flagstat["4Number of properly paired reads"],
        'Nb secondary alignments': secondary,
        'Nb splits'              : global_stat_split["0Number of reads overlapping introns"],
        'Mean len of introns'    : global_stat_split["1Mean length of introns"],
        'Nb detected'            : global_stat_f_detected_introns["0" + definitions['PASS']],
        'Nb detectable'          : global_stat_detectable_features["0" + definitions['PASS']],
        'len detected'           : global_stat_f_detected_introns["6Mean length"],
        'len detectable'         : global_stat_f_detected_introns["6Mean length"],
        'dep detected'           : global_stat_f_detected_introns["9Mean depth"],
        'dep detectable'         : global_stat_detectable_features["9Mean depth"],
        'TP'                     : eval_f_stat["2"+eval_def["TP"]],
        'FN'                     : eval_f_stat["4"+eval_def["FN"]],
        'FP'                     : eval_f_stat["3"+eval_def["FP"]],
        'F1score'                : eval_f_stat["7"+eval_def["F1"]],
        'Se'                     : eval_f_stat["5"+eval_def["Se"]],
        'Sp'                     : eval_f_stat["6"+eval_def["Sp"]],
        'Report'                 : genowebReportPath
    })
    

    # Output path filename comparison json
    output_file = output_path + "_comparison_simulation.json"
    if not force:
        try :
            if os.path.exists(output_file):
                raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file already exists.\n')
            exit(1)

    with open(output_file, 'w') as outfile:
        json.dump(data, outfile)

    # Opening JSON file and loading the data into the variable data 
    with open(output_file) as json_file: 
        data = json.load(json_file) 
    
    stats_data = data[output_path]
    
    # now we will open a file for writing (w) and synthese file without overwriting (a)
    #data_file = open(output_file+'.csv', 'w') 
    #data_file_synthese = open('/home/smaman/Documents/SYNTHESE.csv', 'a') 
    
    # create the csv writer object 
    #csv_writer = csv.writer(data_file) 
    #csv_writer_synthese = csv.writer(data_file_synthese)

    # Counter variable used for writing headers to the CSV file 
    # count = 0
    
    # for i in stats_data: 
    #     if count == 0: 
    
    #         # Writing headers of CSV file 
    #         header = i.keys() 
    #         csv_writer.writerow(header)
    #         count += 1
        
    #     # Write headers in synthesis file only if this file is not empty
    #     if os.stat('/home/smaman/Documents/SYNTHESE.csv').st_size == 0:
    #         header = i.keys()
    #         csv_writer_synthese.writerow(header)

    #     # Writing data of CSV file 
    #     csv_writer.writerow(i.values()) 
    #     csv_writer_synthese.writerow(i.values()) 
    

    # data_file.close()
    # data_file_synthese.close()


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
    parser.add_argument('-m','--modified-fasta', type=argparse.FileType('r'), required=True, dest='mfasta', help="Path to the modified FASTA file.")
    parser.add_argument('-g','--gtf', type=argparse.FileType('r'), required=True, dest='gtf', help="GTF filename which contains the genome annotation.")
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