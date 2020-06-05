#!/usr/bin/env python3

import os
import configparser
import numpy as np  # For Split read signal analysis
import pandas as pd
import argparse
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import pickle
import glob
import time
import configparser # To parse parameters file
import concurrent.futures as prl # For Split read signal analysis
from itertools import repeat     # For Split read signal analysis
# Import all functions from internal modules
from simulation2HTMLparse import *
from simulation2HTMLtags import *
from simulation2HTMLplots import *

config = configparser.RawConfigParser() # On créé un nouvel objet "config"
config.read('parameters') # On lit le fichier de paramètres
# Récupération des parametres dans des variables
fasta        = config.get('MANDATORY','fasta')
mfasta       = config.get('MANDATORY','mfasta')
gtf          = config.get('MANDATORY','gtf')
r1           = config.get('MANDATORY','R1')
output       = config.get('MANDATORY','output')
threads      = config.get('MANDATORY','t')
r2           = config.get('OPTIONNAL','R2')
flagstat     = config.get('OPTIONNAL','flagstat')
ranks        = config.get('OPTIONNAL','ranks')
assemblathon = config.get('OPTIONNAL','assemblathon')
candidat     = config.get('OPTIONNAL','candidat')
split        = config.get('OPTIONNAL','split')
prefix       = config.get('OPTIONNAL','prefix')
force        = config.get('OPTIONNAL','force')

#source activate ISeeker_environment;
#cd scripts/; 
# python3 simulation2HTML.py -m /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs-modified.gtf -o /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML -p test1 -F  -1 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R2.fastq.gz --flagstat /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt -r /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_ranks.txt -s  /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_split_alignments.txt  --assemblathon /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6
# scp  /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML/*.html smaman@genologin.toulouse.inra.fr:/save/smaman/public_html/intronSeeker/.
# See result : http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/report_test1_simulation.html

############
# SUB MAIN #
############
def simulationReport(   fasta:str, mfasta:str, gtf:str, r1:str, r2:str, ranks:str,
                        assemblathon:str, flagstat:str, candidat:str, split:str,
                        output:str, prefix:str, force:bool, threads:int ) :
    
    #Return the value (in fractional seconds) of the sum of the system and user CPU time of the current process. 
    #It does not include time elapsed during sleep. 
    #It is process-wide by definition. 
    #The reference point of the returned value is undefined, 
    #so that only the difference between the results of consecutive calls is valid.
    tmps1=time.process_time()
    
    #Return the value (in fractional seconds) of a performance counter
    #Only the difference between the results of consecutive calls is valid.
    tmps1c=time.perf_counter() 

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

    print("Output dir and report html")        

	### MEMO 
	# fasta  = sequences used to generate reads
	# mfasta = sequences used to align reads
	
	### Build pandas ###
    # df_fasta    : pandas.DataFrame where each line is a seq description from FASTA file
    # df_mfasta   : pandas.DataFrame where each line is a modified seq description from modified FASTA file
    # df_library  : pandas.DataFrame where each line is a read description from R1 (& R2) fastq file(s)
    # df_features : pandas.DataFrame where each line is a simulated features description 
    # df_candidat : pandas.DataFrame where each line is a candidat description
    df_fasta  = parse_fasta(fasta, False)
    df_mfasta = parse_fasta(mfasta, True)

    if r2 :
        df_library = parse_library(r1, r2)
    else :
        df_library = parse_library(r1)

    df_features = parse_gtf(gtf)

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

    print("Build pandas dataframes")  
    print("CPU time = %f" %(time.process_time()-tmps1))
    print("Performance counter = %f\n" %(time.perf_counter()-tmps1c))
    tmps2=time.process_time()
    tmps2c=time.perf_counter() 

    # HEADER
    html = get_html_header()
    
    inputfiles = [
        "Contig FASTA#" + os.path.basename(fasta),
        "Contig FASTA with feature(s)#" + os.path.basename(mfasta),
        "GTF of contig with feature(s)#" + os.path.basename(gtf),
        "Read1 FASTQ#" + os.path.basename(r1)
    ]
    if r2:
        inputfiles.append("Read2 FASTQ#" + os.path.basename(r2))
    ranks_file=""    
    if ranks:
        ranks_file=ranks
        inputfiles.append("Ranks#" + os.path.basename(ranks))
    flagstat_file = ""
    if flagstat :
        flagstat_file = flagstat
        inputfiles.append("Flagstat#" + os.path.basename(flagstat_file))
    assemblathon_file=""    
    if assemblathon:
        assemblathon_file=assemblathon
        inputfiles.append("Assemblathon#" + os.path.basename(assemblathon))
    candidat_file=""    
    if candidat:
        candidat_file=candidat
        inputfiles.append("Candidat#" + os.path.basename(candidat))

    html += get_html_body1(flagstat_file, candidat_file, assemblathon_file)

   
    # INPUT FILES
    html += get_html_inputfiles(inputfiles)
    print("List input files")

    # SEQUENCE STAT
    # Global stat
    nb_distinct_features, nb_ctg_by_feature, ctg_descr = stat_from_gtf(gtf)
    global_stat = dict()
    global_stat["0Contig FASTA - Number of seq."]            = df_fasta.shape[0]
    global_stat["1Contig FASTA - Mean seq. length"]          = int(df_fasta['length'].mean())
    global_stat["2Contig FASTA with feature(s) - Number of seq."]   = df_mfasta.shape[0]
    global_stat["3Contig FASTA with feature(s) - Mean seq. length"] = int(df_mfasta['length'].mean())
    global_stat["4Number of modified sequences"]           = df_mfasta.loc[df_mfasta.index.str.contains(".modif")].shape[0]
    global_stat["5Number of distinct features in GTF"]     = df_features.feature.value_counts().shape[0]
    global_stat["6Number of features in GTF"]              = df_features.shape[0]
    c = 7
    for k, v in (df_features.feature.value_counts()).items() :
        global_stat[str(c)+k] = v
        c+=1

    html += get_html_seq_descr(global_stat, nb_ctg_by_feature, ctg_descr, gtf, df_features['pos_on_contig'], df_fasta, df_mfasta)
    print("Global statistics")
    print("CPU time = %f" %(time.process_time()-tmps2))
    print("Performance counter = %f\n" %(time.perf_counter()-tmps2c))
    tmps3=time.process_time()
    tmps3c=time.perf_counter() 

    # READS STAT
    # Global stat
    global_stat_fastq = dict()
    global_stat_fastq["0Number of reads"] = df_library['contig'].count()
    global_stat_fastq["1Mean coverage"] = 0
    for i,row in df_library.iterrows():
        global_stat_fastq["1Mean coverage"] += row['end'] - row['start'] + 1
    global_stat_fastq["1Mean coverage"] /= (global_stat["1Contig FASTA - Mean seq. length"] * global_stat["0Contig FASTA - Number of seq."])
    global_stat_fastq["1Mean coverage"] = round(global_stat_fastq["1Mean coverage"], 2)
    global_stat_fastq["2Min reads length"] = df_features['length'].min()
    global_stat_fastq["3Max reads length"] = df_features['length'].max()
    global_stat_fastq["4Mean reads length"] = round(df_features['length'].mean())
      
    html += get_html_reads_descr(global_stat_fastq)
    print("Reads statistics")
    print("CPU time = %f" %(time.process_time()-tmps3))
    print("Performance counter = %f\n" %(time.perf_counter()-tmps3c))
    tmps4=time.process_time()
    tmps4c=time.perf_counter() 

    # ABUNDANCE number of reads by contig
    # Build a dataframe with:
    #   ctg
    #   abund_perc => (number of read on this contig / number of reads) * 100
    #   requested  => if grinder ranks output file is given ...
    #   norm       => (((number of read on this contig / contig len) * mean len of all contigs) / number of reads) * 100
    df_tmp = pd.DataFrame((df_library.groupby('contig').size()/len(df_library))*100, columns = ['abund_perc'])
    df_fasta = df_fasta.assign(real=df_tmp.abund_perc.values)
    if ranks:
        df_tmp = parse_rank_file(ranks_file)
        df_fasta = df_fasta.assign(rank=df_tmp['rank'].values)
        df_fasta = df_fasta.assign(waiting=df_tmp.rel_abund_perc.values)
    df_tmp = pd.DataFrame(
        (((df_library.groupby('contig').size()/df_fasta['length'])*(df_fasta['length'].mean()))/df_library.shape[0])*100,
        columns = ['norm'])
    df_fasta = df_fasta.assign(norm=df_tmp['norm'].values)
    del df_tmp
    html += get_html_abundance(df_fasta)
    print("Abundance")
    print("CPU time = %f" %(time.process_time()-tmps4))
    print("Performance counter = %f\n" %(time.perf_counter()-tmps4c))
    tmps5=time.process_time()
    tmps5c=time.perf_counter() 

    ## ALIGNMENT STATS
    if flagstat:
        nbreads, mapped, paired, proper, secondary, singletons=parse_flagstat(flagstat_file)

        global_stat_flagstat = dict()
        global_stat_flagstat["0Number of mapped"] = mapped
        global_stat_flagstat["1Number of properly paired reads"] = proper
        global_stat_flagstat["2Number of unmapped reads (nb QC passed - nb mapped)"] = int(nbreads) - int(mapped)
        global_stat_flagstat["3Secondary"] = secondary
        global_stat_flagstat["4Singletons"] = singletons

        html += get_html_flagstat_descr(global_stat_flagstat)
   
    html += get_html_results()
    print("Mapping statistics")
    print("CPU time = %f" %(time.process_time()-tmps5))
    print("Performance counter = %f\n" %(time.perf_counter()-tmps5c))
    tmps6=time.process_time()
    tmps6c=time.perf_counter() 


    ## SPLITREADSEARCH STAT
    if split:
        df_split=parse_split(split)   
        global_stat_split = dict()
        global_stat_split["0Number of reads overlapping potential retained introns"]= df_split.shape[0]
        global_stat_split["1Mean length of potential retained introns"]= df_split['split_length'].mean()

        nbCanonic = 0
        nbOtherJunctions = 0
        df_split.sort_values(by=['split_borders'])
        c = 2
        n = 0
        for k, v in (df_split['split_borders'].value_counts()).items() :
            if k == "GT_AG" or k == "CT_AC":
                nbCanonic += v
            elif n < 11:
                global_stat_split[str(c)+"Junction "+k] = v
                n += 1
                c += 1
            else:
                nbOtherJunctions += v
        global_stat_split[str(c)+"Other junctions"] = nbOtherJunctions
        global_stat_split[str(c+1)+"Canonical junction (GT_AG or CT_AC)"] = nbCanonic
        
        html += get_html_split_descr(global_stat_split)   
        print("Intron reads") 
        print("CPU time = %f" %(time.process_time()-tmps6))
        print("Performance counter = %f\n" %(time.perf_counter()-tmps6c))
        tmps7=time.process_time()
        tmps7c=time.perf_counter() 


    # Candidat statistics - detected introns
    if candidat:
        df_candidat, mindepth, maxlen = parse_candidat(candidat)

         # Definition dict
        definitions = dict()
        definitions['DP']   = "Number of detected introns depth <= min depth ("+ str(mindepth) +")"
        definitions['LEN']  = "Number of detected introns length >= max len ("+ str(maxlen) +")"
        definitions['SS']   = "Number of features without canonical borders (neither CT_AC nor GT_AG)"
        definitions['PASS'] = "Number of features with canonical borders (CT_AC or GT_AG) and ...."
        
        global_stat_detected_introns = dict()
        global_stat_detected_introns["0Number"] = df_candidat.shape[0]
        global_stat_detected_introns["1Mean length"] = 0
        for i,row in df_candidat.iterrows():
            global_stat_detected_introns["1Mean length"] += row['end'] - row['start'] + 1
        global_stat_detected_introns["1Mean length"] /= (global_stat_detected_introns["0Number"])
        global_stat_detected_introns["1Mean length"] = round(global_stat_detected_introns["1Mean length"], 2)
        global_stat_detected_introns["2Mean depth"]= round(df_candidat['depth'].mean(), 2)


        global_stat_filtred_detected_introns = dict()
        c = 0
        for k, v in (df_candidat['filter'].value_counts()).items() :
            for key, value in definitions.items():
                if k == key:
                    global_stat_filtred_detected_introns[str(c)+k+" ("+value+")"] = v
                    c+=1


        # Comparison between candidats and features from GTF file
        nbTotCandidatsIncludingFeatures, nbSameStartEnd, nbLen, minimumDepth, nonCanonical=candidatsVsFeatures(df_candidat, df_features, mindepth, maxlen)

        global_stat_candidat_vs_gtf = dict()
        global_stat_candidat_vs_gtf["0Number of detected introns"] = global_stat_detected_introns["0Number"]
        c = 1
        for k, v in (df_features.feature.value_counts()).items() :
            global_stat_candidat_vs_gtf[str(c)+'Number of '+k+ ' in GTF'] = v
            c+=1
        global_stat_candidat_vs_gtf[str(c+1)+"Number of detected introns corresponding features (Overlaps)"] = nbSameStartEnd
        global_stat_candidat_vs_gtf[str(c+2)+"Detected introns not found in GTF"]   = global_stat_detected_introns["0Number"]- nbTotCandidatsIncludingFeatures
        

        html += get_html_candidat_descr(global_stat_detected_introns, global_stat_filtred_detected_introns, df_candidat)
        print("Detected introns statistics and histogram")
        print("CPU time = %f" %(time.process_time()-tmps7))
        print("Performance counter = %f\n" %(time.perf_counter()-tmps7c))
        tmps8=time.process_time()
        tmps8c=time.perf_counter() 

    # Assemblathon files
    if assemblathon:
        df_assemblathon_all = parse_assemblathon(assemblathon_file, "title")
        global_stat_assemblathon = dict()
        global_stat_assemblathon["0Number of contigs"]         = df_assemblathon_all.iloc[0,0]
        global_stat_assemblathon["1Total size of contigs"]     = df_assemblathon_all.iloc[1,0]
        global_stat_assemblathon["2Longest contig"]            = df_assemblathon_all.iloc[2,0]
        global_stat_assemblathon["3Shortest contiged"]         = df_assemblathon_all.iloc[3,0]
        nbLongContigs=re.sub(r'([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%)', r" ", df_assemblathon_all.iloc[4,0])
        global_stat_assemblathon["4Number of contigs > 1K nt"] = nbLongContigs
        global_stat_assemblathon["5N50 contig length"]         = df_assemblathon_all.iloc[5,0]
        global_stat_assemblathon["6L50 contig count"]          = df_assemblathon_all.iloc[6,0]
        html += get_html_assemblathon_descr(global_stat_assemblathon)    
        print("Assemblathon statistics")
        print("CPU time = %f" %(time.process_time()-tmps8))
        print("Performance counter = %f\n" %(time.perf_counter()-tmps8c))
        tmps9=time.process_time()
        tmps9c=time.perf_counter() 

    # Precision, recall and F1 score
    # TP is the number of detectable and found features (int value)
    # TN is the number of detectable and not found features (int value)
    # FP is the number of undetectable and found features (int value)
    # FN is the number of undetectable and not found features (int value)
    
    # https://fr.wikipedia.org/wiki/Pr%C3%A9cision_et_rappel"
    ##  les "features" qui ont assez de lectures les couvrant pour être trouvées. Il n'y que celles-ci qui pourront être vues comme T (True).
    # TP = nombre de "features" détectables et trouvées (assez de profondeur et bonnes bornes). 
    #TP = nbOverlap - nonCanonical - minimumDepth
    TP = nbTotCandidatsIncludingFeatures
    # TN = nombre de "features" détectables non couvertes par des lectures et/ou couvertes par un nombre de lectures après alignement sont sous le seuil de profondeur (alors que dans la simulation c'était le cas)
    # PROFONDEUR (DEP)
    #TN =  nbOverlap - df_features.shape[0]
    TN =  minimumDepth
    # FP = nombre de zones hors "feature" ou "features" indétectables avec une couverture insuffisante
    # COUVERTURE (LEN)
    #FP = nbOverlap - df_candidat.shape[0]
    FP = nbSameStartEnd
    # FN = nombre de zones hors "feature" ou "features" indétectables trouvées (passant le seuil de profondeur)
    # nb candidat - nb features
    FN = df_candidat.shape[0] - df_features.shape[0]

    global_stat_precision= dict()
    precision = TP/(FP+TP)
    global_stat_precision["0Precision (between 0 - 1)"]= precision
    recall = TP/(TN+TP)
    global_stat_precision["1Recall or sensitivity (0.0 for no recall, 1.0 for full or perfect recall)"] = recall
    #global_stat_precision["2F1 score (1 for a perfect model, 0 for a failed model)"]  = 2*((global_stat_precision["0Precision"]*global_stat_precision["1Recall"])/(global_stat_precision["0Precision"]+global_stat_precision["1Recall"]))
    global_stat_precision["2F1 score (1 for a perfect model, 0 for a failed model)"]  = 2*((precision*recall)/(precision+recall))

    html += get_html_precision(global_stat_precision, TP, TN, FP, FN, global_stat_candidat_vs_gtf)

    print("Detectability statistics")
    print("CPU time = %f" %(time.process_time()-tmps9))
    print("Performance counter = %f\n" %(time.perf_counter()-tmps9c))
    tmps10=time.process_time()
    tmps10c=time.perf_counter() 

    # GLOSSARY
    html += get_html_glossary()

    # FOOTER
    html += "<br/>"
    html += get_html_footer()

    with open(output_file, "w") as f:
        f.write(html)

    print('DATAFRAMES:\n\n')
    print('fasta head :' , df_fasta.head(5), '\n\n')
    print('mfasta head :' , df_mfasta.head(5), '\n\n')
    print('library head :' , df_library.head(5), '\n\n')
    print('features head :', df_features.head(5), '\n\n')
    if split:
        print('df_split', df_split, '\n\n')
    if candidat:
        print('df_candidat', df_candidat, '\n\n')
    if assemblathon:
        print('df_assemblathon', df_assemblathon_all, '\n\n')


if __name__ == '__main__' :
    # parser = argparse.ArgumentParser(add_help=False)
    # parser.add_argument('-f','--fasta', type=argparse.FileType('r'), required=True, dest='fasta')
    # parser.add_argument('-m','--modifiedfasta', type=argparse.FileType('r'), required=True, dest='mfasta')
    # parser.add_argument('-g','--gtf', type=argparse.FileType('r'), required=True, dest='gtf')
    # parser.add_argument('-1','--R1', type=argparse.FileType('r'), required=True, dest='r1')
    # parser.add_argument('-2','--R2', type=argparse.FileType('r'), required=False, dest='r2')
    # parser.add_argument('--flagstat', type=argparse.FileType('r'), required=False, dest='flagstat')
    # parser.add_argument('--assemblathon', type=argparse.FileType('r'), required=False, dest='assemblathon') 
    # parser.add_argument('-r','--ranksfile', type=argparse.FileType('r'), required=False, dest='ranks')
    # parser.add_argument('-c','--candidat', type=argparse.FileType('r'), required=False, dest='candidat')
    # parser.add_argument('-s','--split', type=argparse.FileType('r'), required=False, dest='split')
    # parser.add_argument('-o','--output', type=str, required=True, dest='output')
    # parser.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    # parser.add_argument('-t','--threads', type=int, default=1, required=False, dest='threads')
    # parser.add_argument('-F', '--force', action='store_true', default=False, dest='force')

    # args = vars(parser.parse_args())
    
    # simulationReport(**args)
    simulationReport(fasta, mfasta, gtf, r1, r2, ranks,
                        assemblathon, flagstat, candidat, split,
                        output, prefix, force, threads)