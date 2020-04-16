#!/usr/bin/env python3

import os
import configparser
import numpy as np  # For Split read signal analysis
import pandas as pd
import argparse
import pysam   # To generate a dataframe from a BAM : pysam and pickle
import pickle
import glob
import concurrent.futures as prl # For Split read signal analysis
from itertools import repeat     # For Split read signal analysis
# Import all functions from internal modules
from simulation2HTMLparse import *
from simulation2HTMLtags import *
from simulation2HTMLplots import *




#step 1  full random simulation : intronSeeker fullRandomSimulation -r -o FRS/ 
#source activate ISeeker_environment;
#cd scripts/; 
# python3 simulation2HTML.py -m /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML -p test1 -F  -1 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R2.fastq.gz --flagstat /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt -r /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_ranks.txt -b /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.bam --assemblathon /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6
# scp  /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML/*.html smaman@genologin.toulouse.inra.fr:/save/smaman/public_html/intronSeeker/.
# See result : http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/report_test1_simulation.html

############
# SUB MAIN #
############
def simulationReport(   fasta:str, mfasta:str, gtf:str, r1:str, r2:str, ranks:str,
                        assemblathon:str, flagstat:str, bam:str, candidat:str,
                        output:str, prefix:str, force:bool, threads:int ) :
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
    
    inputfiles = [
        "Contig FASTA#" + os.path.basename(fasta.name),
        "Contig FASTA with feature(s)#" + os.path.basename(mfasta.name),
        "GTF of modified sequences#" + os.path.basename(gtf.name),
        "Read1 FASTQ#" + os.path.basename(r1.name)
    ]
    if r2:
        inputfiles.append("Read2 FASTQ#" + os.path.basename(r2.name))
    ranks_file=""    
    if ranks:
        ranks_file=ranks.name
        inputfiles.append("Ranks#" + os.path.basename(ranks.name))
    assemblathon_file=""    
    if bam:
        inputfiles.append("Bam#" + os.path.basename(bam.name))
    flagstat_file = ""
    if flagstat :
        flagstat_file = flagstat.name
        inputfiles.append("Flagstat#" + os.path.basename(flagstat_file))
    if assemblathon:
        assemblathon_file=assemblathon.name
        inputfiles.append("Assemblathon#" + os.path.basename(assemblathon.name))
    candidat_file=""    
    if candidat:
        candidat_file=candidat.name
        inputfiles.append("Candidat#" + os.path.basename(candidat.name))

    html += get_html_body1(flagstat_file, candidat_file, assemblathon_file)

   
    # INPUT FILES
    html += get_html_inputfiles(inputfiles)

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

    global_stat_assemblathon = dict()
    #Compare assemblathon files
    if assemblathon:
        df_assemblathon_all = parse_assemblathon(assemblathon_file, "title")
        print('df_assemblathon_all', df_assemblathon_all)
        #global_stat_assemblathon = dict()
        print(df_assemblathon_all.shape[0])#lignes 
        print(df_assemblathon_all.shape[1])#colonnes 
        global_stat_assemblathon["0Number of contigs"]         = df_assemblathon_all.iloc[0,0]
        global_stat_assemblathon["1Total size of contigs"]     = df_assemblathon_all.iloc[1,0]
        global_stat_assemblathon["2Longest contig"]            = df_assemblathon_all.iloc[2,0]
        global_stat_assemblathon["3Shortest contiged"]         = df_assemblathon_all.iloc[3,0]
        nbLongContigs=re.sub(r'([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%)', r" ", df_assemblathon_all.iloc[4,0])
        global_stat_assemblathon["4Number of contigs > 1K nt"] = nbLongContigs
        global_stat_assemblathon["5N50 contig length"]         = df_assemblathon_all.iloc[5,0]
        global_stat_assemblathon["6L50 contig count"]          = df_assemblathon_all.iloc[6,0]
       # html += get_html_assemblathon_descr(global_stat_assemblathon)

    html += get_html_seq_descr(global_stat, nb_ctg_by_feature, ctg_descr, gtf.name, df_features['pos_on_contig'], df_fasta, df_mfasta, global_stat_assemblathon)
  
    '''
    #https://stackoverflow.com/questions/45759966/counting-unique-values-in-a-column-in-pandas-dataframe-like-in-qlik
    print('6Number of features in GTF:', df_features['contig'].nunique())
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

    # READS STAT
    # Global stat
    global_stat_fastq = dict()
    global_stat_fastq["0Number of reads"] = df_library['contig'].count()
    global_stat_fastq["1Mean coverage"] = 0
    for i,row in df_library.iterrows():
        global_stat_fastq["1Mean coverage"] += row['end'] - row['start'] + 1
    global_stat_fastq["1Mean coverage"] /= (global_stat["1FASTA file - Mean seq. length"] * global_stat["0FASTA file - Number of seq."])
    global_stat_fastq["1Mean coverage"] = round(global_stat_fastq["1Mean coverage"], 2)
    #to remove ?
    global_stat_fastq["2Min reads length"] = df_features['length'].min()
    global_stat_fastq["3Max reads length"] = df_features['length'].max()
    global_stat_fastq["4Mean reads length"] = round(df_features['length'].mean())
      
    html += get_html_reads_descr(global_stat_fastq)

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

    ## ALIGNMENT STATS
    if flagstat:
        df_flag_all=parse_flagstat(flagstat_file, len(df_library),"Parse flagstat")
   
        #print(pd.concat([df_flag_all_hisat,df_flag_all_star],axis=1,sort=False).fillna(0))
        #df_flag_all=pd.concat([df_flag_all_hisat,df_flag_all_star],axis=1,sort=False).fillna(0)
        global_stat_flagstat = dict()
        global_stat_flagstat["0Total counts of reads to map"] = str(len(df_library))
        global_stat_flagstat["1Total count"]     = df_flag_all.iloc[0,0]
        global_stat_flagstat["2Secondary"]       = df_flag_all.iloc[1,0]
        global_stat_flagstat["3Mapped"]          = df_flag_all.iloc[2,0]
        global_stat_flagstat["4Properly paired"] = df_flag_all.iloc[3,0]
        global_stat_flagstat["5Singletons"]      = df_flag_all.iloc[4,0]
        html += get_html_flagstat_descr(global_stat_flagstat, flagstat.name, df_flag_all)
        #https://plotly.com/python/sunburst-charts/
    
    #MAPPING
    #Counting table and barplots of mapped covering reads' main characteristics
    if bam:
        bam_file=bam.name
        # Add three columns on df_library : 
        # 1 - one for the intron covering reads (True/False) :   covering  : True/False 
        # 2 - another for the covered intron id (if True)  intron : SEQUENCE685.modif|556|897
        # 3 - and the last for the intron insertion position in read (if True - in term of read length percentage) : pos_on_read : 98.019802
        df_cov=prlz_process_intron(df_features, df_library)
        df_library = df_library.join(df_cov,lsuffix='',rsuffix='_cov').loc[:,df_cov.columns]
        df_library.loc[lambda df : df.covering != True, "covering"] = False
        #mapping_bam=pd.read_pickle(process_bam(parse_BAM(bam_file), df_mfasta, df_features, df_library))
        #html += get_html_split(mapping_bam)
        #html += get_html_mapping_descr(parse_BAM(bam_file))


    ## SPLITREADSEARCH STAT
    if candidat:
        df_candidat = parse_candidat(candidat.name)
        print('candidat head :' , df_candidat)
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
        html += get_html_candidat_descr(global_stat_candidat, df_candidat[df_candidat['filter'] == 'PASS'])
             
    # GLOSSARY
    html += get_html_glossary()

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
    parser.add_argument('--flagstat', type=argparse.FileType('r'), required=False, dest='flagstat')
    parser.add_argument('-b','--bam', type=argparse.FileType('r'), required=False, dest='bam')
    parser.add_argument('--assemblathon', type=argparse.FileType('r'), required=False, dest='assemblathon') 
    parser.add_argument('-r','--ranksfile', type=argparse.FileType('r'), required=False, dest='ranks')
    parser.add_argument('-c','--candidat', type=argparse.FileType('r'), required=False, dest='candidat')
    parser.add_argument('-o','--output', type=str, required=True, dest='output')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="", dest='prefix')
    parser.add_argument('-t','--threads', type=int, default=1, required=False, dest='threads')
    parser.add_argument('-F', '--force', action='store_true', default=False, dest='force')

    args = vars(parser.parse_args())
    
    simulationReport(**args)