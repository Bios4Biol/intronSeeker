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
#$ conda deactivate
#module load system/Miniconda3-4.7.10;
#source activate ISeeker_environment;
#cd scripts/; 
#(ISeeker_environment) sigenae@genologin1 /work/project/sigenae/sarah/intronSeeker/scripts $ python3 simulation2HTML.py -m ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -1 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R2.fastq.gz -o HTML -p tests -F  --frs  ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf  -D ../../archives_intronSeeker/TESTS/
#(ISeeker_environment) sigenae@genologin1 /work/project/sigenae/sarah/intronSeeker/scripts $ python3 simulation2HTML.py -m ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o HTML -p tests -F  -1 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 ../../archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R2.fastq.gz
#python3 simulation2HTML.py -m /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/HTML -p TOTO -F  -1 /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sr_R2.fastq.gz -a /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /work/project/sigenae/sarah/archives_intronSeeker/TESTS/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt
#python3 simulation2HTML.py -m /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML -p test1 -F  -1 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R2.fastq.gz -a /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt -r /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_ranks.txt -S /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -H /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.flagstat.txt -bha /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.bam -bhm /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.bam -bsa /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.bam -bsm /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.bam -aia /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_HISAT2/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6
#python3 simulation2HTML.py -m /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs-modified.fa -f /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_contigs.fa -g /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/frs_sample1_modifications.gtf -o /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HTML -p test1 -F  -1 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R1.fastq.gz -2 /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_R2.fastq.gz -a /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /home/Sarah/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt -r /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sr_ranks.txt --workDirSTAR /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/ --workDirHISAT /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/ -bha /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.bam -bhm /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/HISAT2_alignment/hisat2.sort.bam -bsa /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.bam -bsm /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/STAR_alignment/star.sort.bam -aia /home/smaman/Documents/PROJETS/INTRONSEEKER/FRS/CAS-A/sample1/sample1_splicing_event_HISAT2/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6

############
# SUB MAIN #
############
def simulationReport(fasta:str, mfasta:str, gtf:str, r1:str, r2:str, ranksfile:str, allwiintronsassemblathon:str, flagstat:str, workDirSTAR:str, workDirHISAT:str, bamHISATall:str, bamHISATmix: str, bamSTARall: str, bamSTARmix: str, candidat:str, output:str, prefix:str, force:bool, threads:int) :
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
    # With workdirs as parameters
    flgSTARintrons    =[f for f in glob.glob(workDirSTAR + "*flagstat*.txt", recursive=True)]
    #df_flag_all_star  =parse_flagstat(flgSTARintrons[-1], len(df_library),"All with introns - Star")
    df_flag_all_star  =parse_flagstat(flgSTARintrons[-1], len(df_library),workDirSTAR)
    flgHISAT2introns  =[f for f in glob.glob(workDirHISAT + "*flagstat*.txt", recursive=True)]
    #df_flag_all_hisat =parse_flagstat(flgHISAT2introns[-1], len(df_library),"All with introns - Hisat2")
    df_flag_all_hisat =parse_flagstat(flgHISAT2introns[-1], len(df_library),workDirHISAT)
   
    # If 1 parameter = 1 file
    '''
    df_flag_all_hisat = parse_flagstat(flgHISAT2introns.name, len(df_library),"All with introns - Hisat2")
    #df_flag_half_hisat = parse_flagstat(flgSTAR, len(df_library),"Mix-states contigs - Hisat2")
    df_flag_all_star = parse_flagstat(flgSTARintrons.name, len(df_library),"All with introns - STAR")
    #df_flag_half_star = parse_flagstat(flgHISAT2, len(df_library),"Mix-states contigs - STAR")    
    '''

    print(pd.concat([df_flag_all_hisat,df_flag_all_star],axis=1,sort=False).fillna(0))
    df_flag_all=pd.concat([df_flag_all_hisat,df_flag_all_star],axis=1,sort=False).fillna(0)
    global_stat_flagstat_hisat2 = dict()
    global_stat_flagstat_hisat2["0Total counts of reads to map"] = str(len(df_library))
    global_stat_flagstat_hisat2["1Total count"]     = df_flag_all.iloc[0,0]
    global_stat_flagstat_hisat2["2Secondary"]       = df_flag_all.iloc[1,0]
    global_stat_flagstat_hisat2["3Mapped"]          = df_flag_all.iloc[2,0]
    global_stat_flagstat_hisat2["4Properly paired"] = df_flag_all.iloc[3,0]
    global_stat_flagstat_hisat2["5Singletons"]      = df_flag_all.iloc[4,0]
    global_stat_flagstat_star= dict()
    global_stat_flagstat_star["0Total counts of reads to map"] = str(len(df_library))
    global_stat_flagstat_star["1Total count"]       = df_flag_all.iloc[0,1]
    global_stat_flagstat_star["2Secondary"]         = df_flag_all.iloc[1,1]
    global_stat_flagstat_star["3Mapped"]            = df_flag_all.iloc[2,1]
    global_stat_flagstat_star["4Properly paired"]   = df_flag_all.iloc[3,1]
    global_stat_flagstat_star["5Singletons"]        = df_flag_all.iloc[4,1]
    html += get_html_flagstat_descr(global_stat_flagstat_hisat2, global_stat_flagstat_star, flgSTARintrons[-1], flgHISAT2introns[-1])
    html += get_html_flagstat_pie(df_flag_all)
    

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
    # To decrease number of parameters, just give workdir path to STAR and HISAT files
    parser.add_argument('--workDirSTAR', type=str, required=True, dest='workDirSTAR')
    parser.add_argument('--workDirHISAT', type=str, required=True, dest='workDirHISAT')
    #parser.add_argument('-S','--flgSTARintrons', type=argparse.FileType('r'), required=True, dest='flgSTARintrons')
    #parser.add_argument('-s','--flgSTAR', type=argparse.FileType('r'), required=True, dest='flgSTAR')
    #parser.add_argument('-H','--flgHISAT2introns', type=argparse.FileType('r'), required=True, dest='flgHISAT2introns')
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
