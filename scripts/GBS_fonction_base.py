#!/usr/bin/env python3

import configparser
import argparse
import numpy.random as rd
import tempfile as tmp
import subprocess as sp
import pandas as pd
import os
from collections import OrderedDict,defaultdict
from pprint import pprint


def read_gtf(path_to_file):
    """
    Reads an entire GTF file.
    Returns the list of all feature and a set of transcripts ID.
        
    :param path_to_file: GTF file name.
    :return liste: List of list. Each sublist corresponds to a feature of the GTF file.
    :rtype: list
    :return transcript: Set of transcript IDs. Only transcript annotated as "protein_coding" are selected.
    :rtype: set
    
    :Example:
    >>> gtf,transcripts = read_gtf(os.path.abspath(os.path.dirname(__file__))+"/../data/test_GBS.gtf") # doctest: +SKIP
    >>> gtf[0] # doctest: +SKIP
    ['V', 'WormBase', 'gene', '180', '329', '.', '+', '.', {'gene_id': 'WBGene00197333', 'gene_version': '1', 'gene_name': 'cTel3X.2', 'gene_source': 'WormBase', 'gene_biotype': 'ncRNA'}]
    >>> gtf[:3] # doctest: +SKIP
    [['V', 'WormBase', 'gene', '180', '329', '.', '+', '.', {'gene_id': 'WBGene00197333', 'gene_version': '1', 'gene_name': 'cTel3X.2', 'gene_source': 'WormBase', 'gene_biotype': 'ncRNA'}],
     ['V', 'WormBase', 'transcript', '180', '329', '.', '+', '.', {'gene_id': 'WBGene00197333', 'gene_version': '1', 'transcript_id': 'cTel3X.2', 'gene_name': 'cTel3X.2', 'gene_source': 'WormBase', 'gene_biotype': 'ncRNA', 'transcript_name': 'cTel3X.2', 'transcript_source': 'WormBase', 'transcript_biotype': 'ncRNA'}],
     ['V', 'WormBase', 'exon', '180', '329', '.', '+', '.', {'gene_id': 'WBGene00197333', 'gene_version': '1', 'transcript_id': 'cTel3X.2', 'exon_number': '1', 'gene_name': 'cTel3X.2', 'gene_source': 'WormBase', 'gene_biotype': 'ncRNA', 'transcript_name': 'cTel3X.2', 'transcript_source': 'WormBase', 'transcript_biotype': 'ncRNA', 'exon_id': 'cTel3X.2.e1'}]]
    >>> transcripts # doctest: +SKIP
    {'B0348.4e.1', 'K08E3.5b.1', 'K08E3.7a.2', ... , 'B0348.4u.1', 'K08E3.5e.1', 'K08E3.5a.2'}
        
    .. seealso:: annoToData()
    .. warning::
        Developed to read Ensembl GTF files. A file from another source could 
        induce errors particularly if the additional attributes are different of
        these defined by Ensembl
    """
    
    liste = []
    transcripts = set()
    with open(path_to_file) as gtf:
        for ligne in gtf:
            if ligne.startswith("#"):  # We ignore the file's header
                continue
            feature = ligne.rstrip().split("\t")
            feature.append({
                el.split()[0]: el.split()[1].strip('"') 
                for el in feature.pop().rstrip(";").split("; ")
                })
            # We check if it's a transcript and if produce a protein
            if feature[2] == "transcript" and feature[-1]["transcript_biotype"] == "protein_coding":
                # we add the transcript id to the set of transcript ID in order
                # to pickly random among them
                transcripts.add(feature[-1]["transcript_id"])
            liste.append(feature)
        return liste, transcripts


def make_density_law():
    """
    Reads a distribution among several class from the property file 
    and returns a tuple where each class is represented according to distribution.
    If you want to modify this distribution, modify the property file.
    The distibution have to be expressed in term of effectives or percentages (i.e. only int values)
    
    :Example:
    ---Extract of intronStalker.properties file for this example---
    [Density]
    -1=5
    0=60
    1=20
    2=10
    3=5
    ---------------------------------------------------------------
    >>> law = make_density_law() # doctest: +SKIP
    >>> set(law) # doctest: +SKIP
    {0, 1, 2, 3, -1}
    >>> law.count(0) # doctest: +SKIP
    60
    >>> law.count(-1) # doctest: +SKIP
    5
    >>> law.count(1) # doctest: +SKIP
    20
    >>> law.count(2) # doctest: +SKIP
    10
    >>> law.count(3) # doctest: +SKIP
    5
    
    ..seealso: choose_transcript()
    """
    law = []

    config = configparser.RawConfigParser()
    config.read(os.path.abspath(os.path.dirname(__file__))+"/../config/intronStalker.properties")
    for classe in config["Density"]:
        effectif = int(config["Density"][classe])
        law += [int(classe)] * effectif
    return tuple(law)


def choose_transcripts(transcripts, nb_to_choose):
    """
    Randomly picks a given number of transcripts among a list of transcripts, and for each,
    randomly assign a class according to the distribution given by the make_density_law() function.
    Returns a dict object where the keys correspond to transcripts IDs and the values are their class.
    If the precised number of transcripts is greater than the transcripts list's length,
    all the transcripts IDs with their randombly assigned class are returned.
    
    :param transcripts: List of IDs transcript.
    :type transcripts: list
    :param nb_to_choose: Number of transcripts to pick.
    :type nb_to_choose: int
    :return choosen: Picked transcripts with the corresponding randomly choosen class 
    :rtype: dict
    
    :Example:
    >>> transcripts = ['3R5.1a.1', '3R5.1b.1', 'B0348.1.1', 'B0348.10.1', 'B0348.2a.1', 'B0348.2b.1', 'B0348.4a.1'] # doctest: +SKIP 
    >>> choose_transcripts(transcripts,3) # doctest: +SKIP
    {'B0348.1.1': 0, 'B0348.10.1': 2, '3R5.1b.1': 0}
    >>> choose_transcripts(transcripts,100) # doctest: +SKIP
    {'3R5.1a.1': -1, '3R5.1b.1': 1, 'B0348.1.1': 0, 'B0348.10.1': 1, 'B0348.2a.1': 0, 'B0348.2b.1': 0, 'B0348.4a.1': 0}
    
    .. seealso: make_density_law(), annoToData()
    """
    law = make_density_law()

    if nb_to_choose == 0 or nb_to_choose >= len(transcripts):
        choosen = {t: int(rd.choice(law, 1)) for t in list(transcripts)}
    else:
        choosen = {
            t: int(rd.choice(law,1)) 
            for t in rd.choice(
                list(transcripts),
                int(nb_to_choose),
                replace=False
                )
            }

    return choosen


def parse_gtf_content(gtf_content, choosen, mix):
    """
    This function parse the ouput of the read_gtf() function (a list of each line of the gtf file).
    According to a dict with choosen ones transcripts and their associated class (the ouput of the 
    choose_transcripts() function), it gathers the exons of these transcripts, treats them according 
    to their class (construct_new_transcript() function) and constructs three lists of feature which
    will be the future thre output files : the reference file wich contains the pseudo-assembly (all the 
    longer transcripts, with or without retained introns), the library file from which the reads library
    will be generated (all the shorter transcripts, with or without spliced exons) and the control file 
    which contains all the retained introns and spliced exons.
    A bool parameter, mix, rules if the produced library list is mixed or not, that is, the library
    contains the non-0-class transcripts in two states : normal transcripts (all the exons and no intron)
    and modified transcript (with retained introns or with a spliced exon).
    
    :param gtf_content: All the line of a gtf file. Exactly the ouput of read_gtf().
    :type gtf_content: list
    :param choosen: All the choosen ones transcripts associated to their own class. Output of choose_transcripts()
    :type choosen: dict
    :param mix: Rules the production of a mixed library of reads.
    :type bool:
    :returns: reference, library, control
    :rtype: list, list, list
    
    :Example:
    
    >>> gtf,transcripts = read_gtf("../data/test_GBS.gtf") # doctest: +SKIP
    >>> choosen = choose_transcripts(transcripts,10) # doctest: +SKIP
    >>> reference, library, control = parse_gtf_content(gtf,choosen,False) # doctest: +SKIP
    
    .. seealso: read_gtf(), choose_transcripts(), construct_new_transcript(), annoToData()
    """
    reference = []
    library = []
    control = []
    exons = []
    # We add a false transcript feature at the end in the case of the true
    # last transcript of the file is choosen (in order to oblige the function to
    # processes last exons)
    gtf_content.append(["ref",".","transcript",".",".",".",".",".","."])
    for feature in gtf_content:
        if feature[2] == "exon" and feature[-1]["transcript_id"] in choosen:
            exons.append(feature)
        if feature[2] == "transcript" and exons:
            reference_transcript,library_transcript,feat_interest = construct_new_transcript(list(exons),choosen[exons[0][-1]["transcript_id"]])
            # We add corresponding transcript for each result file
            reference.extend(list(reference_transcript))
            library.extend(list(library_transcript))
            if feat_interest :
                control.extend(list(feat_interest))
            # If --mix-library flag is precised (i.e. the library has to contain the modified transcript and original transcript),
            # we add the modified reference_transcript to library only if it is really modified (if the class is not 0)
<<<<<<< HEAD
            if mix and reference_transcript[0][-1]["classe"] !='0':
=======
            if mix and reference_transcript[0][-1]["classe"] != '0':
>>>>>>> 835952fcc3cc9da4e376db380b58b2203add073c
                library.extend([e[:-1] + [{**e[-1],"transcript_id":e[-1]["transcript_id"]+".ref"}] for e in list(reference_transcript)])
            del choosen[exons[0][-1]["transcript_id"]]
            exons = []
        if not choosen:
            break

    return reference, library, control

def transcript_df(exons) :
    """
    From a transcrispt's exons list, returns a pandas.DataFrame corresponding to this transcript.
    This DataFrame contains all the exons and introns of the transcript(lines) and their GTF attributes (columns).
    It contains also an additional column "in_transcript" corresponding to a bool if a feature is an exon.
    
    :param exons: A list of exons with all Ensembl GTF format's fields which represents a transcript
    :type exons: list
    :return: whole_transcript
    :rtype: pandas.DataFrame
    
    :Example:
    
    >>> transcript = [ # doctest: +SKIP
    ... ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1'}], # doctest: +SKIP
    ... ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2'}]] # doctest: +SKIP 
    >>> transcript_df(transcript) # doctest: +SKIP
    ref        DB feature     start       end score strand frame                                          misc_attr  in_transcript
    0  III  WormBase    exon  13782587  13782934     .      +     .  {'gene_id': 'WBGene00007066', 'transcript_id':...           True
    1  III  WormBase  intron  13782935  13783360     .      +     .  {'gene_id': 'WBGene00007066', 'transcript_id':...          False
    2  III  WormBase    exon  13783361  13783459     .      +     .  {'gene_id': 'WBGene00007066', 'transcript_id':...           True
    
    
    .. seealso:: construct_new_transcript()
    """
    
    introns = []
    for i in range(len(exons) -1) :
        if exons[i][6] == "+" :
            intr_start = int(exons[i][4])+1 ; intr_end = int(exons[i+1][3])-1
        else :
            intr_start = int(exons[i+1][4])+1 ; intr_end = int(exons[i][3])-1
        introns.append([
            exons[i][0],
            exons[i][1],
            "intron",
            str(intr_start),
            str(intr_end),
            exons[i][5],
            exons[i][6],
            exons[i][7],
            {**{key: exons[i][-1][key] for key in exons[i][-1] if not key.startswith("exon")},"intron_id" : "+".join([exons[i][-1]["exon_id"],exons[i+1][-1]["exon_id"]])}
            ])
    whole_transcript = pd.DataFrame(
        [*[l for cpl in zip(exons,introns) for l in cpl],exons[-1]],
        columns=["ref","DB","feature","start","end","score","strand","frame","misc_attr"]
        )
    whole_transcript["in_transcript"] = whole_transcript.apply(lambda df : df.feature == "exon",axis=1)
    whole_transcript["misc_attr"] = whole_transcript.apply(lambda df : dict(df.misc_attr),axis=1) # Just to fix the mutability of misc_attr dict
    return whole_transcript

def construct_new_transcript(exons, classe):
    """
    From a transcript (a list of exons) and its class, modify the transcript according to its class
    by splicing one exon (class -1) or retaining some introns (class greater than 0). It returns two 
    state for a transcript : one for the reference pseudo-assembly and the other one for the library
    file from which the reads will be generated. It returns also the modified feature(s) (spliced exon or retained introns) 
    for the control file.
    If the class is 0, the function returns in lib_t and ref_t the unmodifed transcript.
    In the same vein, if a modification cannot be performed because of transcript construction, the minimum suitable modification 
    is performed and the class is changed. For example if, a transcript contains only one exon, whatever its class,
    its returned without modification and its class is reset to 0. Another example is when the class is greater than the number of intron :
    if a transcript has 2 exons (i.e. only 1 intron) but its class is 3, the intron is retained and its class becomes 1.
    
    :param exons: A collection of exons represents a transcript
    :type exons: list
    :param classe: A number indicates the modification to perform on transcript
    :type classe: int
    :return ref_t: Transcript for reference file
    :rtype ref_t: list
    :return lib_t: Transcript for library file
    :rtype lib_t: list
    :return ft_on_t: Feature of interest (spliced exon or retained intron) for control file
    :rtype ft_on_t: list
    
    :Example:
    >>> transcript = [ # doctest: +SKIP
    ... ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1'}], # doctest: +SKIP
    ... ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2'}]] # doctest: +SKIP 
    >>> reference,library,feature = construct_new_transcript(list(transcript),-1) # doctest: +SKIP 
    >>> reference # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'classe': '-1'}], 
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'classe': '-1'}]]
    >>> library # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'classe': '-1'}]]
    >>> feature # doctest: +SKIP 
    [['3R5.2', 'WormBase', 'spliced_exon', '1', '348', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'classe': '-1'}]]
    >>> reference,library,feature = construct_new_transcript(list(transcript),0) # doctest: +SKIP 
    >>> reference # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'classe': '0'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'classe': '0'}]]
    >>> library # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'classe': '0'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'classe': '0'}]]
    >>> feature # doctest: +SKIP 
    >>> reference,library,feature = construct_new_transcript(list(transcript),3) # doctest: +SKIP 
    >>> reference # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'classe': '1'}],
     ['III', 'WormBase', 'exon', '13782935', '13783360', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'intron_id': '3R5.2.e1+3R5.2.e2', 'classe': '1'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'classe': '1'}]]
    >>> library # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'classe': '1'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'classe': '1'}]]
    >>> feature # doctest: +SKIP 
    [['3R5.2', 'WormBase', 'retained_intron', '349', '774', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'intron_id': '3R5.2.e1+3R5.2.e2', 'classe': '1'}]]
    
    .. seealso:: transcript_df(), parse_gtf_content()
    .. note:: When a intron is retained, is annotated like exon in the feature field to oblige gffread to consider it like an exon and add it in the FASTA transcript. 
    """
    lib_t = ref_t = ft_on_t = [] # a virer une fois fini
    whole_transcript = transcript_df(list(exons))
    if classe == 0 or len(exons) == 1:
        whole_transcript.apply(lambda df : df.misc_attr.update({"classe":str(0)}),axis=1) # We update the class of the transcript for each feature
        
        # From the DataFrame, we take only the feature with True value in "in_transcript" column (here it's only exons because of the class 0)  
        lib_t = ref_t = [ 
            list(line)[1:] # We exclude the index value of DataFrame which is stored in the first position of the record
            for line in whole_transcript.loc[lambda df : df.in_transcript == True,"ref":"misc_attr"].to_records() 
            ]
        ft_on_t = None
    elif classe == -1:
        whole_transcript.apply(lambda df : df.misc_attr.update({"classe":str(classe)}),axis=1) # We update the class of the transcript for each feature
        
        # We add to reference file the not spliced transcript
        ref_t = [ 
            list(line)[1:] # We exclude the index value of DataFrame which is stored in the first position of the record
            for line in whole_transcript.loc[lambda df : df.in_transcript == True,"ref":"misc_attr"].to_records() 
            ]
        
        # We randomly pick an index of the DataFrame wich corresponds to an exon for splicing 
        spliced_exon = int(rd.choice(whole_transcript.loc[lambda df : df.feature == "exon"].index,1))
        # Splice the choosen exon by puting False in "in_transcript" column
        whole_transcript.at[spliced_exon,"in_transcript"] = False
        
        # Construction of the feature which corresponds to spliced exon for the control GFF file
        e_start = sum(whole_transcript.loc[lambda df : (df.in_transcript == True) & (df.index < spliced_exon),"end"].apply(int)-whole_transcript.loc[lambda df : (df.in_transcript == True) & (df.index < spliced_exon),"start"].apply(int)+1)+1
        e_end = e_start + (int(whole_transcript.at[spliced_exon,"end"])-int(whole_transcript.at[spliced_exon,"start"]))
        ft_on_t = [[
            whole_transcript.at[spliced_exon,"misc_attr"]["transcript_id"],
            whole_transcript.at[spliced_exon,"DB"],
            "spliced_exon",
            str(e_start),
            str(e_end),
            whole_transcript.at[spliced_exon,"score"],
            whole_transcript.at[spliced_exon,"strand"],
            whole_transcript.at[spliced_exon,"frame"],
            whole_transcript.at[spliced_exon,"misc_attr"]
            ]]
        
        # We add to library file (from which the read will be generated) the spliced transcript
        lib_t = [ 
            list(line)[1:] # We exclude the index value of DataFrame which is stored in the first position of the record
            for line in whole_transcript.loc[lambda df : df.in_transcript == True,"ref":"misc_attr"].to_records() 
            ]
    elif classe > 0:
        # if the number of retained introns is greater than the number of introns (the class), we take all the available intron and change the class if it is necessary
        if classe >= len(exons) - 1:
            classe = len(exons) - 1
            choosen_introns = list(whole_transcript.loc[lambda df : df.feature == "intron"].index)
        # if not, we random pick the desired number of introns
        else:
            choosen_introns = sorted(list(rd.choice(whole_transcript.loc[lambda df : df.feature == "intron"].index,classe,replace=False)))
        
        whole_transcript.apply(lambda df : df.misc_attr.update({"classe":str(classe)}),axis=1) # We update the class of the transcript for each feature
        
        # We add to library file (from which the read will be generated) the transcript without retained introns
        lib_t = [ 
            list(line)[1:] # We exclude the index value of DataFrame which is stored in the first position of the record
            for line in whole_transcript.loc[lambda df : df.in_transcript == True,"ref":"misc_attr"].to_records() 
            ]
        
        # Intron retention by putting true in "in_transcript" column and 
        # change feature field from exon to intron in order to oblige gffread to consider the retained intron
        ft_on_t=[]
        for retained_intron in choosen_introns :
            whole_transcript.at[retained_intron,"in_transcript"] = True
            whole_transcript.at[retained_intron,"feature"] = "exon"
            
            # Construction of the feature which corresponds to spliced exon for the control GFF file
            i_start = sum(whole_transcript.loc[lambda df : (df.in_transcript == True) & (df.index < retained_intron),"end"].apply(int)-whole_transcript.loc[lambda df : (df.in_transcript == True) & (df.index < retained_intron),"start"].apply(int)+1)+1
            i_end = i_start + (int(whole_transcript.at[retained_intron,"end"])-int(whole_transcript.at[retained_intron,"start"]))
            ft_on_t.append([
                whole_transcript.at[retained_intron,"misc_attr"]["transcript_id"],
                whole_transcript.at[retained_intron,"DB"],
                "retained_intron",
                str(i_start),
                str(i_end),
                whole_transcript.at[retained_intron,"score"],
                whole_transcript.at[retained_intron,"strand"],
                whole_transcript.at[retained_intron,"frame"],
                whole_transcript.at[retained_intron,"misc_attr"]
                ])
        
        # We add to reference file the transcript with retained_introns
        ref_t = [ 
            list(line)[1:] # We exclude the index value of DataFrame which is stored in the first position of the record
            for line in whole_transcript.loc[lambda df : df.in_transcript == True,"ref":"misc_attr"].to_records() 
            ]
    
    return ref_t, lib_t, ft_on_t


def write_gtf_file(content: str, output: str):
    gtf = "\n".join(["\t".join([str(item) if not isinstance(item, dict) else "; ".join(
        ['{k} "{v}"'.format(k=key, v=item[key]) for key in item]) for item in ligne])for ligne in content])
    with open(output, "w") as gtf_file:
        gtf_file.write(gtf + "\n")


def extractFasta(genome: str, ref_file, lib_file, path: str):
    sp.call(["gffread", ref_file, "-g", genome,
             "-w", path + "_reference.fa", "-F"])
    sp.call(["gffread", lib_file, "-g", genome,
             "-w", path + "_library.fa", "-F"])


def annoToData(annotation: str, fasta: str, nb: int, output: str, mix: bool):
    print("GTF reading...")
    gtf_content, transcripts = read_gtf(annotation)
    print("Generate transcripts...")
    choosen = choose_transcripts(transcripts, nb)
    reference, library, control = parse_gtf_content(gtf_content, choosen, mix)

    lib_tmpfile = tmp.NamedTemporaryFile(dir=".",delete=True) ; ref_tmpfile = tmp.NamedTemporaryFile(dir=".",delete=True)
    write_gtf_file(reference,ref_tmpfile.name) ;  write_gtf_file(library,lib_tmpfile.name)
    write_gtf_file(control,output+"_Features-of-interest.gtf")
    print("FASTA files writing with gffread...")
    extractFasta(fasta,ref_tmpfile.name,lib_tmpfile.name,output)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Functions to find and analyse split in alignment file.")
    subparser = parser.add_subparsers(help="sub-command help")

    parser_annotodata = subparser.add_parser(
        "annoToData",
        help="Data simulation (genome assembly and reads library) from a genome annotation (GFF or GTF file)")
    parser_annotodata.add_argument(
        "-i",
        "--annotation",
        help="Filename of GFF file which contains the genome annotation",
        type=str,
        metavar='',
        required=True)
    parser_annotodata.add_argument(
        "-f",
        "--fasta",
        type=str,
        metavar='',
        required=True,
        dest="fasta")
    group_nb = parser_annotodata.add_mutually_exclusive_group(required=True)
    group_nb.add_argument(
        "-n",
        "--nb_genes",
        help="Total number of genes to transcript ",
        type=float,
        required=False,
        default=0,
        dest="nb",
        metavar='')
    group_nb.add_argument(
        "-a",
        "--all",
        help="Flag which says if all genes from GFF have to be transcripted",
        action="store_const",
        const=0,
        default=False,
        dest="nb")
    parser_annotodata.add_argument(
        "-o",
        "--output",
        help="prefix of ouput files",
        type=str,
        required=True,
        dest="output",
        metavar='')
    # ~ parser_annotodata.add_argument("--no-grinder",help="Boolean which rules if let the library at temporary state (i.e. before the reads library generation)",
    # ~ action="store_false",default = True, dest="grinder") ;
    parser_annotodata.add_argument("--mix-library", help="Boolean which rules if the generated library is mixed i.e. if the library contains the transcript in two state when a intron is retained or an exon is spliced",
                                   action="store_true", default=False, dest="mix")
    # ~ parser_annotodata.set_defaults(func=ds.annoToData)

    args = parser.parse_args()

    annoToData(args.annotation, args.fasta, args.nb, args.output, args.mix)
