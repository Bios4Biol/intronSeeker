#!/usr/bin/env python3

# Modules
try :
    import gzip ;
    import os ;
    import sys ;
    import subprocess ;
    import random ;
    import numpy as np ;
    import numpy.random as rd ;
    import configparser ;
    import re ;
    import pandas as pd ;
    import tempfile as tmp
    import subprocess as sp
    from Bio.Seq import Seq ;
    from Bio import SeqIO ;
    from Bio.SeqRecord import SeqRecord ;
    from collections import defaultdict ;
except ImportError as error :
    print(error) ;
    exit(1) ;

from pprint import pprint


##############################
### Full Random Simulation ###
##############################

def random_seq(length: int, letter: str) :
    """
    Create a sequence based on letter and as long as length

    :param length: length of the sequence
    :type length: int
    :param letter: letter used in the sequence
    :type letter: str
    :return: sequence
    :rtype: str
    
    :Example: 
    >>> random_seq(10,"A")
    'AAAAAAAAAA'
    >>> s = random_seq(100,"ATCG")
    >>> len(s)
    100
    >>> s.strip("ATCG")
    ''
    
    .. seealso:: insert_intron(), full_random_simulation()
    """
    sequence = "" ;
    for l in range(0, length) :
        sequence += random.choice(letter) ; # Random choice of each letter 
    return sequence ;

def insert_intron(contig_seq : str, lower : int, upper: int,):
    """
    Insert in a sequence a random intron with a length beetween [lower,upper]. 
    Return the new sequence with the inserted intron and the coordinates of the intron
    in the new sequence. The coordinates are in the BED format : begin at 0 position and 
    the end coordinate is not contained in the intron.

    :param contig_seq: Sequence in which a intron will be inserted.
    :type contig_seq: str
    :param lower: Minimal length of intron
    :type lower: int
    :param upper: Maximal length of intron
    :type upper: int
    :return: 
        * [0] new_seq 
        * [1] insert_pos
        * [2] intron_end
    :rtype : (str,int,int) 
    
    :Example:
    >>> insert_intron("ATCCCGGTGGTGAATGTCGTGATGC",5,8) # doctest: +SKIP
    ('ATCGTCGGCAGCCGGTGGTGAATGTCGTGATGC', 3, 11)
    
    .. seealso:: random_seq(), full_random_simulation()
    """
    rand_length = random.randint(lower, upper) - 4 # Random choice of the intorn length
    intron_seq = "".join(["GT", random_seq(rand_length, "ATCG"),"AG"]) # Generate the random sequence
    insert_pos = random.randint(1, len(contig_seq)) # Random choice of the insertion position of the intron
    new_seq = "".join([contig_seq[:insert_pos], intron_seq, contig_seq[insert_pos:]]) # Insertion of the intron
    intron_end = insert_pos+rand_length+4
    return new_seq, insert_pos, intron_end

def full_random_simulation(nb:int, maxi:int, mini:int, part:int, lower:int, upper:int, prefix:str, output:str, force:bool, mix: bool) :
    """
    Simulate a set of nb contigs (with entirely random sequence) with random length beetween [mini,maxi].
    In each contig, an intron with a random length beetween [lower,upper] is randomly inserted. 
    If part isn't 100%, only a random part of simulated contigs have a intron insertion.  
    Write the results in two files : 
        - output_contigs.fa : FASTA of the simulated contigs
        - output_introns.txt : Tabulate-separated values which define introns (contig,start,end,reverse) 
    
    :param nb: Number of contigs to simulate
    :type nb: int
    :param maxi: Maximum length of contigs
    :type maxi: int
    :param mini: Minimum length of contigs
    :type mini: int
    :param print: int the intron insertion is performed only on radom part% of the contigs
    :type part: int
    :param lower: Minimum length of inserted introns
    :type lower: int
    :param upper: Maximum length of inserted introns
    :type upper: int
    :param output: Basename of the two output files
    :type output: str
    :param mix: Boolean which rules if a mixed library is generated.
    :type mix: bool
    
    :Example:
    >>> full_random_simulation(10,150,1000,False,150,1000,"FullRandomSimulation") # doctest: +SKIP
    
    .. seealso:: random_seq(), insert_intron(), annoToData()
    """
    output_path = output + "/frs";
    if prefix:
        output_path += "_" + prefix;
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.makedirs(output)
    if not force:
        try :
            if os.path.exists(output_path + "_contigs.fa") or os.path.exists(output_path + "_contigs-modified.fa") or os.path.exists(output_path + "_modifications.txt") :
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
        
    # Check if the length intervals are correct
    if mini > maxi :
        print('*WARNING* : Contigs maximal length is inferior to minimal length. Max and Min values are swapped.')
        M = maximum
        maxi = mini
        mini = M
    
    if lower > upper :
        print('*WARNING* : Introns maximal length is inferior to minimal length. lower and upper values are swapped.')
        U = upper
        upper = lower
        lower= U
    
    # Define the distribution of the intron insertion among the contigs
    # Two cases : all the contigs are modified or only a random part 
    distrib = [1]*int(nb*part/100)
    distrib.extend([0]*(nb-len(distrib)))
    random.shuffle(distrib)
    
    #Generate the contigs
    reference_contigs_set = []
    library_contigs_set = []
    introns = []
    for c in range(0,nb) :
        length = random.randint(mini, maxi) ; # Random contig length 
        name = "SEQUENCE" + str(c+1) ; # Contig name creation
        contig_seq = random_seq(length, "ATCG") # Random sequence creation
        library_seq = Seq(contig_seq)
        
        # Insert intron according to the distribution
        if distrib[c] : 
            modified_seq, intron_start, intron_end = insert_intron(contig_seq, lower, upper)
            reference_seq = Seq(modified_seq)
        else :
            intron_start,intron_end = None,None
            reference_seq = Seq(contig_seq)
        
        # Half of the contigs are reversed (i.e. like if they comes from - strand)
        reverse = random.choice([True,False])
        if reverse : 
            reference_seq = Seq.reverse_complement(reference_seq)
            library_seq = Seq.reverse_complement(library_seq)
            if intron_start and intron_end :
                old_start = intron_start
                intron_start = len(reference_seq)-intron_end
                intron_end = len(reference_seq)-old_start

        library_contigs_set.append(SeqRecord(library_seq,id=name.split()[0],description="reverse="+str(reverse)))
        
        description = " ".join(["intron_start="+str(intron_start),"intron_end="+str(intron_end),"reverse="+str(reverse)])
        if intron_start == None :
            reference_contigs_set.append(SeqRecord(reference_seq,id=name.split()[0],description=description))
        else:
            reference_contigs_set.append(SeqRecord(reference_seq,id=name.split()[0]+".modif",description=description))
            if mix :
                reference_contigs_set.append(SeqRecord(reference_seq,id=name.split()[0],description=description))
     
        if distrib[c] :
            r = "+"
            if(reverse) :
                r = "-"
            introns.append("\t".join([name+".modif","frs","retained_intron",str(intron_start),str(intron_end),".",r,".","."]))
        
    SeqIO.write(reference_contigs_set,output_path+"_contigs-modified.fa","fasta")
    SeqIO.write(library_contigs_set,output_path+"_contigs.fa","fasta")
    
    with open(output_path+"_contigs-modified.gtf","w") as out :
        out.write("\n".join(introns))
    


################################
### Reads library Simulation ###
################################

#~ def _grinder(input_file : str, profile_file : str, output_file : str) :

    #~ os.system("grinder -rf {input_file} -pf {profile_file} -bn {output_file}".format(
        #~ input_file=input_file, profile_file=profile_file, output_file=output_file
        #~ )) ;

def grinder(rf: str, pf: str, prefix: str, output: str, force: bool):
    """
    Generate reads for a reference file depending on grinder parameters from the profile file.
    Call the software Grinder.
    :param rf: reference file
    :param pf:profile file containing all grinder parameters
    :param pref: prefix of the output files
    :return:
    """    
    output_path = output + "/sr"
    if prefix:
        output_path += "_" + prefix

    # Create output dir if not exist
    if not os.path.exists(output) :
        os.makedirs(output)
    if not force:
        try :
            if os.path.exists(output_path + "_ranks.txt") or os.path.exists(output_path + "_R1.fastq.gz") or os.path.exists(output_path + "_R2.fastq.gz") :
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)

    os.system("grinder -rf {input_file} -pf {profile_file} -bn {output_file} > {log}".format(
        input_file=rf.name, profile_file=pf.name, output_file=output_path, log=output_path + ".log"
        ))
    os.rename(output_path + "-ranks.txt", output_path + "_ranks.txt")
    split_read(output_path)


def split_read(output_path : str):
    """
    Split a paired-end fasta file in two fasta files.q
    :param input_file: name/path of the fasta file
    :param pref: prefix oh the two output file
    :return:
    """
    outdir = os.path.dirname(output_path) ;
    pref = os.path.basename(output_path) ;
    input_file = outdir + "/" + [n for n in os.listdir(outdir) if n.startswith(pref+"-reads.")].pop() ;
    pile = 0 ;
    if input_file.endswith(".fa") :
        with open(input_file, 'rU') as filer, gzip.open(output_path + "_R1.fa.gz", 'wt') as read1, \
                gzip.open(output_path + "_R2.fa.gz", 'wt') as read2:
            for record in SeqIO.parse(filer, "fasta"):
                if pile % 2 == 0:
                    SeqIO.write(record, read1, "fasta") ;
                else:
                    SeqIO.write(record, read2, "fasta") ;
                pile += 1
    elif input_file.endswith(".fastq") :
        with open(input_file, 'rU') as filer, gzip.open(output_path + "_R1.fastq.gz", 'wt') as read1, \
                gzip.open(output_path + "_R2.fastq.gz", 'wt') as read2:
            for record in SeqIO.parse(filer, "fastq"):
                if pile % 2 == 0:
                    SeqIO.write(record, read1, "fastq") ;
                else:
                    SeqIO.write(record, read2, "fastq") ;
                pile += 1 ; 
    os.remove(input_file) ;


###############################
### Genome-Based Simulation ###
###############################

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
    config.read(os.path.abspath(os.path.dirname(__file__))+"/../config/intronSeeker.properties")
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
            library.extend(list(library_transcript))
            if reference_transcript[0][-1]["class"] != '0':
                reference.extend([e[:-1] + [{**e[-1],"transcript_id":e[-1]["transcript_id"]+".modif"}] for e in list(reference_transcript)])
            else:
                reference.extend(list(reference_transcript))
            if feat_interest :
                control.extend([e[:-1] + [{**e[-1],"transcript_id":e[-1]["transcript_id"]+".modif"}] for e in list(feat_interest)])
            # If --mix-library flag is precised (i.e. the library has to contain the modified transcript and original transcript),
            # we add the modified reference_transcript to library only if it is really modified (if the class is not 0)
            if mix and reference_transcript[0][-1]["class"] != '0':
                library.extend([e[:-1] + [{**e[-1],"transcript_id":e[-1]["transcript_id"]+".modif"}] for e in list(reference_transcript)])
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
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'class': '-1'}], 
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'class': '-1'}]]
    >>> library # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'class': '-1'}]]
    >>> feature # doctest: +SKIP 
    [['3R5.2', 'WormBase', 'spliced_exon', '1', '348', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'class': '-1'}]]
    >>> reference,library,feature = construct_new_transcript(list(transcript),0) # doctest: +SKIP 
    >>> reference # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'class': '0'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'class': '0'}]]
    >>> library # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'class': '0'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'class': '0'}]]
    >>> feature # doctest: +SKIP 
    >>> reference,library,feature = construct_new_transcript(list(transcript),3) # doctest: +SKIP 
    >>> reference # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'class': '1'}],
     ['III', 'WormBase', 'exon', '13782935', '13783360', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'intron_id': '3R5.2.e1+3R5.2.e2', 'class': '1'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'class': '1'}]]
    >>> library # doctest: +SKIP 
    [['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'class': '1'}],
     ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'class': '1'}]]
    >>> feature # doctest: +SKIP 
    [['3R5.2', 'WormBase', 'retained_intron', '349', '774', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'intron_id': '3R5.2.e1+3R5.2.e2', 'class': '1'}]]
    
    .. seealso:: transcript_df(), parse_gtf_content()
    .. note:: When a intron is retained, is annotated like exon in the feature field to oblige gffread to consider it like an exon and add it in the FASTA transcript. 
    """
    lib_t = ref_t = ft_on_t = [] # a virer une fois fini
    whole_transcript = transcript_df(list(exons))
    if classe == 0 or len(exons) == 1:
        whole_transcript.apply(lambda df : df.misc_attr.update({"class":str(0)}),axis=1) # We update the class of the transcript for each feature
        
        # From the DataFrame, we take only the feature with True value in "in_transcript" column (here it's only exons because of the class 0)  
        lib_t = ref_t = [ 
            list(line)[1:] # We exclude the index value of DataFrame which is stored in the first position of the record
            for line in whole_transcript.loc[lambda df : df.in_transcript == True,"ref":"misc_attr"].to_records() 
            ]
        ft_on_t = None
    elif classe == -1:
        whole_transcript.apply(lambda df : df.misc_attr.update({"class":str(classe)}),axis=1) # We update the class of the transcript for each feature
        
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
            whole_transcript.at[spliced_exon,"misc_attr"]["transcript_id"]+".modif",
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
        
        whole_transcript.apply(lambda df : df.misc_attr.update({"class":str(classe)}),axis=1) # We update the class of the transcript for each feature
        
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
                whole_transcript.at[retained_intron,"misc_attr"]["transcript_id"]+".modif",
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
    """
    Write a GTF file from a list of features.
    
    :param content: List of features. Each feature is a list where each element corresponds to each field of GTF (the last element is a dict for the descrption GTF field)
    :type content: list
    :param output: Name of the output file
    :type output: str
    
    :Example:
    >>> content=[["chr1","RefSeq","exon","2","200",".","+",".",{"gene":"g1","exon":"g1.e1"}]] # doctest: +SKIP 
    >>> write_gtf_file(content,"exon.gtf") # doctest: +SKIP 
    -------------   exon.gtf   ---------------------------
    chr1	RefSeq	exon	2	200	.	+	.	gene "g1"; exon "g1.e1"
    ------------------------------------------------------
    
    .. seealso:: gtf_based_simulation(), parse_gtf_content()
    """
    gtf = "\n".join(["\t".join([str(item) if not isinstance(item, dict) else "; ".join(
        ['{k} "{v}"'.format(k=key, v=item[key]) for key in item]) for item in ligne])for ligne in content])
    with open(output, "w") as gtf_file:
        gtf_file.write(gtf + "\n")


def extract_fasta(genome: str, mix: bool, ref_file: str, lib_file: str, path: str):
    """
    Call gffread program to generate a FASTA file from a genome (FASTA format) and a GTF file.
    The output files contains the sequence of each transcript defined by the entry GTF file.
    
    :param genome: Filename of the genome FASTA file
    :type genome: str
    :param ref_file: GTF filename of the reference transcripts. Corresponds to the produced reference.fasta file.
    :type ref_file: str
    :param lib_file: GTF filename of the library transcripts. Corresponds to the produced library.fasta file.
    :type lib_file: str
    :param path: Path of the directory where the two FASTA files will be stored.
    :type path: str
    """
    sp.call(["gffread", ref_file, "-g", genome,
             "-w", path + "_transcripts-modified.fa", "-F"])
    if mix:
        sp.call(["gffread", lib_file, "-g", genome,
             "-w", path + "_transcripts-mix-state.fa", "-F"])
    else:
        sp.call(["gffread", lib_file, "-g", genome,
             "-w", path + "_transcripts.fa", "-F"])


def gtf_based_simulation(annotation: str, fasta: str, nb: int, prefix: str, output: str, force: bool, mix: bool):
    """
    Simulate a RNA-seq pseudo-assembly from a GTF file with retained introns or spliced exons. This procedure produces 3 files :
        - output_reference.fasta : pseudo-assembly where the contigs have potentially retained introns. 
        - output_library.fasta : pseudo-assembly where the contigs have potentially spliced exons. This file have to be used to generate reads library with simulateReads module
        - output_Features_of_interest.gtf : gtf which conatains all simulated retained introns and spliced exons.
    To choose the number of pseudo-contigs, you have to use the nb argument. If nb=0, all the transcripts contained in the GTF will be pseudo-transcripted.
    The mix arguments rules if a mixed library is genrated : by default, the library file contains the shorter contigs (without intron or with spliced exon) compared to reference file ; 
    a mixed library file contains both of the peudo-contig states  :the shorter one (without intron and spliced exon) and the longer on (like in reference file : with retained intron or with all exons).
    The goal is to produce some reads with the intron or the exon (which will not split in the reference/library alignment) and reads without the intron or the exon (which will split during the alignment).
    
    :param annotation: GTF filename which contains genome annotation
    :type annotation: str
    :param fasta: FASTA file which contains genome sequence
    :type fasta: str
    :param nb: Number of pseudo-contigs to generate. If 0, all the transcripts in GTF will give pseudo contig.
    :type nb: str
    :param output: Basename of the output files
    :type output: str
    :param mix: Boolean which rules if a mixed library is generated.
    :type mix: bool
    
    """
    output_path = output + "/gbs";
    if prefix:
        output_path += "_" + prefix;
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.makedirs(output)
    if not force:
        try :
            if os.path.exists(output_path + "_transcripts-modified.fa") or os.path.exists(output_path + "_transcripts.fa") or os.path.exists(output_path + "_transcripts-modified.gtf") or  os.path.exists(output_path + "_transcripts-mix-state.fa"):
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    #print("GTF reading...")
    gtf_content, transcripts = read_gtf(annotation.name)
    #print("Generate transcripts...")
    choosen = choose_transcripts(transcripts, nb)
    reference, library, control = parse_gtf_content(gtf_content, choosen, mix)

    lib_tmpfile = tmp.NamedTemporaryFile(dir=output,delete=True)
    ref_tmpfile = tmp.NamedTemporaryFile(dir=output,delete=True)
    write_gtf_file(reference,ref_tmpfile.name)
    write_gtf_file(library,lib_tmpfile.name)
    write_gtf_file(control,output_path + "_transcripts-modified.gtf")
    #print("FASTA files writing with gffread...")
    extract_fasta(fasta.name, mix, ref_tmpfile.name, lib_tmpfile.name, output_path)
    
