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

def full_random_simulation(nb : int, maxi : int, mini : int, half : bool, lower : int, upper : int, output : str) :
    """
    Simulate a set of nb contigs (with entirely random sequence) with random length beetween [mini,maxi].
    In each contig, an intron with a random length beetween [lower,upper] is randomly inserted. 
    If half is True, only a random half of simulated contigs have a intron insertion.  
    Write the results in two files : 
        - output_contigs.fa : FASTA of the simulated contigs
        - output_introns.txt : Tabulate-separated values which define introns (contig,start,end,reverse) 
    
    :param nb: Number of contigs to simulate
    :type nb: int
    :param maxi: Maximum length of contigs
    :type maxi: int
    :param mini: Minimum length of contigs
    :type mini: int
    :param half: Boolean which rules if the intron insertion is performed on a random half of the contigs
    :type half: bool
    :param lower: Minimum length of inserted introns
    :type lower: int
    :param upper: Maximum length of inserted introns
    :type upper: int
    :param output: Basename of the two output files
    :type output: str
    
    :Example:
    >>> full_random_simulation(10,150,1000,False,150,1000,"FullRandomSimulation") # doctest: +SKIP
    
    .. seealso:: random_seq(), insert_intron(), annoToData()
    """
    outdir = output.split(".")[0] + "_FullRandomSimulation" ; # Output directory name (where this function will write 
                                                # all the results and tmp files).
    if not os.path.exists(outdir) :
        os.mkdir(outdir) 
    output_path = "/".join([outdir,output]) ; # Path to output file
    
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
    # Two cases : all the contigs are modified or only a random half 
    if half :
        distrib = [0]*int(nb/2)
        distrib.extend([1]*(nb-len(distrib)))
        random.shuffle(distrib)
    else :
        distrib = [1]*nb
    
    #Generate the contigs
    contigs_set = []
    introns = ["Contig\tStart\tEnd\tReverse"]
    for c in range(0,nb) :
        length = random.randint(mini, maxi) ; # Random contig length 
        name = "SEQUENCE" + str(c+1) ; # Contig name creation
        contig_seq = random_seq(length, "ATCG") # Random sequence creation
        
        # Insert intron according to the distribution
        if distrib[c] : 
            contig_seq, intron_start, intron_end = insert_intron(contig_seq, lower, upper)
            description = " ".join(["intron_start="+str(intron_start),"intron_end="+str(intron_end)])
        else :
            description = " ".join(["intron_start="+str(None),"intron_end="+str(None)])
        contig_seq = Seq(contig_seq)
        
        # Half of the contigs are reversed (i.e. like if they comes from - strand)
        reverse = random.choice([True,False])
        if reverse : 
            contig_seq = Seq.reverse_complement(contig_seq)
        description = " ".join([description,"reverse="+str(reverse)])
        
        contigs_set.append(SeqRecord(contig_seq,id=name.split()[0],description=description))
        if distrib[c] :
            introns.append("\t".join([name,str(intron_start),str(intron_end),str(reverse)]))
        
    SeqIO.write(contigs_set,output_path+"_contigs.fa","fasta")
    
    with open(output_path+"_introns.txt","w") as out :
        out.write("\n".join(introns))
    


################################
### Reads library Simulation ###
################################

def _grinder(input_file : str, profile_file : str, output_file : str) :

    os.system("grinder -rf {input_file} -pf {profile_file} -bn {output_file}".format(
        input_file=input_file, profile_file=profile_file, output_file=output_file
        )) ;

def grinder(rf: str, pf: str, pref: str):
    """
    Generate reads for a reference file depending on grinder parameters from the profile file.
    Call the software Grinder.
    :param rf: reference file
    :param pf:profile file containing all grinder parameters
    :param pref: prefix of the output files
    :return:
    """
    outdir = pref + "_grinder" ;
    os.system("mkdir " + outdir ) ;
    output_path = outdir + "/" + pref ;

    _grinder(rf, pf, output_path) ;
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
        with open(input_file, 'rU') as filer, gzip.open(output_path + "_read_1.fa.gz", 'wt') as read1, \
                gzip.open(output_path + "_read_2.fa.gz", 'wt') as read2:
            for record in SeqIO.parse(filer, "fasta"):
                if pile % 2 == 0:
                    SeqIO.write(record, read1, "fasta") ;
                else:
                    SeqIO.write(record, read2, "fasta") ;
                pile += 1
    elif input_file.endswith(".fastq") :
        with open(input_file, 'rU') as filer, gzip.open(output_path + "_read_1.fastq.gz", 'wt') as read1, \
                gzip.open(output_path + "_read_2.fastq.gz", 'wt') as read2:
            for record in SeqIO.parse(filer, "fastq"):
                if pile % 2 == 0:
                    SeqIO.write(record, read1, "fastq") ;
                else:
                    SeqIO.write(record, read2, "fastq") ;
                pile += 1 ; 
    os.remove(input_file) ;



# ~ def intron_unique(pref: str):
    # ~ """
    # ~ Create the gff file of the intron from star alignment
    # ~ :param pref: output prefix
    # ~ :return:
    # ~ """
    # ~ data = pd.read_table(pref + ".txt", names=['ref', 'start', 'end', 'query', 'qstart', 'qend'])
    # ~ data = data.drop_duplicates(subset=('ref', 'start', 'end', 'qstart', 'qend'))
    # ~ data = pd.DataFrame(
        # ~ {'col0_ref': data.ref, 'col1': ".", "col2_class": "Intron", "col3_start": data.qstart, "col4_end": data.qend,
         # ~ 'col5_quality': "5000", 'col6': ".", 'col7': ".", 'col8': "."})
    # ~ data.to_csv(pref+".gff", header=False, sep="\t", index=False)


###############################
### Genome-Based Simulation ###
###############################

def countGenes(gff : str) :
    """
    Count the protein-coding genes in the GFF file in order to randomly
     pick among them for simulate the transcripts.
    
    :param gff: Name of the genome reference GFF file.
    :return: Number of genes. 
    """
    
    grep = subprocess.Popen(["grep", "gene.*=protein_coding;",gff], stdout=subprocess.PIPE) ;
    wcline = subprocess.Popen(["wc", "-l"], stdin=grep.stdout, stdout=subprocess.PIPE) ;
    # Allow grep to receive a SIGPIPE if wcline exits before grep
    grep.stdout.close() ; 
    # Run the commands 
    sdo = wcline.communicate()[0] ;
    # Extract the result from standard output
    nb_genes = int(re.search(r'(\d+)',str(sdo)).group(1)) ;
    
    return nb_genes ;


def chooseGenes(gff : str, nb_genes = int) :
    """
    Randomly pick 
    """
    try :
        total_genes = countGenes(gff) ;
        if nb_genes < 0 :
            raise ValueError ;
        elif nb_genes >= total_genes or nb_genes == 0 :
            choosen = range(total_genes);
        else :
            choosen = rd.choice(range(total_genes),int(nb_genes),replace=False) ;
    except ValueError :
        print("*Value Error* : Number of genes to select can't be negative.") ;
        exit(1);

    return choosen ;


def generateTranscripts(gff_file : str, choosen_genes : list, output_path : str, mix : bool) :
    classes_transcripts = makeDensityLaw() ;
    reference = "" ;
    library = "" ;
    all_features = "" ;
    with open(gff_file,"r") as gff :
        ligne = gff.readline().rstrip() ;
        
        i = 0 ;
        num_transcript = 0 ;
        while ligne :
            if ligne.find("\tgene\t") != -1 and ligne.find("=protein_coding;") != -1 :
                if i in choosen_genes :
                    introns = None ; exon = None ;
                    gene = [] ;
                    nb_mRNA = 0 ;
                    ligne = gff.readline().rstrip() ;
                    while ligne and (ligne.find("\tgene\t") == -1 and ligne.find("\tpseudogene\t") == -1 ) :
                        if ligne.find("\ttranscript\t") != -1 or ligne.find("\ttRNA\t") != -1 or ligne.find("\trRNA\t") != -1:
                            while ligne and (ligne.find("\tgene\t") == -1 and ligne.find("\tmRNA\t") == -1 ) :
                                ligne = gff.readline().rstrip() ;
                        if ligne.find("\tmRNA\t") != -1  :
                            nb_mRNA += 1 ;
                        gene.append(ligne) ;
                        ligne = gff.readline().rstrip() ;
                    transcript = "" ;
                    num_transcript += 1 ;
                    mRNA,nb_exons = parseGene(gene,nb_mRNA,num_transcript) ;
                    classe = int(rd.choice(classes_transcripts, 1)) ;
                    if classe == -1 :
                        fragments,exon,real_class = constructTranscript(list(tuple(mRNA)),nb_exons, classe) ;
                        library += "\n".join(fragments) + "\n" ;
                        transcript += "\n".join(addClasse(mRNA,real_class)) ;
                        if mix :
                            mRNA_copy = [] ;
                            for ele_fea in list(mRNA) :
                                new = ele_fea.split("\t") ;
                                attr = new[-1].split(";") ;
                                attr[1] += ".1" ;
                                new[-1] = ";".join(attr) ;
                                mRNA_copy.append("\t".join(new)) ; 
                            library += "\n".join(mRNA_copy) + "\n" ;
                        transcript += "\n".join(addClasse(mRNA,real_class)) ;
                    elif classe == 0 :
                        library += "\n".join(addClasse(mRNA,classe)) + "\n"; 
                        transcript = "\n".join(addClasse(mRNA,classe)) ;
                    else :
                        fragments,introns,real_class = constructTranscript(list(tuple(mRNA)), nb_exons, classe) ;
                        transcript = "\n".join(fragments) ;
                        library += "\n".join(addClasse(mRNA,real_class)) + "\n"; 
                        if mix :
                            mRNA_copy = [] ;
                            for ele_fea in list(fragments) :
                                new = ele_fea.split("\t") ;
                                attr = new[-1].split(";") ;
                                attr[1] += ".1" ;
                                new[-1] = ";".join(attr) ;
                                mRNA_copy.append("\t".join(new)) ; 
                            library += "\n".join(mRNA_copy) + "\n" ;
                    reference += transcript +"\n";
                    if introns :
                        all_features += "\n".join(introns) + "\n" ; 
                    elif exon :
                        all_features += "\n".join(exon) + "\n" ;
                
                else :
                    ligne = gff.readline().rstrip() ;
                i += 1 ;
                
            else :
                ligne = gff.readline().rstrip() ;
                
    writeGFF(all_features,output_path + "_Features_of_interest") ;
    writeGFF(reference,output_path + "_reference_transcripts.tmp") ;
    writeGFF(library, output_path + "_library_transcripts.tmp") ;


def makeDensityLaw() :
    law = [] ;
    config_path = os.path.abspath(os.path.dirname(sys.argv[0]) + "/../config/intronStalker.properties")

    config = configparser.RawConfigParser() ;
    config.read(config_path) ;
    for classe in config["Density"] :
        effectif = int(config["Density"][classe]) ;
        law += [int(classe)]*effectif ;
    return tuple(law) ;


def parseGene(gene : list, nb_mRNA : int ,num_transcript : int) :
    # We randomly pick one mRNA
    choosen_mRNA = int(rd.choice(range(1,nb_mRNA+1),1)) ;
    j = 0 ;
    current_mRNA = [] ;
    nb_exons = 0 ;
    while j <= choosen_mRNA :
        ligne = gene.pop(0) ;
        if ligne.find("\tmRNA\t") != -1 or not gene:
            j += 1 ;
            old_mRNA = current_mRNA ;
            old_nb = nb_exons ;
            current_mRNA = [] ;
            nb_exons = 0 ;
        elif ligne.find("\texon\t") != -1  :
            nb_exons += 1 ;
            exon = changeAttributes(ligne, str(num_transcript), nb_exons) ;
            current_mRNA.append(exon) ;
    return old_mRNA, old_nb ;


def changeAttributes(exon : str, num_transcript : str, num_exon : int) :
    transcript_id = "T"+num_transcript ;
    exon_id = "exon"+str(num_exon) ;
    items = exon.split("\t") ;
    old_attributes = items[-1].split(";") ;
    new_attributes = "ID="+exon_id+";Parent="+transcript_id+";"+(";".join(old_attributes[2:])) ;
    return "\t".join(items[:-1]) + "\t" + new_attributes ;


def addClasse(mRNA : list, classe : int) :
    for e in range(len(mRNA)) :
        mRNA[e] += ";transcript_class="+str(classe) ;
    return mRNA ;


def constructTranscript(mRNA : list, nb_exons : int, classe : int) :
    features_interest=[] ;
    if classe == -1 and nb_exons > 2:
        spliced_exon = int(rd.choice(range(1,nb_exons-1),1)) ;
        mRNA,exon = spliceExon(mRNA,spliced_exon) ; 
        features_interest.append(exon) ;
    elif classe >= 1 :
        nb_introns = nb_exons -1 ;
        if nb_introns <= classe :
            retained_introns = range(nb_introns) ;
            classe = nb_introns ;
        else :
            retained_introns = rd.choice(range(nb_introns), classe,replace=False) ;
        # We'll take each retained intron in decreasing order to prevent exon index error
        for intron in sorted(retained_introns,reverse=True) :
            coord_start = calcCoord(mRNA[0:intron+1])
            mRNA[intron],feature_interest = keepIntron(mRNA[intron],mRNA[intron+1],coord_start) ;
            del mRNA[intron+1] ;
            features_interest.append(feature_interest) ;
    else :
        classe = 0 ;
    return addClasse(mRNA,classe),features_interest,classe ;


def keepIntron(exon1 : str, exon2 : str , coord_start : int) :
    items1 = exon1.split("\t") ;
    items2 = exon2.split("\t") ;
    
    # We check the strand 
    if items1[6] == "+" :
        # New exon's start is the first exon's start
        new_start = items1[3] ;
        intron_start = items1[4] ;
        # New exon's end is the second exon's end
        intron_end = items2[3] ;
        new_end = items2[4] ;
    else :
        # New exon's start is the second exon's start
        new_start = items2[3] ;
        intron_start = items2[4] ;
        # New exon's end is the first exon's end
        intron_end = items1[3] ;
        new_end = items1[4] ;
    
    intron_length = int(intron_end) - 1 - (int(intron_start) + 1) +1 ;
    coord_end = coord_start + intron_length -1

    # We update the new exon's attributes
    attributes1 = items1[-1].split(";") ;
    attributes2 = items2[-1].split(";") ;
    new_ID = attributes1[0] + "+" + attributes2[0].lstrip("ID=") ;
    new_attributes = ";".join([
        new_ID,
        attributes1[1],
        ";".join(attributes1[2:])
        ]) ;

    new_exon = "\t".join([
        items1[0],
        items1[1],
        "exon",
        new_start,
        new_end,
        ".",
        items1[6],
        ".",
        new_attributes
        ]) ;
    
    intron = "\t".join([
        attributes1[1].lstrip("Parent="),
        ".",
        "retained_intron",
        str(coord_start),
        str(coord_end),
        ".",
        items1[6],
        ".",
        ";".join([
            "ID=intron" + attributes1[0].lstrip("ID=exon"),
            "length="+str(intron_length)
        ])
        ]) ;
    
    return new_exon, intron ; 


def spliceExon(mRNA : list, spliced_exon : int) :
    feature = mRNA.pop(spliced_exon) ;
    items = feature.split("\t") ;
    attributes = items[-1].split(";") ;

    new_start = calcCoord(mRNA[0:spliced_exon]) ;
    exon_length = int(items[4]) - int(items[3]) +1
    new_end = int(new_start) + exon_length -1 ;

    exon = "\t".join([
        attributes[1].lstrip("Parent="),
        ".",
        "spliced_exon",
        str(new_start),
        str(new_end),
        ".",
        ".",
        ".",
        ";".join([attributes[0],"length="+str(exon_length)])
        ]) ;
    return mRNA,exon ;


def calcCoord(features_list : list) :
    total_length = 0 ;
    for feature in features_list :
        items = feature.split("\t") ;
        start,end = items[3],items[4] ;
        total_length += int(end) - int(start) +1 ;

    return total_length + 1 ;
    


def annoToData(annotation : str, fasta : str, nb : int, output : str, grinder : bool, mix : bool) :



    outdir = output + "_annoToData" ; # Output directory name (where this function will write 
                                 # all the results and tmp files).
    os.system("mkdir " + outdir) ;
    output_path = outdir + "/" + output ;

    choosen_genes = chooseGenes(annotation, int(nb)) ;
    generateTranscripts(annotation, choosen_genes, output_path, mix) ;
    extractFasta(fasta,output_path) ;

    if grinder :
        profile_path = os.path.abspath(os.path.dirname(sys.argv[0]) + "/../config/profile_file.txt")
        _grinder(output_path+"_library.fa", profile_path, output_path)
        split_read(output_path)
        os.system("rm {path}-reads.fa".format(path=output_path)) ;
        os.system("rm {path}_library.fa".format(path=output_path)) ;

