#!/usr/bin/env python3

# Modules
try :
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Blast import NCBIXML as bx ;
    from collections import defaultdict
    from itertools import repeat
    import concurrent.futures as prl
    import pandas as pd
    import numpy as np
    import pysam
    import os
    import re
    import subprocess as sp
    import time
except ImportError as error :
    print(error) ;
    exit(1) ;

def change_key(dictionary: dict, key_to_clean: dict):
    """
    Change key dict name

    :param dictionary: dictionary to change keys
    :param key_to_clean: dictionary of new keys associate to previous keys
    :return: clean dictionary
    """
    for clef, value in key_to_clean.items():
        if clef in dictionary:
            dictionary[value] = dictionary.pop(clef)
    return dictionary


#######################
# Extract split reads #
#######################

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

def find_split(ref_id_list, bamfile, fastafile):
    """
    For an alignment file, list all split reads.

    :param fastafilename: reference fasta file
    :param bamfilename: AlignmentFile object with all reads alignment information
    :return: list of split. For each split, save its reference, name, start, stop, length and flanking sequences.
    """
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    try :
        ref_dict = pysam.FastaFile(fastafile)
    except OSError as e :
        if str(e).startswith("error when opening file") :
            print("\nIndexingError : it's impossible to write in Reference fasta directory to indexing file.")
            print("Please move the reference file in a directory where it's possible to write in and retry.\n")
        else :
            print(e)
        exit(1)
    
    split_alignments = []
    candidates=[]
    for id in ref_id_list:
        aligned = bamfile.fetch(id, multiple_iterators=True)
        split_reads=[]
        for read in aligned:
            reference = ref_dict.fetch(id)
            if read.cigartuples is not None:
            
                if '(3,' in str(read.cigartuples):
                    split = limit_from_cigar(read.cigartuples, read.reference_start, reference)
                    split['read'] = read.query_name
                    split['contig'] = read.reference_name
                    if read.is_reverse :
                        split['strand'] = '-'
                    else :
                        split['strand'] = '+'
                    split_reads.append(split)
        df_split_reads = pd.DataFrame(split_reads)
        if not df_split_reads.empty :
            candidates.append(merge_split(df_split_reads,id))
            split_alignments.append(df_split_reads)
        
    return pd.concat(candidates), pd.concat(split_alignments)

def merge_split(contig_reads,contig_name) :
    candidates = []
    split_alignments = contig_reads.sort_values(by=['start_split','end_split']).reset_index(drop=True)
    while not split_alignments.empty :
        current = split_alignments.loc[0,['start_split','end_split']] 
        split_alignments['start_split']=pd.to_numeric(split_alignments['start_split'])
        split_alignments['end_split']=pd.to_numeric(split_alignments['end_split'])
        selected_reads = split_alignments.query(
            '(-5 <= start_split-{current_start} <= 5) and (-5 <= end_split-{current_end} <= 5)'.format(
                current_start=float(current.start_split),
                current_end=float(current.end_split)
                )
            )
        old_size = 1
        while len(selected_reads) != old_size :
            current.start_split = round(selected_reads['start_split'].mean())
            current.end_split = round(selected_reads['end_split'].mean())
            old_size = len(selected_reads)
            selected_reads = split_alignments.query(
                '(-5 <= start_split-{current_start} <= 5) and (-5 <= end_split-{current_end} <= 5)'.format(
                    current_start=current.start_split,
                    current_end=current.end_split
                    )
                )
        candidates.append(pd.Series(
                name= '|'.join([
                    contig_name,
                    str(int(current.start_split)),
                    str(int(current.end_split))
                    ]),
                data= [contig_name,int(current.start_split),int(current.end_split),len(selected_reads)],
                index=['contig','start','end','depth']
                ))
        split_alignments = split_alignments.drop(selected_reads.index).reset_index(drop=True)
    return pd.DataFrame(candidates)

def splitReadSearch(bamfile, fastafile, output, prefix, force, threads) :
    """
    Search the split reads and write two output files : the first with all the spliced events and the second where all the identical spliced events are merged
    :param bamfilename: name of the input alignment file
    :param fastafilename: name of the reference fasta file (contains the contigs) 
    :param basename: prefix for the name of output files 
    :return: nothing
    """
    output_path = output + "/srs";
    if prefix:
        output_path += "_" + prefix;
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.mkdir(output)
    if not force:
        try :
            if os.path.exists(output_path + "_candidates.txt") or os.path.exists(output_path + "_split_alignments.txt") :
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    ref_id_list = [x.split("\t")[0] for x in pysam.idxstats(bamfile.name).split("\n")[:-2]]

    with prl.ProcessPoolExecutor(max_workers=threads) as ex :
        # ~ s_t = time.time()
        ref_id_array = np.array_split(ref_id_list,ex._max_workers)
        out = list(zip(*list(ex.map(
            find_split,
            ref_id_array,
            repeat(bamfile.name,ex._max_workers),
            repeat(fastafile.name,ex._max_workers)
            ))))
        candidates= pd.concat(out[0])
        split_alignments = pd.concat(out[1])
        candidates['selected'] = 1
        # ~ print(time.time()-s_t)
    
    # ~ candidates = find_split(ref_id_list,bamfile.name,fastafile.name)
    # ~ candidates['selected'] = 1
    # ~ print(candidates)
    
    header_sa= list(split_alignments.columns.values)
    header_sa[0] = '#'+header_sa[0]
    split_alignments.to_csv(output_path+'_split_alignments.txt',header=header_sa,sep='\t',index=False)
    
    header_cand = ["#ID"] + list(candidates.columns.values)
    candidates.reset_index().to_csv(output_path+'_candidates.txt', header=header_cand, sep='\t', index=False)

#########################
# write truncated fasta #
#########################

def _Truncate(fasta_file : str, gff_feature : str, output_path : str) :
    print("\nTruncate FASTA file : "+fasta_file+"\n");
    with open(fasta_file,"r") as input_f :
        list_seqs = input_f.read().lstrip(">").rstrip().split("\n>") ;
        seqs = { seq.split("\n",1)[0].split()[0] : seq.split("\n",1)[1].replace("\n","") for seq in list_seqs} ;
        with open(gff_feature,"r") as features_f :
            list_features = features_f.read().rstrip().split("\n") ;
            fasta_truncated = "" ;
            num_feature = 1 ;
            old_contig = "" ;
            for feature in list_features :
                items = feature.split("\t") ;
                if old_contig != items[0] :
                    num_feature = 1 ;
                else :
                    num_feature += 1 ;
                seq_truncated = seqs[items[0]][0:int(items[3])] + seqs[items[0]][int(items[4]):]+"\n"
                seq_head = ">"+items[0]+".e"+str(num_feature)+"\t"+items[0]+"\t"+items[3]+"-"+items[4]+"\t"+items[-1]+"\n"
                print(seq_head)
                fasta_truncated += seq_head + seq_truncated ;
                old_contig = items[0] ;
        with open(output_path + "-trq.fa","w") as fa :
            fa.write(fasta_truncated.rstrip()) ;

def truncate(fasta_file : str, gff_feature : str, output : str) :
    """
    From a fasta file and a ggf file, write a new fasta file which contains expurgated sequences of provided features
    :param fasta_file: original sequences file
    :param gff_feature: gff file which contains the features to delete
    :param output: name of the output fasta file
    """
    outdir = output + "_truncate" ;
    os.system("mkdir " + outdir) ;
    output_path = outdir + "/" + output ;
    
    _Truncate(fasta_file,gff_feature,output_path) ;



#######################
# search ORF on fasta #
#######################

def runTransDecoder(fasta_file : str, tmp_directory : str, refine : bool) :
    """
    Run TransDecoder.
    :param fasta_file: Fasta file on which the ORF prediction will be perform
    :param output_directory: directory name where the TransDecoder files will be stored
    :return: nothing
    """
    if refine :
        os.system("TransDecoder.LongOrfs -t {fasta} -O {tempdir} ; TransDecoder.Predict -t {fasta} -O {tempdir}".format(fasta=fasta_file,tempdir=tmp_directory)) ;
    else :
        os.system("TransDecoder.LongOrfs -t {fasta} -O {tempdir} ; TransDecoder.Predict -t {fasta} -O {tempdir} --no_refine_starts".format(fasta=fasta_file,tempdir=tmp_directory)) ;

def writeORFasTsv(transdecoder_orfs : str) : 
    """
    From the TransDecoder output, write a tsv file which contains information about each predicted ORF.
    :param transdecoder_orfs: Transdecoder output which contains all the predicte ORFs (nucleic or proteic fasta format)
    :param output_tsv: name of the tsv ouput file
    :return: nothing
    """
    print("Write TSV file from "+transdecoder_orfs) ;
    tsv="\t".join(["Sequence","Start","End","ORFid","length","type","strand","frame","score"]) ;
    grep_out = sp.getoutput("grep 'CDS' {filename} | cut -f 8,9".format(filename=transdecoder_orfs.replace(".pep",".gff3"))) ;
    frames = { cds.split("Parent=")[1] : cds.split("\t")[0] for cds in grep_out.split("\n") }
    with open(transdecoder_orfs,"r") as ORFfile :
        lines = [orf.split("\n",1)[0] for orf in ORFfile.read().lstrip(">").split("\n>")] ;
        for line in lines :
            vals = line.split() ;
            ORFid = vals[0] ; typ=vals[3].lstrip("type:") ; length=vals[4].lstrip("len:") ;
            strand= vals[5].split(",")[0].lstrip("(").rstrip(")") ;
            score = vals[5].split(",")[1].lstrip("score=") ;
            seq = vals[-1].split(":")[0] ;
            start = vals[-1].split(":")[1].split("-")[0] ; end = vals[-1].split(":")[1].split("-")[1].split("(")[0] ;
            tsv += "\n"+"\t".join([seq,start,end,ORFid,length,typ,strand,frames[ORFid],score]) ;
    with open(transdecoder_orfs+".cds.tsv","w") as tsvfile :
        tsvfile.write(tsv) ;

def predictORF(fasta_file, trunc, gff, rm, output : str,refine) :

    outdir = output + "_orf" ;
    os.system("mkdir " + outdir) ;
    output_path = outdir + "/" + output ;
    fasta_name = os.path.basename(fasta_file)
    

    runTransDecoder(fasta_file,output_path + "_intermediate",refine) ;
    os.system("mv *.transdecoder.* {outdir}".format(outdir=outdir)) ;
    writeORFasTsv(outdir +"/"+ fasta_name + ".transdecoder.pep") ;

    if trunc and gff:
        _Truncate(fasta_file, gff, output_path) ;
        runTransDecoder(output_path+"-trq.fa",output_path + "_intermediate-tronque",refine) ;
        os.system("mv *.transdecoder.* {outdir}".format(outdir=outdir)) ;
        writeORFasTsv(output_path+"-trq.fa.transdecoder.pep") ;

    if rm :
        os.system("rm pipeliner* ; rm -r {outdir}_intermediate* ".format(outdir=output_path))


#################################################
# Search long gap on contigs-proteins alignment #
#################################################

def runDiamond(fasta_file : str, protein_db : str, output : str,threads : int) :
    """
    Run Diamond v.0.9.9 to align nucleic sequences against proteic database.
    :param fasta_file: Filename of the fasta which contains the nucleic sequences.
    :param protein_db: Name of the Diamond proteic database.
    :param output: Name of the output file.
    :param threads: Number of threads used for the alignment.
    """
    os.system("diamond blastx -q {fasta} -d {proteindb} -o {output} -f 5 -p {threads} -e 0.05 --max-target-seqs 0".format(
        fasta=fasta_file,proteindb=protein_db, output=output+".tmp.xml",threads=threads)) ;

def searchGapsOfInterest(alignseq_subject : str , alignseq_query : str) :
    """
    From query eand susject aligned sequences (extracted from Diamond XML output), find long gaps and return a list of them.
    :param alignseq_subject: String which contains the subject sequence.
    :param alignseq_query: String which contains the query sequence.
    :return list_gaps: List of long Gaps found in the alignment.
    """
    gap_pattern = re.compile("-+") ;
    begin_idx = 0 ;
    gaps = gap_pattern.findall(alignseq_subject)
    list_gaps = [] ;
    if len(gaps) < 4 and len(gaps) > 0 :
        for gap_str in gap_pattern.findall(alignseq_subject) :
            if len(gap_str) > 5 :
                gap_start = alignseq_subject.index(gap_str,begin_idx) ;
                gap_length = len(gap_str) 
                
                # We want the gap coordinates on query sequence, so we calculate the total length of all the query gaps located before the gap of interest 
                translation = len("".join(gap_pattern.findall(alignseq_query,0,gap_length)))
                query_AA_start = gap_start-translation ;
                query_AA_end = gap_start-translation+gap_length-1 ;
                list_gaps.append((query_AA_start,query_AA_end)) ;
        return list_gaps

def parseDiamondXML(diamXml : str,output : str) :
    """
    Parse Diamond XML output and write a GFF file which contains the Long gaps.
    :param diamXml: Name of the Diamond ouput file (XML format).
    :param output: Name of the ouput file.
    """
    print("\nParsing of Diamond XML output to produce Long Gaps GFF \n")
    gff = [] ;
    align_gff = [] ;
    with open(diamXml) as resDiamond:
        blast_records = bx.parse(resDiamond) ;
        for record in blast_records :
            for align in record.alignments :
                num_hsp = 1 ;
                for hsp in align.hsps :
                    contig = record.query ;
                    align_start = hsp.query_start ;
                    align_end = hsp.query_end ;
                    frame = hsp.frame[0] ;
                    if frame > 0 :
                        strand = "+" ;
                    else :
                        strand = "-" ;
                    prot_align = align.title.split()[0]+".hsp"+str(num_hsp) ;
                    ID_align = contig.split()[0] +"_"+align.title.split()[0]+".hsp"+str(num_hsp) ;
                    features_align = ";".join(["ID="+ID_align,"PROT="+prot_align,"LEN="+str(hsp.align_length),"BITSCORE="+str(hsp.bits),"EVAL="+str(hsp.expect)])
                    align_gff.append("\t".join(map(str,[contig.split()[0],"Diamond","Alignment",align_start,align_end,".",strand,frame,features_align])))

                    list_gaps = searchGapsOfInterest(hsp.sbjct,hsp.query)
                    if list_gaps :
                        num_gap = 1
                        for gap in list_gaps :
                            # We multiply the gap coordinates 3 times because we want the cordinates on nucleic sequence.
                            start = gap[0]*3 ;
                            end = gap[1]*3 ;
                            prot = align.title.split()[0]+".gap"+str(num_gap) ;
                            ID = contig.split()[0] +"_"+align.title.split()[0]+".gap"+str(num_gap) ;
                            features = ";".join(["ID="+ID,"PROT="+prot,"LEN="+str(end-start),"BITSCORE="+str(hsp.bits),"EVAL="+str(hsp.expect)])
                            gff.append("\t".join(map(str,[contig.split()[0],"Diamond","Gap",start,end,".",strand,frame,features])))
                            num_gap += 1 ;
                    num_hsp += 1
    with open(output+".long_gaps.gff","w") as ofile :
        ofile.write("\n".join(gff)+"\n") ;
    with open(output+".alignments.gff","w") as o2file :
        o2file.write("\n".join(align_gff)+"\n") ;

def searchProtein(fasta : str, dbprotein : str, output : str, threads : int, rm = False) :
    """
    Run all the Protein analysis to search potential introns (i.e. long gaps beetween contigs and proteins
    :param fasta: Name of the fasta file which contains contigs.
    :param dbprotein: Name of the Diamond proteic database.
    :param output: Name of the ouput file.
    :param threads: Number of threads used for the alignment.
    :param rm : Boolean which rules the temporary files erasure.
    """
    outdir = output + "_protein" ;
    os.system("mkdir " + outdir) ;
    path_to_output = outdir + "/" + output ;
    
    runDiamond(fasta, dbprotein, path_to_output, threads) ;
    parseDiamondXML(path_to_output+".tmp.xml", path_to_output) ;
    if not rm :
        os.system("rm {outdir}/*.tmp.*".format(outdir=outdir)) ;

