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
    split_start = start + limit_from_start[0] + 1
    split_end = start + limit_from_start[0] + limit_from_start[1]
    length = limit_from_start[1]
    flank_left = ref_seq[split_start - 1: split_start + 1]
    flank_right = ref_seq[split_end - 2: split_end]
    return pd.Series([int(split_start),int(split_end),length,flank_left+"_"+flank_right],
                    index = ["start_split","end_split","split_length","split_borders"])

def find_split(ref_id_list, bamfile, fastafile, mindepth, maxlen, minfootsize):
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
    for ref_id in ref_id_list:
        aligned = bamfile.fetch(ref_id, multiple_iterators=True)
        split_reads=[]
        contig_seq = ref_dict.fetch(ref_id)
        for read in aligned:
            #sarah
            if read.cigarstring is not None and "N" in read.cigarstring and "I" not in read.cigarstring:
                #print('cigar string', read.cigarstring)
                cigar_pattern = "([0-9]+)M([0-9]+)N([0-9]+)M"  
                cigar         = re.search(cigar_pattern, read.cigarstring)
                try:
                    cigarM1 = int(cigar.group(1))
                except:
                    cigarM1 = 0
                try : 
                    cigarM2 = int(cigar.group(3))
                except:
                    cigarM2 = 0   
                  # Remove reads with cigarM1 or cigarM2 < minfootsize       
                if min(cigarM1,cigarM2) < minfootsize:
                    foot = False
                    #print('READ TO REMOVE :', min(cigarM1,cigarM2), 'read: ', read.query_name,' and cigar', cigar )
                else:
                    foot = True    
      
            if read.cigartuples is not None and read.mapping_quality >= 2 and foot == True :
                if '(3,' in str(read.cigartuples):
                    split = limit_from_cigar(read.cigartuples, read.reference_start, contig_seq)
                    split['read'] = read.query_name
                    split['reference'] = read.reference_name
                    if read.is_reverse :
                        split['strand'] = '-'
                    else :
                        split['strand'] = '+'
                    split_reads.append(split)
                    #print('Keep ', foot, ' read : ', read.query_name, ' and cigar ', read.cigarstring)
        df_split_reads = pd.DataFrame(split_reads)
        if not df_split_reads.empty :
            contig_len = len(ref_dict.fetch(ref_id))
            candidates.append(merge_split(df_split_reads, ref_id, contig_seq, contig_len, mindepth, maxlen))
            split_alignments.append(df_split_reads)
    return pd.concat(candidates), pd.concat(split_alignments)

# Work per contig
def merge_split(contig_reads,contig_name,contig_seq,contig_len,mindepth,maxlen) :
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
        left_border  = contig_seq[int(current.start_split)-1: int(current.start_split) + 1]
        right_border = contig_seq[int(current.end_split) - 2: int(current.end_split)]
        
        #Flag "filter"
        flag = ""
        if len(selected_reads) <= mindepth:
            flag = "DP;"
        if ((current.end_split - current.start_split)/contig_len)*100 > maxlen:
            flag += "LEN;"
        if left_border+'_'+right_border != "CT_AC" and left_border+'_'+right_border != "GT_AG":
            flag += "SS;"
        if flag == "":
            flag = "PASS"
        else:
            flag = flag[:-1]
        
        candidates.append(pd.Series(
            name = '|'.join([
                        contig_name,
                        str(int(current.start_split)),
                        str(int(current.end_split))
                    ]),
            data = [
                    contig_name,
                    int(current.start_split),
                    int(current.end_split),
                    len(selected_reads),
                    left_border+'_'+right_border,
                    flag
                    ],
            index =  ['reference','start','end','depth','split_borders','filter']
        ))
        split_alignments = split_alignments.drop(selected_reads.index).reset_index(drop=True)
    return pd.DataFrame(candidates)

def splitReadSearch(bamfile, fastafile, mindepth, maxlen, output, prefix, force, threads, minfootsize) :
    """
    Search the split reads and write two output files : the first with all the spliced events and the second where all the identical spliced events are merged
    :param bamfilename: name of the input alignment file
    :param fastafilename: name of the reference fasta file (contains the contigs) 
    :param basename: prefix for the name of output files 
    :return: nothing
    """
    output_path = output + "/srs"
    if prefix:
        output_path += "_" + prefix
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.makedirs(output)
    if not force:
        try :
            if os.path.exists(output_path + "_candidates.txt") or os.path.exists(output_path + "_split_alignments.txt") :
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    ref_id_list = [x.split("\t")[0] for x in pysam.idxstats(bamfile.name).split("\n")[:-2]]
    
    if not os.path.exists(fastafile.name + '.fai'):
        pysam.faidx(fastafile.name)

    with prl.ProcessPoolExecutor(max_workers=threads) as ex :
        # ~ s_t = time.time()
        ref_id_array = np.array_split(ref_id_list,ex._max_workers)
        out = list(zip(*list(ex.map(
            find_split,
            ref_id_array,
            repeat(bamfile.name,ex._max_workers),
            repeat(fastafile.name,ex._max_workers),
            repeat(mindepth,ex._max_workers),
            repeat(maxlen,ex._max_workers),
            repeat(minfootsize,ex._max_workers)
            ))))
        candidates= pd.concat(out[0])
        split_alignments = pd.concat(out[1])
        
        # ~ print(time.time()-s_t)

    # After filtering, focus on flagged PASS candidates and modify
    # filter field (from PASS to OI) if two "retained introns"
    # are overlapping
    prevCtg   = ""
    prevStart = 0
    prevEnd   = 0
    prevFilter= ""
    prevIndex = 0
    for k, v in candidates.iterrows():
        if "PASS" in str(v['filter']): 
            if (prevCtg == "" or prevCtg != v['reference']):  
                prevCtg   = v['reference']
                prevStart = v['start']
                prevEnd   = v['end']
                prevFilter= v['filter']
                prevIndex = k
            elif (v['start'] < prevEnd):
                candidates.loc[k, 'filter'] = "OI"
                candidates.loc[prevIndex, 'filter'] = "OI" 
            prevCtg   = v['reference']
            prevStart = v['start']
            if prevEnd > v['end']:
                prevEnd = v['end']
            prevFilter= v['filter']
            prevIndex = k
    # ~ candidates = find_split(ref_id_list,bamfile.name,fastafile.name)
    # ~ candidates['selected'] = 1
    # ~ print(candidates)
    
    split_alignments=split_alignments[['reference','read','start_split','end_split','split_length','split_borders','strand']] #re-arrange the columns order of split_alignments output
    header_sa= list(split_alignments.columns.values)
    header_sa[0] = '#'+header_sa[0]
    split_alignments.to_csv(output_path+'_split_alignments.txt',header=header_sa,sep='\t',index=False)
    
    f = open(output_path+'_candidates.txt', 'w')
    f.write('##mindepth:' + str(mindepth) + '\n')
    f.write('##maxlen:'   + str(maxlen) + '\n') 
    header_cand = ["#ID"] + list(candidates.columns.values)
    candidates.reset_index().to_csv(f, header=header_cand, sep='\t', index=False)
    f.close()

#########################
# write truncated fasta #
#########################

def single_trim(candidate, sequences) :
    contig = sequences[candidate.reference]
    new_seq = contig[0:candidate.start-1].seq+contig[candidate.end:].seq
    new_record= SeqRecord(
        new_seq,
        id=contig.id+'##'+candidate.ID,
        name=contig.name+'##'+candidate.ID,
        description=contig.description+' trimmed_candidate='+candidate.ID
        )
    return new_record

def multi_trim(candidates,sequences) :
    contig = sequences[candidates.name]
    new_seq = contig.seq
    for i, row in candidates.sort_values('start',ascending=False).iterrows() :
        new_seq = new_seq[0:row.start-1]+new_seq[row.end:]
    new_record= SeqRecord(
        new_seq,
        id=contig.id+'.trimmed',
        name=contig.name+'.trimmed',
        description=contig.description+' trimmed_candidate='+','.join(candidates['ID'].values)
        )
    return new_record

def trimFastaFromTXT(reference, cand_file, output, prefix, force, multi) :
    """
    From a fasta file and a candidates file, write a new fasta file which contains expurgated sequences of provided features
    :param fasta_file: original sequences file
    :param gff_feature: candidates file which contains the features to delete
    """
    output_path = output + "/tf";
    if prefix:
        output_path += "_" + prefix;
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.makedirs(output)
    if not force:
        try :
            if os.path.exists(output_path + '_trimmed.fa'):
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    candidates = pd.read_csv(cand_file.name,sep='\t',skiprows=2).rename(columns={'#ID':'ID','filter':'Filter'})
    sequences  = SeqIO.to_dict(SeqIO.parse(reference.name,'fasta'))
    
    trimmed_records=[]
    if candidates.loc[lambda df : df.Filter == "PASS"].size :
        if not multi :
            trimmed_records = list(candidates.loc[lambda df : df.Filter == "PASS"].apply(single_trim,axis=1,sequences = sequences))
        else :
            trimmed_records = list(candidates.loc[lambda df : df.Filter == "PASS"].groupby('reference',sort=False).apply(multi_trim,sequences = sequences))
    
    no_trimmed_ids = set(sequences.keys()) - set(candidates.loc[lambda df : df.Filter == "PASS",'reference'].values)
    no_trimmed_records=[]
    for c in no_trimmed_ids :
        no_trimmed_records.append(sequences[c])
    SeqIO.write(trimmed_records+no_trimmed_records,output_path+'_trimmed.fa','fasta')



#######################
# search ORF on fasta #
#######################

def df_lg_ref_ORF(getorf_file, candidates) :
    """
    Search longest getorf overlapping start of each candidate
    Return new df candidates with one more col "lg_ref_ORF"
    """
    grep_out = sp.run(['grep','^>', getorf_file],stdout=sp.PIPE)
    lines = grep_out.stdout.decode('utf-8').rstrip().split('\n')
    # orfs = dict of array of found orfs by sequence
    orfs={}
    for line in lines :
        vals = line.lstrip('>').split()
        tmp_id = vals[0].split("_")
        tmp_id.pop()
        id    = "_".join(tmp_id)
        start = int(vals[1].lstrip('['))
        end   = int(vals[3].rstrip(']'))
        if id not in orfs:
            orfs[id] = []
        if start<end:
            orfs[id].append([start, end])
        else:
            orfs[id].append([end, start])
    df_new_candidates = candidates[candidates["filter"] == "PASS"].copy()
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    #print('New_cand:', df_new_candidates.head(5), '\n\n')
    longestORF = []
    for index, row in df_new_candidates.iterrows(): 
        current_longestORF = 0
        if row["reference"] in orfs:
            for s,e in orfs[row["reference"]]:
                #print("Deb:",s," End:",e)
                if(row["start"]>=s and row["start"]<=e):
                    #print("Deb:",s," End:",e, "<======DANS LA ZONE")
                    t = e-s+1
                    if(t>current_longestORF):
                        current_longestORF = t
        longestORF.append(current_longestORF) 
        #print(row["#ID"], " ",row["reference"]," ", str(row["start"])," ",str(row["end"]))
    df_new_candidates['lgORF_ref'] = longestORF
    print('New_cand:', df_new_candidates.head(5), '\n\n')
    return df_new_candidates

def df_lg_trim_ORF(getorf_file, candidates) :
    """
    Search longest getorf overlapping start of each candidate
    Return new df candidates with one more col "lg_ref_ORF"
    """
    grep_out = sp.run(['grep','^>', getorf_file],stdout=sp.PIPE)
    lines = grep_out.stdout.decode('utf-8').rstrip().split('\n')
    # orfs = dict of array of found orfs by sequence
    orfs={}
    for line in lines :
        vals = line.lstrip('>').split()
        id    = vals[-1].split('=')[1]
        start = int(vals[1].lstrip('['))
        end   = int(vals[3].rstrip(']'))
        if id not in orfs:
            orfs[id] = []
        if start<end:
            orfs[id].append([start, end])
        else:
            orfs[id].append([end, start])
    df_new_candidates = candidates[candidates["filter"] == "PASS"].copy()
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    longestORF = []
    for index, row in df_new_candidates.iterrows(): 
        current_longestORF = 0
        if row["#ID"] in orfs:
            for s,e in orfs[row["#ID"]]:
                if(row["start"]>=s and row["start"]<=e):
                    t = e-s+1
                    if(t>current_longestORF):
                        current_longestORF = t
        longestORF.append(current_longestORF)
        
    df_new_candidates['lgORF_trim'] = longestORF
    print('New_cand:', df_new_candidates.head(5), '\n\n')
    return df_new_candidates

def df_parseDiamond(diamond_file, candidates) :
    file = open(diamond_file, "r")
    line = file.readline()
    prots={}
    while line:
        vals = line.split('\t')
        id    = vals[0]
        start = int(vals[2])
        end   = int(vals[3])
        if id not in prots:
            prots[id] = []
        if start<end:
            prots[id].append([start, end])
        else:
            prots[id].append([end, start]) 
        line = file.readline()
    file.close()
 
    protoverlap = []
    protbefore = []
    protafter = []
    protbilan = []
    for index, row in candidates.iterrows(): 
        nboverlap = 0
        nbbefore = 0
        nbafter = 0

        #ENSGALT00000097829.modif|478|1157      1       441     405     12      12      1       --

        if row["reference"] in prots:
            for s,e in prots[row["reference"]]:
                #if(s+10<=row["end"] and e-10>=row["start"]):
                #    nboverlap += 1
                
                #if((row["start"]+15>s and row["start"]+15<e) or ((row["end"]-15>s and row["end"]-15<e))):
                #    nboverlap += 1
                center = row["start"] + ((row["end"]-row["start"]+1)/2)
                if(s<=center and e>=center):
                    nboverlap += 1
                elif(abs(row["start"]-e) < 10):
                    nbbefore += 1
                elif(abs(s-row["end"]) < 10):
                    nbafter += 1
        
        protoverlap.append(nboverlap)
        protbefore.append(nbbefore)
        protafter.append(nbafter)
        if(nbafter+nbbefore > nboverlap):
            protbilan.append("OK")
        else:
            protbilan.append("--")
        
    candidates['protoverlap'] = protoverlap
    candidates['before'] = protbefore
    candidates['after'] = protafter
    candidates['bilan'] = protbilan

    print('New_cand:', candidates.head(5), '\n\n')
    return candidates


def analyzeORF(reference, trim_ref, db_file, cand_file, output, force, prefix, no_refine, rm) :
    candidates = pd.read_csv(cand_file.name,sep='\t',skiprows=2)
    output_path = output + "/orf"
    if prefix:
        output_path += "_" + prefix
    output_tmp_dir = output_path+'_tmp'
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.makedirs(output)
    if not os.path.exists(output_tmp_dir) :
        os.makedirs(output_tmp_dir)
    if not force:
        try :
            if os.path.exists(output_path + '.txt'):
                raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    cmd_getorf = ['getorf','-sequence',os.path.abspath(reference.name),'-outseq',output_tmp_dir+'/reference.getorf']
    cmd_getorf_trim = ['getorf','-sequence',os.path.abspath(trim_ref.name),'-outseq',output_tmp_dir+"/trimmed_reference.getorf"]

    with open(output_path+'.log', 'w') as log :
        log.write('COMMANDS launched by intronSeeker :\n')
        log.write(' '.join(cmd_getorf)+'\n')
        getorf_prcs = sp.run(cmd_getorf,stdout=log,stderr=sp.STDOUT)
        log.write(' '.join(cmd_getorf_trim)+'\n\n')
        getorf_trim_prcs = sp.run(cmd_getorf_trim,stdout=log,stderr=sp.STDOUT)
    
    df_candidates_new = df_lg_ref_ORF(output_tmp_dir+'/reference.getorf', candidates)
    df_candidates_new = df_lg_trim_ORF(output_tmp_dir+'/trimmed_reference.getorf', df_candidates_new)

    if rm :
        os.system('rm -r '+output_tmp_dir)

    #df_candidates_new.to_csv(output_path+'.txt',sep='\t',index=False)
    #orfs_trim.to_csv(output_path+'_trim.txt',sep='\t',index=False)

    ## Prot
    db = output_tmp_dir+"/"+os.path.basename(db_file.name)
    cmd_diamond_makedb = ['diamond','makedb','--in',os.path.abspath(db_file.name),'--db',db]
    makedb_prcs = sp.run(cmd_diamond_makedb)

    os.system("diamond blastx -q {fasta} -d {db} -o {output} -f 6 qseqid qlen qstart qend sseqid slen sstart send pident length evalue bitscore -p 1 -e 0.05 --max-target-seqs 0".format(
        fasta=os.path.abspath(reference.name),db=db+".dmnd",output=db+"_ref_output.tsv"))
    print(db+"_ref_output.tsv")
    df_candidates_new = df_parseDiamond(db+"_ref_output.tsv", df_candidates_new)
    df_candidates_new.to_csv(output_path+'.txt',sep='\t',index=False)


    #os.system("diamond blastx -q {fasta} -d {db} -o {output} -f 6 qseqid qlen sseqid slen pident length qstart qend evalue bitscore -p 1 -e 0.05 --max-target-seqs 0".format(
    #    fasta=os.path.abspath(trim_ref.name),db=db+".dmnd",output=db+"_trim_output.tsv"))

    
    # `diamond makedb --in $diamond --db $outputDir/private_tmp_ 2> $outputDir/private.diamond_makedb_tmp_.stderr`;
    # `diamond blastx --db $outputDir/private_tmp_ --query $prokkaTFA --threads $threads --outfmt 6 qseqid sseqid pident nident length mismatch gaps gapopen qstart qend sstart send evalue bitscore stitle qcovhsp  -k 1000 --query-cover $diamond_queryCover -e $diamond_evalue --out $outputDir/private.diamond.tsv 2> $outputDir/private.diamond_tmp_.stderr`;


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
    os.system("makedirs " + outdir) ;
    path_to_output = outdir + "/" + output ;
    
    runDiamond(fasta, dbprotein, path_to_output, threads) ;
    parseDiamondXML(path_to_output+".tmp.xml", path_to_output) ;
    if not rm :
        os.system("rm {outdir}/*.tmp.*".format(outdir=outdir)) ;

