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
                    index = ["start_split","end_split","split_length","split_borders"])

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
    for ref_id in ref_id_list:
        aligned = bamfile.fetch(ref_id, multiple_iterators=True)
        split_reads=[]
        contig_seq = ref_dict.fetch(ref_id)
        for read in aligned:
            if read.cigartuples is not None:
            
                if '(3,' in str(read.cigartuples):
                    split = limit_from_cigar(read.cigartuples, read.reference_start, contig_seq)
                    split['read'] = read.query_name
                    split['reference'] = read.reference_name
                    if read.is_reverse :
                        split['strand'] = '-'
                    else :
                        split['strand'] = '+'
                    split_reads.append(split)
        df_split_reads = pd.DataFrame(split_reads)
        if not df_split_reads.empty :
            candidates.append(merge_split(df_split_reads,ref_id,contig_seq))
            split_alignments.append(df_split_reads)
    return pd.concat(candidates), pd.concat(split_alignments)

def merge_split(contig_reads,contig_name,contig_seq) :
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
        left_border = contig_seq[int(current.start_split): int(current.start_split) + 2]
        right_border = contig_seq[int(current.end_split) - 2: int(current.end_split)]
        candidates.append(pd.Series(
                name= '|'.join([
                    contig_name,
                    str(int(current.start_split)),
                    str(int(current.end_split))
                    ]),
                data= [
                    contig_name,
                    int(current.start_split),
                    int(current.end_split),
                    len(selected_reads),
                    left_border+'_'+right_border],
                index=['reference','start','end','depth','split_borders']
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
    
    # The assemblathon ouput will be named with the basename of the fasta file + '_saaemblathon.txt' as suffix
    assemblathon_name = output_path + os.path.splitext(os.path.basename(fastafile.name))[0] + '_assemblathon.txt'
    with open(assemblathon_name,'w') as assemblathon :
        sp.run(['assemblathon_stats.pl',fastafile.name],stdout=assemblathon)
    
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
    
    
    split_alignments=split_alignments[['reference','read','start_split','end_split','split_length','split_borders','strand']] #re-arrange the columns order of split_alignements output
    header_sa= list(split_alignments.columns.values)
    header_sa[0] = '#'+header_sa[0]
    split_alignments.to_csv(output_path+'_split_alignments.txt',header=header_sa,sep='\t',index=False)
    
    header_cand = ["#ID"] + list(candidates.columns.values)
    candidates.reset_index().to_csv(output_path+'_candidates.txt', header=header_cand, sep='\t', index=False)

#########################
# write truncated fasta #
#########################

def single_trim(candidate, sequences) :
    contig = sequences[candidate.reference]
    new_seq = contig[0:candidate.start].seq+contig[candidate.end:].seq
    new_record= SeqRecord(
        new_seq,
        id=contig.id+'::'+candidate.ID,
        name=contig.name+'::'+candidate.ID,
        description=contig.description+' trimmed_candidate='+candidate.ID
        )
    return new_record

def multi_trim(candidates,sequences) :
    contig = sequences[candidates.name]
    new_seq = contig.seq
    for i, row in candidates.sort_values('start',ascending=False).iterrows() :
        new_seq = new_seq[0:row.start]+new_seq[row.end:]
    new_record= SeqRecord(
        new_seq,
        id=contig.id+'.trimmed',
        name=contig.name+'.trimmed',
        description=contig.description+' trimmed_candidate='+','.join(candidates['ID'].values)
        )
    return new_record

def trimFastaFromTXT(reference, cand_file, output, prefix, force, multi) :
    """
    From a fasta file and a ggf file, write a new fasta file which contains expurgated sequences of provided features
    :param fasta_file: original sequences file
    :param gff_feature: gff file which contains the features to delete
    :param output: name of the output fasta file
    """
    
    print(type(cand_file.name))
    output_path = output + "/tf";
    if prefix:
        output_path += "_" + prefix;
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.mkdir(output)
    if not force:
        try :
            if os.path.exists(output_path + '_trimmed.fa'):
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    candidates=pd.read_csv(cand_file.name,sep='\t').rename(columns={'#ID':'ID'})
    sequences = SeqIO.to_dict(SeqIO.parse(reference.name,'fasta'))
    if not multi :
        trimmed_records = list(candidates.loc[lambda df : df.selected == 1].apply(single_trim,axis=1,sequences = sequences))
    else :
        trimmed_records = list(candidates.loc[lambda df : df.selected == 1].groupby('reference',sort=False).apply(multi_trim,sequences = sequences))
    
    
    no_trimmed_ids = set(sequences.keys()) - set(candidates.loc[lambda df : df.selected == 1,'reference'].values)
    no_trimmed_records=[]
    for c in no_trimmed_ids :
        no_trimmed_records.append(sequences[c])
    
    SeqIO.write(trimmed_records+no_trimmed_records,output_path+'_trimmed.fa','fasta')



#######################
# search ORF on fasta #
#######################


def df_ORF(pep_file,candidates) : 
    """
    From the TransDecoder output, write a tsv file which contains information about each predicted ORF.
    :param transdecoder_orfs: Transdecoder output which contains all the predicte ORFs (nucleic or proteic fasta format)
    :param output_tsv: name of the tsv ouput file
    :return: nothing
    """
    grep_out = sp.run(['grep','^>', pep_file],stdout=sp.PIPE)
    lines = grep_out.stdout.decode('utf-8').rstrip().split('\n')
    orfs=[]
    for line in lines :
        vals = line.lstrip('>').split()
        orf = {
            '#orf_id' : vals[0], #orf_id
            'reference' : vals[-1].split(':')[0], #reference
            'start' : int(vals[-1].split(':')[1].split('-')[0])-1, #start-1 because of O-based coord
            'end' : int(vals[-1].split(':')[1].split('-')[1].rstrip('(+-)')), #end 
            'length' : int(vals[4].split(':')[1]), #length
            'type' : vals[3].split(':')[1], #type
            'score' : float(vals[5].split('=')[1]), #score
            'strand' : vals[5].split(')')[0].lstrip('(')
        }
        # ~ print((candidates.reference==orf['reference'])&(~((candidates.end<orf['start'])|(candidates.start >= orf['end']))))
        c = candidates.loc[lambda c :(c.reference==orf['reference'])&~((c.end <= orf['start'])|(c.start >= orf['end'])),:]
        
        if not c.empty :
            print(orf)
            print(c)
            print()
        orfs.append(pd.Series(orf))
    return pd.DataFrame(orfs)
    
def analyzeORF(reference, cand_file, output, force, prefix, no_refine, rm) :
    
    candidates = pd.read_csv(cand_file.name,sep='\t')
    output_path = output + "/orf";
    if prefix:
        output_path += "_" + prefix;
        inter_dir = 'orf_'+prefix+'_intermediate'
    else :
        inter_dir = 'orf_intermediate'
    
    # Create output dir if not exist
    if not os.path.exists(output) :
        os.mkdir(output)
    if not force:
        try :
            if os.path.exists(output_path + '.txt'):
                   raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    cmd_long_orf = ['TransDecoder.LongOrfs','-t',os.path.abspath(reference.name),'-O',inter_dir]
    cmd_predict = ['TransDecoder.Predict','-t',os.path.abspath(reference.name),'-O',inter_dir]
    if no_refine :
        cmd_predict.append('--no_refine_starts')
        
    with open(output_path+'.log', 'w') as log :
        here = os.path.abspath('.')
        os.chdir(output)
        log.write('COMMANDS launched by intronSeeker :\n')
        log.write(' '.join(cmd_long_orf)+'\n')
        log.write(' '.join(cmd_predict)+'\n\n')
        long_prcs = sp.run(cmd_long_orf,stdout=log,stderr=sp.STDOUT)
        predict_prcs = sp.run(cmd_predict,stdout=log,stderr=sp.STDOUT)
        os.system('mv *.transdecoder.* {intermediate} ; mv pipeliner* {intermediate}'.format(intermediate=inter_dir))

    
    transdecoder_out = inter_dir+'/'+os.path.basename(reference.name)+".transdecoder.pep"
    orfs = df_ORF(transdecoder_out,candidates)
    
    if rm :
        os.system('rm -r '+inter_dir+'*')

    os.chdir(here)
    orfs.to_csv(output_path+'.txt',sep='\t',index=False)


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

