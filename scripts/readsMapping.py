#!/usr/bin/env python3

try :
    import os
    import pysam
    import subprocess as sp
    import sys
except ImportError as error :
    print(error)
    exit(1)

def star_bam_fix(input_file: str, output_file):
    """
    Rename each read from STAR alignment according to which mate it is (mate 1 or mate 2)
    :param input_file: name of the bam which contains STAR alignment 
    :param output_file: name of the ouput file 
    :return: nothing
    """
    bamfile = pysam.AlignmentFile(input_file, "rb",check_sq=False)
    outfile = pysam.AlignmentFile(output_file, "wb", template=bamfile)
    for record in bamfile.fetch(until_eof=True):
        if record.is_read1:
            record.query_name += "/1"
        if record.is_read2:
            record.query_name += "/2"
        outfile.write(record)

def bam_indexing(bamfile, cmdlog) :
    bam_index_cmd = ['samtools','index',
                    bamfile,
                    bamfile.rstrip('.bam') + '.bai'
                    ]
    cmdlog.write("\n# BAM indexing:\n") 
    cmdlog.write(' '.join(bam_index_cmd) + "\n")
    sp.run(bam_index_cmd)

def flagstat(bamfile, cmdlog, threads = 1) :
    bam_flagstat_cmd = ['samtools','flagstat',
                        '--threads',str(threads),
                        bamfile
                        ]
    cmdlog.write("\n# BAM flagstat computation:\n") 
    cmdlog.write(' '.join(bam_flagstat_cmd) + "\n")
    stdout = sp.run(bam_flagstat_cmd,stdout=sp.PIPE).stdout
    with open(bamfile.rstrip('.bam')+'.flagstat.txt','w') as out_flag :
        out_flag.write(stdout.decode('utf-8'))

def star(reference, r1, r2, output, prefix, force, rm, threads):
    """
    Run STAR.
    Call the aligner STAR. 
    :param ref: reference fasta file
    :param r1: reads files. If it's paired-end, precise the second reads files with r2
    :param r2: second read file for paired-end
    :param pref: output prefix
    :return:
    """
    output_path = output + "/star";
    if prefix:
        output_path += "_" + prefix;
    
    genomedir = output + "/star_genomeRef" # Path to the directory will contain reference fasta file and index files 
    
    # Create output dir if not exist
    if not os.path.exists(genomedir) :
        os.makedirs(genomedir)
    if not force:
        try :
            if os.path.exists(output_path + ".sort.bam") or os.path.exists(output_path + ".sort.flagstat.txt"):
               raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    cmdlog = open(output_path+".log", "w")
    
    # Reference File indexing
    index_command = ['STAR',
                    '--runMode', 'genomeGenerate',
                    '--runThreadN', str(threads),
                    '--genomeSAindexNbases', str(6),
                    '--genomeDir', genomedir, 
                    '--genomeFastaFiles', reference.name,
                    '--outFileNamePrefix', output_path+'.star.'
                    ]
    cmdlog.write('\n# Fasta indexing:\n')
    cmdlog.write(" ".join(index_command) + "\n")
    sp.run(index_command) 
    
    # ~ # Reads Mapping and ouput files writing 
    star_command = ['STAR',
                    '--genomeDir', genomedir,
                    '--runThreadN', str(threads),
                    '--outSAMstrandField', 'intronMotif',
                    '--outSAMtype', 'BAM', 'SortedByCoordinate',
                    '--outSAMunmapped', 'Within',
                    '--outFileNamePrefix', output_path+'.star.'
                    ]
    if r2 :
        star_command.extend(['--readFilesIn', r1.name, r2.name])
    else :
        star_command.extend(['--readFilesIn', r1.name])
    if r1.name.endswith('.gz') :
        star_command.extend(['--readFilesCommand', 'zcat'])
    cmdlog.write('\n# STAR alignement:\n')
    cmdlog.write(' '.join(star_command) + "\n")
    sp.run(star_command)
    star_bam_fix(output_path+'.star.Aligned.sortedByCoord.out.bam', output_path+'.sort.bam')
    sp.run(['rm',output_path+'.star.Aligned.sortedByCoord.out.bam']) # Temporary files erasure
    
    bam_indexing(output_path + '.sort.bam', cmdlog)
    flagstat(output_path + '.sort.bam', cmdlog, threads)
    cmdlog.write("\n")
    cmdlog.close()
    
    if not rm:
        os.system("rm {path}*progress.out {path}*final.out {path}*SJ.out.tab".format(path=output_path))

def hisat2(reference, r1, r2, output, prefix, force, threads):
    """
    Run HiSat2.
    Call the aligner HiSat2
    :param ref: reference fasta file
    :param r1: reads files. If it's paired-end, precise the second reads files with r2
    :param r2: second read file for paired-end
    :param pref: output prefix
    :param threads: number of threads used to perform the alignment
    :return:
    """
    output_path = output + "/hisat2";
    if prefix:
        output_path += "_" + prefix;
    
    genomedir = output + "/hisat2_genomeRef" # Path to the directory will contain reference fasta file and index files 
    
    # Create output dir if not exist
    if not os.path.exists(genomedir) :
        os.makedirs(genomedir)
    if not force:
        try :
            if os.path.exists(output_path + ".sort.bam") or os.path.exists(output_path + ".sort.flagstat.txt"):
               raise FileExistsError
        except FileExistsError as e :
            print('\nError: output file(s) already exists.\n')
            exit(1)
    
    ref_path = "/".join([genomedir, os.path.basename(reference.name)]) ; # Path to reference fasta file and index files

    cmdlog = open(output_path+".log", "w")
    
    # Reference File indexing
    if not os.path.exists(ref_path + ".1.ht2") or force:
        index_command = ['hisat2-build', reference.name, ref_path]
        cmdlog.write('\n# Fasta indexing:\n')
        cmdlog.write(" ".join(index_command) + "\n")
        log = sp.check_output(index_command, stderr=sp.STDOUT)
        with open(output_path + '_build.log','w') as log_file :
            log_file.write(log.decode('utf-8'))

    # Reads Mapping and ouput files writing 
    hisat_command = ['hisat2',
                    '-x', ref_path,
                    '-p', str(threads)
                    ]
    
    if r2 :
        hisat_command.extend([
                    '-1', r1.name,
                    '-2',r2.name
                    ])
    else :
        hisat_command.extend([
                    '-U', r1.name
                    ])
    try :
        if '.fq' in r1.name or '.fastq' in r1.name :
            hisat_command.append('-q')
        elif '.fa' in r1.name or '.fasta' in r1.name :
            hisat_command.append('-f')
        else :
            raise OSError
    except OSError :
        print(r1.name+': File format not supported to call Hisat2. Only regular, gzipped or bzipped fastq and fasta files are supported.')
        exit(1)
    samtools_view_cmd = ['samtools', 'view', '-bS']
    samtools_sort_cmd = ['samtools','sort', '-o',output_path+'.sort.bam'] 
    cmdlog.write('\n# HiSat2 Alignement:\n')
    cmdlog.write(' | '.join([' '.join(hisat_command),' '.join(samtools_view_cmd), ' '.join(samtools_sort_cmd)]) + "\n")
    with open(output_path+'_aln.log','w') as log : 
        align = sp.Popen(hisat_command,stdout=sp.PIPE,stderr=log)
        view = sp.Popen(samtools_view_cmd, stdin = align.stdout, stdout=sp.PIPE)
        sort = sp.Popen(samtools_sort_cmd, stdin = view.stdout)
        align.stdout.close()
        view.stdout.close()
        sort.wait()
   
    bam_indexing(output_path+'.sort.bam', cmdlog)
    flagstat(output_path+'.sort.bam', cmdlog, threads)
    cmdlog.write("\n")
    cmdlog.close()
