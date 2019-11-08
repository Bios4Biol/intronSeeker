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

def bam_indexing(bamfile) :
    
    bam_index_cmd = ['samtools','index',
                    bamfile,
                    bamfile.rstrip('.bam') + '.bai'
                    ]
    print("\nBAM indexing :") 
    print(' '.join(bam_index_cmd))
    sp.run(bam_index_cmd)

def flagstat(bamfile,threads = 1) :

    bam_flagstat_cmd = ['samtools','flagstat',
                        '--threads',str(threads),
                        bamfile
                        ]
    print("\nBAM flagstat computation :") 
    print(' '.join(bam_flagstat_cmd))
    stdout = sp.run(bam_flagstat_cmd,stdout=sp.PIPE).stdout
    with open(bamfile.rstrip('.bam')+'.flagstat.txt','w') as out_flag :
        out_flag.write(stdout.decode('utf-8'))

def star(reference, r1, r2, prefix, threads):
    """
    Run STAR.
    Call the aligner STAR. 
    :param ref: reference fasta file
    :param r1: reads files. If it's paired-end, precise the second reads files with r2
    :param r2: second read file for paired-end
    :param pref: output prefix
    :return:
    """
    outdir = prefix + "_starAlignement/" ; # Output directory name (where this function will write 
                                 # all the results and tmp files).
    outfile = outdir + prefix ; # Final output file path
    genomedir = outdir + "GenomeRef" ; # Path to the directory will contain reference fasta file and index files 

    # Directories creation 
    try :
        os.makedirs(genomedir)
    except FileExistsError :
        print('WARNING: Output directories already exist. The potential existing output files will be overwritten.')
    sp.run(['cp', reference.name, genomedir])
    ref_path = "/".join([genomedir, os.path.basename(reference.name)]) ; # Path to reference fasta file and index files
    
    # Reference File indexing
    index_command = ['STAR',
                    '--runMode', 'genomeGenerate',
                    '--runThreadN', str(threads),
                    '--genomeSAindexNbases', str(6),
                    '--genomeDir', genomedir, 
                    '--genomeFastaFiles', ref_path
                    ]
    print('\nFasta indexing :')
    print(" ".join(index_command))
    sp.run(index_command) 
    sp.run(['mv', 'Log.out', outdir]) ; # Move the index log file in output directory

    # ~ # Reads Mapping and ouput files writing 
    star_command = ['STAR',
                    '--genomeDir', genomedir,
                    '--runThreadN', str(threads),
                    '--outSAMstrandField', 'intronMotif',
                    '--outSAMtype', 'BAM', 'SortedByCoordinate',
                    '--outSAMunmapped', 'Within',
                    '--outFileNamePrefix', outfile+'.star.'
                    ]
    if r2 :
        star_command.extend(['--readFilesIn', r1.name, r2.name])
    else :
        star_command.extend(['--readFilesIn', r1.name])
    if r1.name.endswith('.gz') :
        star_command.extend(['--readFilesCommand', 'zcat'])
    print('\nSTAR alignement : ')
    print(' '.join(star_command))
    sp.run(star_command)
    star_bam_fix(outfile+'.star.Aligned.sortedByCoord.out.bam', outfile+'.Aligned.sortedByCoord.out.bam')
    sp.run(['rm',outfile+'.star.Aligned.sortedByCoord.out.bam']) # Temporary files erasure
    
    bam_indexing(outfile + '.Aligned.sortedByCoord.out.bam')
    
    flagstat(outfile + '.Aligned.sortedByCoord.out.bam',threads)


def hisat2(reference, r1, r2, prefix, threads):
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
    outdir = prefix + "_hisat2Alignement/" # Output directory name (where this function will write)
                                  # all the results and tmp files).
    outfile = outdir + prefix # Final output file path
    genomedir = outdir + "GenomeRef" # Path to the directory will contain reference fasta file and index files 

    # Directories creation 
    try :
        os.makedirs(genomedir)
    except FileExistsError :
        print('WARNING: Output directories already exist. The potential existing output files will be overwritten.')
    sp.run(['cp', reference.name, genomedir])
    ref_path = "/".join([genomedir, os.path.basename(reference.name)]) ; # Path to reference fasta file and index files

    # Reference File indexing
    index_command = ['hisat2-build',reference.name ,ref_path]
    print('\nFasta indexing :')
    print(" ".join(index_command))
    log = sp.check_output(index_command) 
    with open(outdir+'hisat2-build.log','w') as log_file :
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
        print(r1.name+' : File format not supported to call Hisat2. Only regular, gzipped or bzipped fastq and fasta files are supported.')
        exit(1)
    samtools_view_cmd = ['samtools', 'view', '-bS']
    samtools_sort_cmd = ['samtools','sort', '-o',outfile+'.Aligned.sortedByCoord.out.bam'] 
    print('\nHiSat2 Alignement : ')
    print(' | '.join([' '.join(hisat_command),' '.join(samtools_view_cmd), ' '.join(samtools_sort_cmd)]))
    align = sp.Popen(hisat_command,stdout=sp.PIPE,stderr=sp.PIPE)

    view = sp.Popen(samtools_view_cmd, stdin = align.stdout, stdout=sp.PIPE)
    sort = sp.Popen(samtools_sort_cmd, stdin = view.stdout)
    
    bam_indexing(outfile+'.Aligned.sortedByCoord.out.bam')
    
    flagstat(outfile+'.Aligned.sortedByCoord.out.bam')
    # ~ if ".fq" in r1 or ".fastq" in r1 :
        # ~ os.system("hisat2 -x {ref_path} -1 {reads1} -2 {reads2} -q -p {threads} | samtools view -bS | samtools sort -o {outfilename}".format(
            # ~ ref_path=ref_path, reads1=r1, reads2=r2, threads=threads, outfilename=outfile + ".Aligned.sortedByCoord.out.bam")) ;
    # ~ else :
        # ~ os.system("hisat2 -x {ref_path} -1 {reads1} -2 {reads2} -f -p {threads} | samtools view -bS | samtools sort -o {outfilename}".format(
            # ~ ref_path=ref_path, reads1=r1, reads2=r2, threads=threads, outfilename=outfile + ".Aligned.sortedByCoord.out.bam")) ;
    
    # ~ os.system("samtools index {bamfile} {indexfile}".format(bamfile=outfile+".Aligned.sortedByCoord.out.bam",
                                                            # ~ indexfile=outfile+".Aligned.sortedByCoord.out.bai")) ;


