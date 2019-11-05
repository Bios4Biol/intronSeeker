#!/usr/bin/env python3

try :
	import os ; 
	import pysam ;
except ImportError as error :
	print(error) ;
	exit(1) ;

def star_bam_fix(input_file: str, output_file):
	"""
	Rename each read from STAR alignment according to which mate it is (mate 1 or mate 2)
	:param input_file: name of the bam which contains STAR alignment 
	:param output_file: name of the ouput file 
	:return: nothing
	"""
	bamfile = pysam.AlignmentFile(input_file, "rb",check_sq=False)
	outfile = pysam.AlignmentFile(output_file, "wb", template=bamfile)
	for alignment in bamfile:
		if read.is_read1:
			alignment.query_name += "/1"
		if read.is_read2:
			alignment.query_name += "/2"
		outfile.write(read)


def star(reference: str, r1: str, r2: str, prefix: str, threads: int):
	"""
	Run STAR.
	Call the aligner STAR. 
	:param ref: reference fasta file
	:param r1: reads files. If it's paired-end, precise the second reads files with r2
	:param r2: second read file for paired-end
	:param pref: output prefix
	:return:
	"""
	outdir = prefix + "_star/" ; # Output directory name (where this function will write 
								 # all the results and tmp files).
	outfile = outdir + prefix ; # Final output file path
	genomedir = outdir + "GenomeRef" ; # Path to the directory will contain reference fasta file and index files 

	# Directories creation 
	os.system("mkdir " + outdir) ;
	os.system("mkdir " + genomedir) ;
	os.system("cp {reference} {directory}".format(reference = reference, directory = genomedir)) ;
	ref_path = "/".join([genomedir, reference.split("/")[-1]]) ; # Path to reference fasta file and index files

	# Reference File indexing
	os.system("STAR --runMode genomeGenerate --runThreadN {threads} --genomeSAindexNbases 6 --genomeDir {genomedir} "
			  "--genomeFastaFiles {ref_path}".format(threads=threads, genomedir=genomedir, ref_path=ref_path)) ; 
	os.system("mv Log.out ./" +outdir) ; # Move the index log file in output directory

	# Reads Mapping and ouput files writing 
	if r1.endswith(".gz") : # check if reads files are zipped
		os.system("STAR --genomeDir {genomedir} --runThreadN {threads} --readFilesIn {reads1} {reads2} "
				"--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --readFilesCommand zcat "
				"--outFileNamePrefix {outfile}".format(genomedir=genomedir, threads=threads, reads1=r1, reads2=r2, outfile = outdir + "temp")) ;
	else :
		os.system("STAR --genomeDir {genomedir} --runThreadN {threads} --readFilesIn {reads1} {reads2} "
				"--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within "
				"--outFileNamePrefix {outfile}".format(genomedir=genomedir, threads=threads, reads1=r1, reads2=r2, outfile = outdir + "temp")) ;
	star_bam_fix(outdir + "tempAligned.sortedByCoord.out.bam", outfile + ".Aligned.sortedByCoord.out.bam") ;
	os.system("rm " + outdir + "temp*") ; # Temporary files erasure
	print("BAM indexing...") 
	os.system("samtools index {inputfile} {outputfile}".format(inputfile=outfile + ".Aligned.sortedByCoord.out.bam",
															   outputfile=outfile + ".Aligned.sortedByCoord.out.bai")) ;


def hisat2(reference: str, r1: str, r2: str, prefix: str, threads: int):
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
	outdir = prefix + "_hisat2" ; # Output directory name (where this function will write)
								  # all the results and tmp files).
	outfile = outdir + "/" + prefix ; # Final output file path
	genomedir = outdir + "/GenomeRef" ; # Path to the directory will contain reference fasta file and index files 

	# Directories creation 
	os.system("mkdir " + outdir) ;
	os.system("mkdir " + genomedir) ;
	os.system("cp {reference} {genomedir}".format(reference=reference, genomedir=genomedir)) ;
	ref_path = "/".join([genomedir, reference.split("/")[-1] + ".index"]) ; # Path to reference fasta file and index files

	# Reference File indexing
	os.system("hisat2-build  {reference} {ref_path} > {directory}/hisat2-build.log".format(reference=reference, ref_path=ref_path, directory=outdir)) ;

	# Reads Mapping and ouput files writing 
	if ".fq" in r1 or ".fastq" in r1 :
		os.system("hisat2 -x {ref_path} -1 {reads1} -2 {reads2} -q -p {threads} | samtools view -bS | samtools sort -o {outfilename}".format(
			ref_path=ref_path, reads1=r1, reads2=r2, threads=threads, outfilename=outfile + ".Aligned.sortedByCoord.out.bam")) ;
	else :
		os.system("hisat2 -x {ref_path} -1 {reads1} -2 {reads2} -f -p {threads} | samtools view -bS | samtools sort -o {outfilename}".format(
			ref_path=ref_path, reads1=r1, reads2=r2, threads=threads, outfilename=outfile + ".Aligned.sortedByCoord.out.bam")) ;
	os.system("samtools index {bamfile} {indexfile}".format(bamfile=outfile+".Aligned.sortedByCoord.out.bam",
															indexfile=outfile+".Aligned.sortedByCoord.out.bai")) ;


