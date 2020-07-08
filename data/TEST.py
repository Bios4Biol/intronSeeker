#!/usr/bin/env python3

import os
import argparse 
from argparse import ArgumentParser
import configparser # To parse parameters file

# cd data
# python3 TEST.py -m FRS -g "/home/Sarah/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/FRS_CAS-B/config/grinder_frs_testB-3.cfg" -p "TOTO" --settingsFRS  " -n 1000 -r 50 --mix " -o "/home/smaman/Documents/" --settingsHISAT2 " " --settingsSTAR " " --settingsSRS " " --settingsTF " "
# python3 TEST.py -m FRS -g "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/FRS_CASC_sample6_pw_r20_mix/grinder_frs_testC-pw05.cfg" -p "FRS_CASC_sample6_pw_r20_mix" --settingsFRS  " -n 1000 -r 20 --mix " -o "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/FRS_CASC_sample6_pw_r20_mix" --settingsHISAT2 " " --settingsSTAR " " --settingsSRS " " --settingsTF " "
# python3 TEST.py -m GBS --gtf  "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/GBS/Cele/Caenorhabditis_elegans.WBcel235.97.gtf" --fasta "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/GBS/Cele/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa" -g "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/GBS/Cele/grinder_analysis_homogen.txt" -p "GBS_Cele" --settingsGBS  " -n 1000 " -o "/home/smaman/Documents/" --settingsHISAT2 " " --settingsSTAR " " --settingsSRS " " --settingsTF " "


############
# SUB MAIN #
############

def write_cgf_file(alignment:str, testsDir: str, pipelineName:str, mode:str ):
	print('Step 8/8 : Simulation Report with '+alignment)
	lowerAlignment = alignment.lower()
	if (mode == "FRS"):
		mfasta= testsDir+'/frs_'+pipelineName+'_contigs-modified.fa'
	if (mode == 'GBS'):
		mfasta = pipelineName+'/gbs_'+pipelineName+'_transcripts-modified.fa'		
	with open(testsDir+'/'+pipelineName+'_'+alignment+'.cfg', 'w+') as f:
			print('[Defaults]\nfasta: '+testsDir+'/frs_'+pipelineName+'_contigs.fa\nmfasta: '+mfasta+'\
				\ngtf: '+testsDir+'/frs_'+pipelineName+'_contigs-modified.gtf\nr1: '+testsDir+'/sr_'+pipelineName+'_R1.fastq.gz\
				\nr2: '+testsDir+'/sr_'+pipelineName+'_R2.fastq.gz\nflagstat: '+testsDir+'/'+lowerAlignment+'_'+pipelineName+'.sort.flagstat.txt\
				\ncandidat: '+testsDir+'/srs_'+pipelineName+'_'+alignment+'_candidates.txt\nranks: '+testsDir+'/sr_'+pipelineName+'_ranks.txt\
				\nsplit: '+testsDir+'/srs_'+pipelineName+'_'+alignment+'_split_alignments.txt\nprefix: '+pipelineName+'_'+alignment+'\nthreads: 6\
				\noutput: '+testsDir+'/HTML\nforce: -F', file=f)			
	return f


def mkdirWorkDir(workDir: str, pipelineName:str):
    os.chdir(workDir) 
    testsDir= workDir+"/"+pipelineName
    try :
        if os.path.exists(testsDir):
            raise FileExistsError
    except FileExistsError as e :
        print('\nError: output file already exists.\n')
        exit(1)

    # mkdir tests outputs directory
    testsDir=workDir+"/"+pipelineName
    os.system('mkdir '+ testsDir)
    os.system('chmod 777 '+ testsDir)
    print('Your work directory is: ', workDir+"/"+pipelineName)

    return testsDir
	



def simulationTests(mode: str, grinder: str, pipelineName: str, workDir:str, gtf:str, fasta:str, settingsFRS: str, settingsGBS:str, settingsSTAR: str, settingsHISAT2: str, settingsSRS: str, settingsTF: str):

	if (mode == "InstallIntronSeeker"):
		#Installation IntronSeeker
		os.system('./setup.sh')

	#Tester l'installation
	#Vérification de l'installation d'IntronSeeker et des ses dépendances (présence & version) avec le programme intronSeeker checkInstallEdit
	if (mode == "TestInstallIntronSeeker"):
		os.system('intronSeeker checkInstall')	

	#Full-random data simulation   ==> TODO : re-ecrire https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/blob/master/doc/SIMULATION.md
	if (mode == "FRS"):
		print('FRS mode')
		testsDir = mkdirWorkDir(workDir, pipelineName)
		print('Step 1/8 : Generate contigs fasta')
		os.system('intronSeeker fullRandomSimulation '+settingsFRS+' -o '+pipelineName+' -p '+pipelineName)
		print('Step 2/8 : Generate reads with intronSeeker simulateReads from reference fasta')
		os.system('intronSeeker simulateReads -o '+pipelineName+' -p '+pipelineName+' -f '+pipelineName+'/frs_'+pipelineName+'_contigs.fa -c '+grinder)
		print('Step 3/8 : starAlignment')
		os.system('intronSeeker starAlignment '+settingsSTAR+' -r '+pipelineName+'/frs_'+pipelineName+'_contigs-modified.fa -1 '+pipelineName+'/sr_'+pipelineName+'_R1.fastq.gz \
			    -2 '+pipelineName+'/sr_'+pipelineName+'_R2.fastq.gz -o '+pipelineName+' -p '+pipelineName)
		print('Step 4/8 : splitReadSearch post STAR')
		os.system('intronSeeker splitReadSearch '+settingsSRS+' -a '+pipelineName+'/star_'+pipelineName+'.sort.bam \
			    -r '+pipelineName+'/frs_'+pipelineName+'_contigs-modified.fa \
		        -o '+pipelineName+' -p '+pipelineName+'_STAR')
		print('Step 5/8 : nOK - trimFastaFromTXT post STAR')
		print('Step 6/8 : hisat2Alignment')
		os.system('intronSeeker hisat2Alignment '+settingsHISAT2+' -r '+pipelineName+'/frs_'+pipelineName+'_contigs-modified.fa \
		        -1 '+pipelineName+'/sr_'+pipelineName+'_R1.fastq.gz -2 '+pipelineName+'/sr_'+pipelineName+'_R2.fastq.gz -o '+pipelineName+' -p '+pipelineName)
		print('Step 7/8 : splitReadSearch post HISAT2')
		os.system('intronSeeker splitReadSearch '+settingsSRS+' -a '+pipelineName+'/hisat2_'+pipelineName+'.sort.bam \
		        -r '+pipelineName+'/frs_'+pipelineName+'_contigs-modified.fa -o '+pipelineName+' -p '+pipelineName+'_HISAT2')		
		write_cgf_file('STAR', testsDir, pipelineName, mode)
		write_cgf_file('HISAT2', testsDir, pipelineName, mode)
		os.system('python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file '+ testsDir+'/'+pipelineName+'_STAR.cfg')
		os.system('python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file '+ testsDir+'/'+pipelineName+'_HISAT2.cfg')


	#GBS
	if (mode == "GBS"):
		print('GBS mode')
		testsDir=mkdirWorkDir(workDir, pipelineName)
		print('Step 1/8 : Generate contigs fasta')
		os.system('intronSeeker.py GTFbasedSimulation '+settingsGBS+' -a '+gtf+' -r '+fasta+' -p '+pipelineName+' -o '+pipelineName)
		print('Step 2/8 : Generate reads with intronSeeker simulateReads from reference fasta')
		os.system('intronSeeker.py simulateReads -f '+pipelineName+'/gbs_'+pipelineName+'_transcripts.fa -c '+grinder+' -p '+pipelineName+' -o '+pipelineName)
		print('Step 3/8 : starAlignment')
		os.system('intronSeeker starAlignment -r '+pipelineName+'/gbs_'+pipelineName+'_transcripts-modified.fa -1 '+pipelineName+'/sr_'+pipelineName+'_R1.fastq.gz \
			       -2 '+pipelineName+'/sr_'+pipelineName+'_R2.fastq.gz -o  -o '+pipelineName+' -t 6   -p '+pipelineName)
		print('Step 4/8 : splitReadSearch post STAR')
		os.system('intronSeeker splitReadSearch -a '+pipelineName+'/star_'+pipelineName+'.sort.bam -r '+pipelineName+'/gbs_'+pipelineName+'_transcripts.fa -o starSplit  -o '+pipelineName+' -p '+pipelineName+' -t 6')
		print('Step 5/8 : nOK - trimFastaFromTXT post STAR')
		print('Step 6/8 : hisat2Alignment')
		os.system('intronSeeker hisat2Alignment '+settingsHISAT2+' -r '+pipelineName+'/gbs_'+pipelineName+'_transcripts-modified.fa \
		        -1 '+pipelineName+'/sr_'+pipelineName+'_R1.fastq.gz -2 '+pipelineName+'/sr_'+pipelineName+'_R2.fastq.gz -o '+pipelineName+' -p '+pipelineName)
		print('Step 7/8 : splitReadSearch post HISAT2')
		os.system('intronSeeker splitReadSearch '+settingsSRS+' -a '+pipelineName+'/hisat2_'+pipelineName+'.sort.bam \
		        -r '+pipelineName+'/gbs_'+pipelineName+'_transcripts-modified.fa -o '+pipelineName+' -p '+pipelineName+'_HISAT2')	
		write_cgf_file('STAR', testsDir, pipelineName, mode)
		write_cgf_file('HISAT2', testsDir, pipelineName, mode)
		os.system('python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file '+ testsDir+'/'+pipelineName+'_STAR.cfg')
		os.system('python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file '+ testsDir+'/'+pipelineName+'_HISAT2.cfg')


if __name__ == '__main__' :
    
    parser = ArgumentParser()
   	
    parser.add_argument('-m','--mode', type=str, required=True, dest='mode', help="Tests on FRS or GBS or REAL data.")
    parser.add_argument('-g','--grinder', type=str, required=True, dest='grinder', help="Path to the config grinder file (grinder*.cfg).")
    parser.add_argument('-p','--pipelineName', type=str, required=True, dest='pipelineName', help="Name of your test pipeline.")
    parser.add_argument('-o','--workDir', type=str, required=True, dest='workDir', help="Path to outputs dir")
    parser.add_argument('--fasta', type=str, required=False, dest='fasta', help="Reference fasta file.")
    parser.add_argument('--gtf', type=str, required=False, dest='gtf', help="GTF file.")
    parser.add_argument('--settingsFRS', type=str, required=False, dest='settingsFRS', help="FRS settings.")
    parser.add_argument('--settingsGBS', type=str, required=False, dest='settingsGBS', help="GBS settings.")
    parser.add_argument('--settingsSTAR', type=str, required=False, dest='settingsSTAR', help="STAR settings.")
    parser.add_argument('--settingsHISAT2', type=str, required=False, dest='settingsHISAT2', help="HISAT2 settings.")
    parser.add_argument('--settingsSRS', type=str, required=False, dest='settingsSRS', help="SRS settings.")
    parser.add_argument('--settingsTF', type=str, required=False, dest='settingsTF', help="TF settings.")

    args = parser.parse_args()
    simulationTests(**vars(args))

