#!/usr/bin/env python3

import os
import argparse 
import subprocess
from argparse import ArgumentParser
import configparser # To parse parameters file

# cd data
# python3 TEST.py -m FRS -g "/home/Sarah/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/FRS_CAS-B/config/grinder_frs_testB-3.cfg" -p "TOTO" --settingsFRS  " -n 1000 -r 50 --mix " -o "/home/smaman/Documents/" --settingsHISAT2 " " --settingsSTAR " " --settingsSRS " " --settingsTF " "
# python3 TEST.py -m FRS -g "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/FRS_CASC_sample6_pw_r20_mix/grinder_frs_testC-pw05.cfg" -p "FRS_CASC_sample6_pw_r20_mix" --settingsFRS  " -n 1000 -r 20 --mix " -o "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/FRS_CASC_sample6_pw_r20_mix" --settingsHISAT2 " " --settingsSTAR " " --settingsSRS " " --settingsTF " "
# python3 TEST.py -m GBS --gtf  "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/GBS/Cele/Caenorhabditis_elegans.WBcel235.97.gtf" --fasta "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/GBS/Cele/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa" -g "/home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/data/GBS/Cele/grinder_analysis_homogen.txt" -p "GBS_Cele" --settingsGBS  " -n 1000 " -o "/home/smaman/Documents/" --settingsHISAT2 " " --settingsSTAR " " --settingsSRS " " --settingsTF " "


############
# SUB MAIN #
############

def write_cgf_file(alignment:str, testsDir: str, pipelineName:str, mode:str, fasta="", r1= "", r2= ""):
    print('Report Step : Simulation Report with '+alignment)
    lowerAlignment = alignment.lower()
    if (mode == 'FRS'):
        fa    = 'fasta: '+testsDir+'/frs_'+pipelineName+'_contigs.fa'
        mfasta= 'mfasta: '+testsDir+'/frs_'+pipelineName+'_contigs-modified.fa'
        gtf   = 'gtf: '+testsDir+'/frs_'+pipelineName+'_contigs-modified.gtf'
    if (mode == 'GBS'):
        fa     = 'fasta: '+testsDir+'/gbs_'+pipelineName+'_transcripts.fa'
        mfasta = 'mfasta: '+testsDir+'/gbs_'+pipelineName+'_transcripts-modified.fa'		
        gtf    = 'gtf: '+testsDir+'/gbs_'+pipelineName+'_transcripts-modified.gtf'
    if ((mode == 'GBS') or (mode== 'FRS')):
        fastq1 = 'r1: '+testsDir+'/sr_'+pipelineName+'_R1.fastq.gz'
        fastq2 = 'r2: '+testsDir+'/sr_'+pipelineName+'_R2.fastq.gz'
    if (mode == 'REAL'):
        fa     = fasta
        mfasta = ""
        fastq1 = r1
        fastq2 = r2    
    with open(testsDir+'/'+pipelineName+'_'+alignment+'.cfg', 'w+') as f:
                print('[Defaults]\n'+fa+'\n'+mfasta+'\
                \n'+gtf+'\n'+fastq1+'\
                \n'+fastq2+'\nflagstat: '+testsDir+'/'+lowerAlignment+'_'+pipelineName+'.sort.flagstat.txt\
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
    os.system('chmod -R 777 '+ testsDir)
    print('Your work directory is: ', workDir+"/"+pipelineName)

    return testsDir
	



def simulationTests(mode: str, grinder: str, pipelineName: str, workDir:str, gtf:str, fasta:str, r1: str, r2: str, settingsFRS: str, settingsGBS:str, settingsSTAR: str, settingsHISAT2: str, settingsSRS: str, settingsTF: str):

    if (mode == "InstallIntronSeeker"):
        #Installation IntronSeeker
        os.system('./setup.sh')

    #Check installation
    if (mode == "TestInstallIntronSeeker"):
        os.system('intronSeeker checkInstall')

    #Full-random data simulation   ==> TODO : re-ecrire https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/blob/master/doc/SIMULATION.md
    if (mode == "FRS"):
        print('FRS mode')
        testsDir = mkdirWorkDir(workDir, pipelineName)
        subprocess.check_call(["chmod 777 "+testsDir], shell=True)
        print('Step 1/8 : Generate contigs fasta')
        subprocess.check_call(["intronSeeker fullRandomSimulation "+settingsFRS+" -o "+pipelineName+" -p "+pipelineName],shell=True)
        print('Step 2/8 : Generate reads with intronSeeker simulateReads from reference fasta')
        subprocess.check_call(["intronSeeker simulateReads -o "+pipelineName+" -p "+pipelineName+" -f "+pipelineName+"/frs_"+pipelineName+"_contigs.fa -c "+grinder],shell=True)
        print('Step 3/8 : starAlignment')
        subprocess.check_call(["intronSeeker starAlignment "+settingsSTAR+" -r "+pipelineName+"/frs_"+pipelineName+"_contigs-modified.fa \
            -1 "+pipelineName+"/sr_"+pipelineName+"_R1.fastq.gz -2 "+pipelineName+"/sr_"+pipelineName+"_R2.fastq.gz -o "+pipelineName+" -p "+pipelineName],shell=True)
        print('Step 4/8 : splitReadSearch post STAR')
        subprocess.check_call(["intronSeeker splitReadSearch "+settingsSRS+" -a "+pipelineName+"/star_"+pipelineName+".sort.bam \
            -r "+pipelineName+"/frs_"+pipelineName+"_contigs-modified.fa \
            -o "+pipelineName+" -p "+pipelineName+"_STAR"],shell=True)
        print('Step 5/8 : trimFastaFromTXT post STAR')
        subprocess.check_call(["intronSeeker.py trimFastaFromTXT -r "+pipelineName+"/frs_"+pipelineName+"_contigs-modified.fa \
            -c "+testsDir+"/srs_"+pipelineName+"_STAR_candidates.txt -o "+pipelineName+" -p "+pipelineName],shell=True)
        print('Step 6/8 : hisat2Alignment')
        subprocess.check_call(["intronSeeker hisat2Alignment "+settingsHISAT2+" -r "+pipelineName+"/frs_"+pipelineName+"_contigs-modified.fa \
            -1 "+pipelineName+"/sr_"+pipelineName+"_R1.fastq.gz -2 "+pipelineName+"/sr_"+pipelineName+"_R2.fastq.gz -o "+pipelineName+" -p "+pipelineName],shell=True)
        print('Step 7/8 : splitReadSearch post HISAT2')
        subprocess.check_call(["intronSeeker splitReadSearch "+settingsSRS+" -a "+pipelineName+"/hisat2_"+pipelineName+".sort.bam \
            -r "+pipelineName+"/frs_"+pipelineName+"_contigs-modified.fa -o "+pipelineName+" -p "+pipelineName+"_HISAT2"],shell=True)
        print('Step 8/8 : trimFastaFromTXT post STAR')
        subprocess.check_call(["intronSeeker.py trimFastaFromTXT -r "+pipelineName+"/frs_"+pipelineName+"_contigs-modified.fa \
            -c "+testsDir+"/srs_"+pipelineName+"_HISAT2_candidates.txt -o "+pipelineName+" -p "+pipelineName],shell=True)		
        write_cgf_file('STAR', testsDir, pipelineName, mode)
        write_cgf_file('HISAT2', testsDir, pipelineName, mode)
        subprocess.check_call(["python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file "+ testsDir+"/"+pipelineName+"_STAR.cfg"],shell=True)
        subprocess.check_call(["python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file "+ testsDir+"/"+pipelineName+"_HISAT2.cfg"],shell=True)


    #GBS
    if (mode == "GBS"):
        print('GBS mode')
        testsDir=mkdirWorkDir(workDir, pipelineName)
        subprocess.run(["chmod -R 777 "+workDir], shell=True)
        subprocess.run(["chmod -R 777 "+testsDir], shell=True)
        print('Step 1/8 : Generate contigs fasta')
        print("intronSeeker.py GTFbasedSimulation "+settingsGBS+" -a "+gtf+" -r "+fasta+" -p "+pipelineName+" -o "+testsDir)
        subprocess.run(["intronSeeker.py GTFbasedSimulation "+settingsGBS+" -a "+gtf+" -r "+fasta+" -p "+pipelineName+" -o "+testsDir],shell=True)
        print('Step 2/8 : Generate reads with intronSeeker simulateReads from reference fasta')
        print("intronSeeker.py simulateReads -f "+testsDir+"/gbs_"+pipelineName+"_transcripts.fa -c "+grinder+" -p "+pipelineName+" -o "+pipelineName)
        subprocess.run(["intronSeeker.py simulateReads -f "+testsDir+"/gbs_"+pipelineName+"_transcripts.fa -c "+grinder+" -p "+pipelineName+" -o "+pipelineName],shell=True)
        print('Step 3/8 : starAlignment')
        print("intronSeeker starAlignment -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa -1 "+testsDir+"/sr_"+pipelineName+"_R1.fastq.gz  -2 "+testsDir+"/sr_"+pipelineName+"_R2.fastq.gz -o "+testsDir+" -t 6   -p "+pipelineName)
        subprocess.run(["intronSeeker starAlignment -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa -1 "+testsDir+"/sr_"+pipelineName+"_R1.fastq.gz \
                -2 "+testsDir+"/sr_"+pipelineName+"_R2.fastq.gz -o "+testsDir+" -t 6   -p "+pipelineName],shell=True,  timeout=None)
        print('Step 4/8 : splitReadSearch post STAR')
        print("intronSeeker splitReadSearch -a "+testsDir+"/star_"+pipelineName+".sort.bam -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa -o "+testsDir+" -p "+pipelineName+"_STAR -t 6")
        subprocess.run(["intronSeeker splitReadSearch -a "+testsDir+"/star_"+pipelineName+".sort.bam -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa -o "+testsDir+" -p "+pipelineName+"_STAR -t 6"],shell=True)
        print('Step 5/8 : trimFastaFromTXT post STAR')
        subprocess.run(["intronSeeker.py trimFastaFromTXT -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa \
                -c "+testsDir+"/srs_"+pipelineName+"_STAR_candidates.txt -o "+testsDir+"/STAR_trim/ -p "+pipelineName],shell=True)
        print('Step 6/8 : hisat2Alignment')
        subprocess.run(["intronSeeker hisat2Alignment "+settingsHISAT2+" -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa \
                -1 "+testsDir+"/sr_"+pipelineName+"_R1.fastq.gz -2 "+testsDir+"/sr_"+pipelineName+"_R2.fastq.gz -o "+testsDir+" -p "+pipelineName],shell=True)
        print('Step 7/8 : splitReadSearch post HISAT2')
        subprocess.run(["intronSeeker splitReadSearch "+settingsSRS+" -a "+testsDir+"/hisat2_"+pipelineName+".sort.bam \
                -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa -o "+testsDir+" -p "+pipelineName+"_HISAT2"],shell=True)
        print('Step 8/8 : trimFastaFromTXT post HISAT2')
        subprocess.run(["intronSeeker.py trimFastaFromTXT -r "+testsDir+"/gbs_"+pipelineName+"_transcripts-modified.fa \
                -c "+testsDir+"/srs_"+pipelineName+"_HISAT2_candidates.txt -o "+testsDir+"/HISAT2_trim/ -p "+pipelineName],shell=True)		
        write_cgf_file('STAR', testsDir, pipelineName, mode)
        write_cgf_file('HISAT2', testsDir, pipelineName, mode)
        subprocess.run(["python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file "+ testsDir+"/"+pipelineName+"_STAR.cfg"],shell=True)
        subprocess.run(["python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file "+ testsDir+"/"+pipelineName+"_HISAT2.cfg"],shell=True)

    #REAL
    if (mode == "REAL"):
        print('REAL mode')
        testsDir=mkdirWorkDir(workDir, pipelineName)
        subprocess.check_call(["chmod -R 777 "+testsDir], shell=True)
        print('Step 1/6 : starAlignment')
        subprocess.check_call(["intronSeeker starAlignment -r "+fasta+" -1 "+r1+" \
                -2 "+r2+" -o "+testsDir+" -t 6   -p "+pipelineName],shell=True)
        print("intronSeeker starAlignment -r "+fasta+" -1 "+r1+" -2 "+r2+" -o "+testsDir+" -t 6 -p "+pipelineName)        
        print('Step 2/6 : splitReadSearch post STAR')
        subprocess.check_call(["intronSeeker splitReadSearch -a "+pipelineName+"/star_"+pipelineName+".sort.bam -r "+fasta+" \
                 -o "+testsDir+" -p "+pipelineName+" -t 6"], shell=True)
        print("intronSeeker splitReadSearch -a "+pipelineName+"/star_"+pipelineName+".sort.bam -r "+fasta+" -o "+testsDir+" -p "+pipelineName+" -t 6")         
        print('Step 3/6 : trimFastaFromTXT post STAR')
        subprocess.check_call(["intronSeeker.py trimFastaFromTXT -r "+fasta+" \
                -c "+testsDir+"/"+pipelineName+"_STAR_candidates.txt -o "+testsDir+"/STAR_trim -p "+pipelineName], shell=True)	
        print("intronSeeker.py trimFastaFromTXT -r "+fasta+" -c "+testsDir+"/"+pipelineName+"_STAR_candidates.txt -o "+testsDir+"/STAR_trim -p "+pipelineName)
        print('Step 4/6 : hisat2Alignment')
        subprocess.check_call(["intronSeeker hisat2Alignment "+settingsHISAT2+" -r "+fasta+" \
                -1 "+r1+"   -2 "+r2+" -o "+testsDir+" -p "+pipelineName], shell=True)
        print("intronSeeker hisat2Alignment "+settingsHISAT2+" -r "+fasta+" -1 "+r1+"   -2 "+r2+" -o "+testsDir+" -p "+pipelineName)
        print('Step 5/6 : splitReadSearch post HISAT2')
        subprocess.check_call(["intronSeeker splitReadSearch "+settingsSRS+" -a "+pipelineName+"/hisat2_"+pipelineName+".sort.bam \
                -r "+fasta+" -o "+testsDir+" -p "+pipelineName+"_HISAT2"], shell=True)	
        print("intronSeeker splitReadSearch "+settingsSRS+" -a "+pipelineName+"/hisat2_"+pipelineName+".sort.bam -r "+fasta+" -o "+testsDir+" -p "+pipelineName+"_HISAT2")
        print('Step 6/6 : trimFastaFromTXT post HISAT2')
        subprocess.check_call(["intronSeeker.py trimFastaFromTXT -r "+fasta+" \
                -c "+testsDir+"/"+pipelineName+"_HISAT2_candidates.txt "+testsDir+"/HISAT2_trim -p "+pipelineName], shell=True)
        write_cgf_file('STAR', testsDir, pipelineName, mode, fasta, r1, r2)
        write_cgf_file('HISAT2', testsDir, pipelineName, mode)
        subprocess.check_call(["python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file "+ testsDir+"/"+pipelineName+"_STAR.cfg"], shell= True)
        subprocess.check_call(["python3 /home/smaman/Documents/PROJETS/INTRONSEEKER/intronSeeker/scripts/simulation2HTML.py -F --config_file "+ testsDir+"/"+pipelineName+"_HISAT2.cfg"], shell = True)




if __name__ == '__main__' :
    
    parser = ArgumentParser()
   	
    parser.add_argument('-m','--mode', type=str, required=True, dest='mode', help="Tests on FRS or GBS or REAL data.")
    parser.add_argument('-g','--grinder', type=str, required=False, dest='grinder', help="Path to the config grinder file (grinder*.cfg).")
    parser.add_argument('-p','--pipelineName', type=str, required=True, dest='pipelineName', help="Name of your test pipeline.")
    parser.add_argument('-o','--workDir', type=str, required=True, dest='workDir', help="Path to outputs dir")
    parser.add_argument('--fasta', type=str, required=False, dest='fasta', help="Reference fasta file.")
    parser.add_argument('--gtf', type=str, required=False, dest='gtf', help="GTF file.")
    parser.add_argument('--R1', type=str, required=False, dest='r1', help="Reads 1 fastq file.")
    parser.add_argument('--R2', type=str, required=False, dest='r2', help="Reads 2 fastq file.")
    parser.add_argument('--settingsFRS', type=str, required=False, dest='settingsFRS', help="FRS settings.")
    parser.add_argument('--settingsGBS', type=str, required=False, dest='settingsGBS', help="GBS settings.")
    parser.add_argument('--settingsSTAR', type=str, required=False, dest='settingsSTAR', help="STAR settings.")
    parser.add_argument('--settingsHISAT2', type=str, required=False, dest='settingsHISAT2', help="HISAT2 settings.")
    parser.add_argument('--settingsSRS', type=str, required=False, dest='settingsSRS', help="SRS settings.")
    parser.add_argument('--settingsTF', type=str, required=False, dest='settingsTF', help="TF settings.")

    args = parser.parse_args()
    simulationTests(**vars(args))

