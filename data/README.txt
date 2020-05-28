h1. Tests sur données FRS

h2. CAS-A/ Statistiques d'alignement par défaut de STAR et HISAT2

h3. 1 - Step 1 : Full Random Simulation

source activate ISeeker_environment;

cd /yourWorkDir/

intronSeeker fullRandomSimulation -n 1000 -r  -o FRS_CAS-A/sample1/ -p sample1;
intronSeeker fullRandomSimulation -n 1000 -r  -o FRS_CAS-A/sample2/ -p sample2;
intronSeeker fullRandomSimulation -n 1000 -r  -o FRS_CAS-A/sample3/ -p sample3;
intronSeeker fullRandomSimulation -n 1000 -r  -o FRS_CAS-A/sample4/ -p sample4;
intronSeeker fullRandomSimulation -n 1000 -r  -o FRS_CAS-A/sample5/ -p sample5;

Results:
(ISeeker_environment) $ ls FRS_CAS-A/sample2/
frs_sample2_contigs.fa  frs_sample2_contigs-modified.fa  frs_sample2_contigs-modified.gtf

h3. 2 - Step 2 : Simulate Reads

cd /yourWorkDir/FRS_CAS-A/sample1/
intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testA.cfg


cd /yourWorkDir/FRS_CAS-A/sample1/;
intronSeeker simulateReads -o ./ -f frs_sample2_contigs.fa -c ../config/grinder_frs_testA.cfg
cd /yourWorkDir/FRS_CAS-A/sample2/;
intronSeeker simulateReads -o ./ -f frs_sample3_contigs.fa -c ../config/grinder_frs_testA.cfg
cd /yourWorkDir/FRS_CAS-A/sample3/;
intronSeeker simulateReads -o ./ -f frs_sample4_contigs.fa -c ../config/grinder_frs_testA.cfg
cd /yourWorkDir/FRS_CAS-A/sample4/;
intronSeeker simulateReads -o ./ -f frs_sample5_contigs.fa -c ../config/grinder_frs_testA.cfg


h3. 3 - Step 3 : STAR and HISAT2 alignment and Split Read Search

cd /yourWorkDir/FRS_CAS-A/sample1;
(ISeeker_environment) $  intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment;
(ISeeker_environment) $ ls 
frs_sample1_contigs.fa  frs_sample1_contigs-modified.fa  frs_sample1_contigs-modified.gtf  sr.log  sr_R1.fastq.gz  sr_R2.fastq.gz  sr_ranks.txt  STAR_alignment

(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR ;
(ISeeker_environment) $ ls
frs_sample1_contigs.fa           frs_sample1_contigs-modified.fa.fai  sample1_splicing_event_STAR  sr_R1.fastq.gz  sr_ranks.txt
frs_sample1_contigs-modified.fa  frs_sample1_contigs-modified.gtf     sr.log                       sr_R2.fastq.gz  STAR_alignment

(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment;
(ISeeker_environment) $ ls
frs_sample1_contigs.fa           frs_sample1_contigs-modified.fa.fai  HISAT2_alignment             sr.log          sr_R2.fastq.gz  STAR_alignment
frs_sample1_contigs-modified.fa  frs_sample1_contigs-modified.gtf     sample1_splicing_event_STAR  sr_R1.fastq.gz  sr_ranks.txt

(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2;
(ISeeker_environment) $ ls
frs_sample1_contigs.fa           frs_sample1_contigs-modified.fa.fai  HISAT2_alignment               sample1_splicing_event_STAR  sr_R1.fastq.gz  sr_ranks.txt
frs_sample1_contigs-modified.fa  frs_sample1_contigs-modified.gtf     sample1_splicing_event_HISAT2  sr.log                       sr_R2.fastq.gz  STAR_alignment


cd FRS_CAS-A/sample2;
intronSeeker starAlignment -r frs_sample2_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample2_contigs-modified.fa -o sample2_splicing_event_STAR
intronSeeker hisat2Alignment -r frs_sample2_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample2_contigs-modified.fa -o sample2_splicing_event_HISAT2


cd FRS_CAS-A/sample3
intronSeeker starAlignment -r frs_sample3_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample3_contigs-modified.fa -o sample3_splicing_event_STAR
intronSeeker hisat2Alignment -r frs_sample3_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample3_contigs-modified.fa -o sample3_splicing_event_HISAT2


cd FRS_CAS-A/sample4;
intronSeeker starAlignment -r frs_sample4_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample4_contigs-modified.fa -o sample4_splicing_event_STAR
intronSeeker hisat2Alignment -r frs_sample4_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample4_contigs-modified.fa -o sample4_splicing_event_HISAT2


cd FRS_CAS-A/sample5;
intronSeeker starAlignment -r frs_sample5_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample5_contigs-modified.fa -o sample5_splicing_event_STAR
intronSeeker hisat2Alignment -r frs_sample5_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample5_contigs-modified.fa -o sample5_splicing_event_HISAT2

h3. IntronSeeker simulation Report 

For each sample:

python3 simulation2HTML.py -m /path/to/FRS_CAS-A/sample1/frs_sample1_contigs-modified.fa -f /path/to/FRS_CAS-A/sample1/frs_sample1_contigs.fa -g /path/to/CAS-A/sample1/frs_sample1_contigs-modified.gtf -o /path/to/FRS_CAS-A/sample1/HTML -p FRS_CASA_sample1_n1000_r_STAR -F  -1 /path/to/FRS_CAS-A/sample1/sr_R1.fastq.gz -2 /path/to/FRS_CAS-A/sample1/sr_R2.fastq.gz --flagstat /path/to/FRS_CAS-A/sample1/STAR_alignment/star.sort.flagstat.txt -c /path/to/FRS_CAS-A/sample1/sample1_splicing_event_STAR/srs_candidates.txt -r /path/to/FRS_CAS-A/sample1/sr_ranks.txt -s  /path/to/FRS_CAS-A/sample1/sample1_splicing_event_STAR/srs_split_alignments.txt  --assemblathon /path/to/FRS_CAS-A/sample1/sample1_splicing_event_STAR/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6

h2. CAS-B/ Impact de la couverture

h3. Test several coverages

For sample 1 (only), without -r fullRandomSimulation option (1 contig = 1 intron) nor -n, run several tests with covage 3, 6, 9, 12, 15, 18, 21 :
Remarque : Without fullRandomSimulation -n option because we add -cf 3 in grinder config file.

$ more config/grinder_frs_testB-3.cfg
-cf 3                             
#-tr 2000000
-rd 101
-id 400 normal 100
-mo FR
-md poly4 3e-3 3.3e-8 
-mr 90 10
-am uniform
-fq 1
-ql 32 2

Warning :  « Coverage influence directly intron detection. In fact, the lower the coverage is, the harder the detection is. »

cd /yourWorkDir/
(ISeeker_environment) $ intronSeeker fullRandomSimulation  -o FRS_CAS-B/sample1/cf3/ -p sample1
cd /yourWorkDir/FRS_CAS-B/sample1/cf3/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testB-3.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2

cd /yourWorkDir/
(ISeeker_environment) $ intronSeeker fullRandomSimulation  -o FRS_CAS-B/sample1/cf6/ -p sample1
cd /yourWorkDir/FRS_CAS-B/sample1/cf6/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testB-6.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2

cd /yourWorkDir/
(ISeeker_environment) $ intronSeeker fullRandomSimulation  -o FRS_CAS-B/sample1/cf9/ -p sample1
cd /yourWorkDir/FRS_CAS-B/sample1/cf9/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testB-9.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2

cd /yourWorkDir/
(ISeeker_environment) $ intronSeeker fullRandomSimulation  -o FRS_CAS-B/sample1/cf12/ -p sample1
cd /yourWorkDir/FRS_CAS-B/sample1/cf12/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testB-12.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2

cd /yourWorkDir/
(ISeeker_environment) $ intronSeeker fullRandomSimulation  -o FRS_CAS-B/sample1/cf15/ -p sample1
cd /yourWorkDir/FRS_CAS-B/sample1/cf15/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testB-15.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2

cd /yourWorkDir/ 
(ISeeker_environment) $ intronSeeker fullRandomSimulation  -o FRS_CAS-B/sample1/cf18/ -p sample1
cd /yourWorkDir/FRS_CAS-B/sample1/cf18/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testB-18.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2


(ISeeker_environment) $ intronSeeker fullRandomSimulation  -o FRS_CAS-B/sample1/cf21/ -p sample1
cd /yourWorkDir/FRS_CAS-B/sample1/cf21/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testB-21.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2


h3. IntronSeeker simulation Report 

For each coverage:

python3 simulation2HTML.py -m /path/to/FRS_CAS-B/sample1/cf3/frs_sample1_contigs-modified.fa -f /path/to/FRS_CAS-B/sample1/cf3/frs_sample1_contigs.fa -g /path/to/FRS_CAS-B/sample1/cf3/frs_sample1_contigs-modified.gtf -o /path/to/FRS_CAS-B/sample1/cf3/HTML -p FRS_CASB_sample1_cf3 -F  -1 /path/to/FRS_CAS-B/sample1/cf3/sr_R1.fastq.gz -2 /path/to/FRS_CAS-B/sample1/cf3/sr_R2.fastq.gz --flagstat /path/to/FRS_CAS-B/sample1/cf3/STAR_alignment/star.sort.flagstat.txt -c /path/to/FRS_CAS-B/sample1/cf3/sample1_splicing_event_STAR/srs_candidates.txt -r /path/to/FRS_CAS-B/sample1/cf3/sr_ranks.txt -s  /path/to/FRS_CAS-B/sample1/cf3/sample1_splicing_event_STAR/srs_split_alignments.txt  --assemblathon /path/to/FRS_CAS-B/sample1/cf3/sample1_splicing_event_STAR/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6


h2. CAS-C/ Relative Abundance impact

Objective : Count number of input reference sequences by expression class with deux differents abundance models : uniform and powerlaw.

h3. Test with relative abundance model uniform

FRS_CAS-C/config/grinder_frs_testC-am-uniform.cfg
-tr 2000000
-rd 100 normal 7.5
-id 500 normal 250
-mo FR
-md poly4 3e-3 3.3e-8 
-mr 80 20 
-am uniform
-fq 1
-ql 32 2

cd /yourWorkDir/
(ISeeker_environment) $ intronSeeker fullRandomSimulation -n 1000 -r -o FRS_CAS-C/sample1/am/ -p sample1
cd /yourWorkDir/FRS_CAS-C/sample1/am/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testC-am-uniform.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2

h3. IntronSeeker simulation Report 


python3 simulation2HTML.py -m /path/to/FRS_CAS-C/sample1/am/frs_sample1_contigs-modified.fa -f /path/to/FRS_CAS-C/sample1/am/frs_sample1_contigs.fa -g /path/to/FRS_CAS-C/sample1/am/frs_sample1_contigs-modified.gtf -o /path/to/FRS_CAS-C/sample1/am/HTML -p FRS_CASC_sample1_am -F  -1 /path/to/FRS_CAS-C/sample1/am/sr_R1.fastq.gz -2 /path/to/FRS_CAS-C/sample1/am/sr_R2.fastq.gz --flagstat /path/to/FRS_CAS-C/sample1/am/STAR_alignment/star.sort.flagstat.txt -c /path/to/FRS_CAS-C/sample1/am/sample1_splicing_event_STAR/srs_candidates.txt -r /path/to/FRS_CAS-C/sample1/am/sr_ranks.txt -s  /path/to/FRS_CAS-C/sample1/am/sample1_splicing_event_STAR/srs_split_alignments.txt  --assemblathon /path/to/FRS_CAS-C/sample1/am/sample1_splicing_event_STAR/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6


h3. Test with relative abundance model powerlaw 0.5


FRS_CAS-C/config/grinder_frs_testC-pw05.cfg
-tr 2000000
-rd 100 normal 7.5
-id 500 normal 250
-mo FR
-md poly4 3e-3 3.3e-8 
-mr 80 20 
-am powerlaw 0.5
-fq 1
-ql 32 2

cd /yourWorkDir/
(ISeeker_environment) $ intronSeeker fullRandomSimulation -n 1000 -r -o FRS_CAS-C/sample1/pw05/ -p sample1
cd /yourWorkDir/FRS_CAS-C/sample1/pw05/
(ISeeker_environment) $ intronSeeker simulateReads -o ./ -f frs_sample1_contigs.fa -c ../../config/grinder_frs_testC-pw05.cfg
(ISeeker_environment) $ intronSeeker starAlignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o STAR_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a STAR_alignment/star.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_STAR
(ISeeker_environment) $ intronSeeker hisat2Alignment -r frs_sample1_contigs-modified.fa -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -o HISAT2_alignment
(ISeeker_environment) $ intronSeeker splitReadSearch -a HISAT2_alignment/hisat2.sort.bam -r frs_sample1_contigs-modified.fa -o sample1_splicing_event_HISAT2

h3. IntronSeeker simulation Report 

python3 simulation2HTML.py -m /path/to/FRS_CAS-C/sample1/pw05/frs_sample1_contigs-modified.fa -f /path/to/FRS_CAS-C/sample1/pw05/frs_sample1_contigs.fa -g /path/to/FRS_CAS-C/sample1/pw05/frs_sample1_contigs-modified.gtf -o /path/to/FRS_CAS-C/sample1/pw05/HTML -p FRS_CASC_sample1_pw05 -F  -1 /path/to/FRS_CAS-C/sample1/pw05/sr_R1.fastq.gz -2 /path/to/FRS_CAS-C/sample1/pw05/sr_R2.fastq.gz --flagstat /path/to/FRS_CAS-C/sample1/pw05/STAR_alignment/star.sort.flagstat.txt -c /path/to/FRS_CAS-C/sample1/pw05/sample1_splicing_event_STAR/srs_candidates.txt -r /path/to/FRS_CAS-C/sample1/pw05/sr_ranks.txt -s  /path/to/FRS_CAS-C/sample1/pw05/sample1_splicing_event_STAR/srs_split_alignments.txt  --assemblathon /path/to/FRS_CAS-C/sample1/pw05/sample1_splicing_event_STAR/srs_frs_sample1_contigs-modified_assemblathon.txt -t 6


