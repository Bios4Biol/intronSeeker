# Material and methods

Several analysis have been performed while developing IntronSeeker in order to select an aligner able to splice align reads even if both
contigs with and without intron are present in the reference, to set the correct parameters to collect only correct retained intron, reduce false positive rate, check impact of removing retained on protein completeness and finally to verify that public data sets present retained intron candidates. This section will first present how public data validation was performed before entering in the details of aligner selection and parameter setting.

## Finding retained introns in public datasets

NCBI TSA (Transcriptome Shotgun Assembly) contains over six thousand publicly available RNA-Seq de novo assembled contig sets for various species from different reigns. It also includes the link to the SRA (Short Read Archive) read set which were used for these de novo assemblies. We selected twenty different TSA contig sets and one of their SRA read set corresponding to the species presented in table 2 and ran IntronSeeker with default parameters on these sets. A command line example can be found in Supp. 1 (Supplementary material appendix 1). The corresponding result files can be downloaded from (data inrae : URL TO COMPLETE
http://genoweb.toulouse.inra.fr/ smaman/intronSeeker/DATA-INRAE/). The rate of contigs containing at least an intron candidate was extracted from the result file and collected in table 2.

## Measuring intron removal impact on protein completeness 

For four species *Caenorhabditis elegans*, *Drosophila melanogaster*, *Arabidopsis thaliana* and *Gallus gallus* we simulated thousand contigs and half a million read pairs from the Ensembl (version) genome FASTA and annotation GTF files with the IntronSeeker GTF Based Simulation for the contig FASTA generation and the reads simulation to generate FASTQ files(XXXsimulate function). This genome based simulation strategy is called GBS (Genome Based Simulation). The module outputs a bed file storing which intron of the transcript is retained enabling to check if IntronSeeker was able to detect it. For all transcripts in which the retained intron was overlaping the coding sequence (CDS), we compared protein completeness of the contigs before and after cleaning them with
IntronSeeker.

## Selecting a splice aligner to discover retained introns in de novo transcriptome contigs 

Two full random simulations (FRS) were performed. In this simulation mode contigs are simulated by random picking nucleotides in a uniform distribution. In the first simulation, called no-mix, ten sets containing each thousand contigs with introns were randomly generated with the IntronSeeker FullRandomSimulation function. This function outputs two contig fasta files, the first without intron and the second with introns, as well as a bed file with the intron locations in the with-intron file. The without-intron file is used to simulate reads in the following step. For both sets, contig metrics where calculated with the assemblathon_stats.py script. hundred fold read pair coverage was simulated for each set using the simulateReads function with default parameters. Reads pairs were then aligned to both contig sets with the IntronSeeker Alignment function using both STAR and HiSAT2 options. The resulting bam files were parsed to count reads having a splice pattern (block of N in the cigarline) corresponding to the intron location.

In the second simulation, called mix, the same protocol was applied but this time the contigs with and whitout introns, five hundred of each, were mixed in the reference set. Meaning that each contig was present in its with-intron and without-intron forms. This to check if the aligner is able to produce a splice alignment even if the reference without intron is available.

Results of both simulations were gathered in Figure 2.


|Cele  |  Dmel  |  Ggal |   Athal|
| :---        |    :----:   |     :----:   |          :--- |
|  Nb True with non canonical junctions  |  3     |  0    |   1  |     1|
|  Nb True with canonical junctions      |  6 018 |  5 242 |  5 545  | 2 711|
|  Nb False with non canonical junctions  | 250   |  36   |   419 |    18|
|  Nb False with canonical junctions     |  835    | 692    | 1 575 |  535|

Table: Filter on canonical junctions.


Maximum intron length.\
![GBS length impact](doc/GBS_len_impact.png)

## Defining an retained intron correction strategy from real data sets 

Public de novo transcriptome projects include often several read sets corresponding to different conditions or replicates. Searching retained introns with several independent sets will give partially overlapping lists of candidates which can be merged in different ways such as only considering candidates found in all lists or found at least in two lists. We call this a correction strategy.

Four public TSA contigs sets of *Salvia miltiorrhizafor*, *Anguilla anguilla*, *Gadus morhua* and *Drosophila miranda* having from three to seven SRA reads sets have been processed with IntronSeeker using default parameters and there candidats have been compared in Venn diagrams produced with jvenn  [@bardou2014].

# Results

## Searching retained introns in public datasets 

The twenty public TSA contigs sets processed by IntronSeeker are classified in three Fungi, six Plantae, nine Animalia kingdoms and two
Eukaryote superkingdoms. The number of contig ranged from thirty thousand to three hundred thousand. The number of reads ranged from six
millions to three hundreed thirty millions. The proportion of retained intron candidates ranged from 0.02 to 11.96%. The figures are presented in Table 1.

|  Species          |   (super)kingdom   |  TSA id Nb contigs  | SRA reports                                                                               |         Nb reads.  |     \% cwi|
| :---        |    :----:   |     :----:   |     :----:   |     :----:   |          :--- |
|   Salvia m.    |       Plantae   |         GJJN01 69 705  |     [SRR15718805](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/) |     23 086 599 |  **11.96%**|
|  Platichthys s.   |   Animalia    |       GAPK01 30 630   |    [SRR1023744](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)     | 516 791 904 |  **10,71%**|
|  Rigidoporus m.   |   Fungi     |         GDMN01 34 441    |   [SRR2187438](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)  |     75 600 628 |   **5,97%**|
|  Isatis t.      |     Plantae     |       GARR01 33 238   |    [SRR1051997](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)   |   113 134 348 |   **5,41%**|
|  Goniomonas a.  |     Eukaryota    |      GGUN01 48 844   |    [SRR7601946](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)   |    82 416 944 |   **5,27%**|
|  Vriesea c.  |        Plantae   |         GHCB01 41 228    |   [SRR500874](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)     |   85 726 288  |  **3,94%**|
|  Graminella n.   |    Animalia       |    GAQX01 37 537  |     [SRR857257](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)     |   43 693 708   | **2,88%**|
|  Carassius g.     |   Animalia     |      GJKR01 109 966 |     [SRR12596368](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)  |    21 661 960  |  **2.72%**|
|  Rhizopus a.  |       Fungi     |         GDUK01 30 601  |     [SRR2104505](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)    |   64 801 576   | **2,51%**|
|  Diplonema p.    |    Eukaryota     |     GJNJ01 114 037     | [SRR14933372](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)     | 63 775 926 |   **2.35%**|
|  Azolla f.   |        Plantae     |       GBTV01 36 091 |      [SRR1618559](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)   |   122 059 452 |   **1,97%**|
|  Tripidium r. |       Plantae    |        GJDA01 106 494   |   [SRR14143372](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)    |  15 357 748  |  **1.50%**|
|  Cimex l.         |   Animalia    |       GBYH01 39 124   |    [SRR1660436](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)    |  329 875 624 |   **1,24%**|
|  Trichoplax sp. H2 |  Animalia    |       GFSF01 43 376   |    [SRR5819939](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)     | 128 665 904   | **1,24%**|
|  Phytophthora c.  |   Plantae   |         GBGX01 21 662   |    [SRR1206033](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)  |     14 346 946  |  **1.15%**|
| Piromyces sp.  |     Fungi     |         GGXH01 124 096    |  [SRR7819335](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)   |    15 615 535 |   **0.50%**|
|  Sander l.       |    Animalia       |    GJIW01 56 196     |  [SRR15372351](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)   |   48 694 199  |  **0.36%**|
|  Rhodnius n.     |    Animalia        |   GJJI01 67 217 |      [SRR15602387](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)     |  6 775 534  |  **0.11%**|
|  Choromytilus c.   |  Animalia       |    GJJD01 106 298    |  [SRR15058678](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)     | 10 490 833 |   **0.10%**|
|  Bombus t.   |        Animalia  |         GHFS01 48 241   |    [SRR6148374](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)    |   15 547 444  |  **0.02%**|

Table: Rates of public contigs with introns, data from NCBI Transcript   Shotgun Archive (TSA) and Short Read Archive (SRA).

NCBI taxonomy browser, superkingdom used when kingdom not provided.
Clic on SRR ids to display intronSeeker reports.
fraction of contigs with introns (cwi) = Nb of contigs PASS / Nb total of contigs

## Selecting a splice aligner to discover retained introns in de novo transcriptome contigs 

![Mapper detection accuracy](doc/mapper_detection_accuracy.png)\

Splice aligners enable to find introns from RNA-Seq de novo contigs Mapping detection accuracy to find retained introns when only one
unspliced form exists in the contigs, on simulated FRS data. Both mappers in the no-mix state align well, no noticeable difference between the two.\
The aligners are designed to map to genomes but we notice that they also splice correctly to transcript contigs.

Not all splice aligners have the same results : selected the best for intron detection on simulated data. When there are sequences with and without introns in an assembly, are the mappers able to identify the retained introns ? Compared to secondary alignments, HiSat2 only generates secondary alignments in the mixed state, and STAR generates more secondary alignments than HiSat2. Compared to splits, quite logically, we observe more splits with the mix option since the FASTA file contains the transcript in both states (with and without retained intron). In the mixed state, STAR splits much less than HiSat2.
According to sensibility, specificity and F1 score, HiSat2 detect more introns than STAR. HiSat2 splits better, especially when the percentage of retained intron increases. If a read overlaps introns, then the aligner does not split. Thus, if the sequences are not split, the features will not be mapped and therefore not detected.

Taking into account these criteria (number of splits, sensibility, specificity, number of features detected), HiSat2 is better than STAR
for intron detection on simulated data. More introns detected with HiSat2, especially when the number of sequences with retained introns
increase.


## Defining an retained intron detection strategy from real data sets 


 | Species (Clades)      |      Nb contigs (TSA)  | Samples     |  Nb Seq    |   %Ir|
 | :---        |    :----:   |     :----:   |     :----:     |          :--- |
 | Salvia m.(Magnoliophyta)  |  69 705 (GJJN)  |    SRR13996120  | 22 503 120  | **12.52**%|
 |||                                                SRR15718796  | 20 169 850 |  **11.73**%|
 |||                                                 SRR15718799   18 291 366   **11.55**%
 |  Anguilla a. (Chordata)  |    60 179 (GFIC)    |  SRR1532756  |  19 688 393  | **2.91**%|
 |||                                                SRR1532757   | 35 234 589 |  **1.69**%
 |||                                                SRR1532758   | 19 061 539|   **2.04**%
 |||                                                SRR1532760   | 21 815 169|   **2.09**%
 |||                                                SRR1532761   | 17 971 382 |  **1.72**%
 | Gadus m.(Chordata)      |    50 607 (GFIX)   |   SRR2045415   | 22 436 517  | **3.39**%|
 |||                                                SRR2045416   | 36 457 070|   **1.51**%
 |||                                                SRR2045417   | 35 480 430 |  **1.80**%
 | Drosophila m.(Arthropoda) |  12 521 (GALP)  |    SRR364798   |  27 272 458 |  **0.08**%|
 |||                                                SRR364799    | 13 522 671 |  **0.11**%
 |||                                                SRR364800    | 25 615 616 |  **0.08**%
 |||                                                SRR364801    | 32 719 471|   **0.25**%
 |||                                                SRR364802    | 28 794 285 |  **0.17**%
 |||                                                SRR364803   |  13 682 311 |  **0.13**%
 |||                                                SRR364804   |  33 564 923  | **0.25**%

Table: Description of the different NCBI samples used.

Intron retention rate (%Ir) = Nb tot PASS candidates / Nb tot of contigs in the fasta.
Salvia miltiorrhiza BioProject PRJNA759933
Anguilla anguilla BioProject PRJNA256923
Gadus morhua BioProject PRJNA256972
Drosophila miranda BioProject PRJNA208862

Different tissues give different introns. To compare tissues, all samples are reduced to 20 million sequences of 75 nucleotides.

|  Anguilla a.                            |                   Salvia m.|
| :---        |            :--- |
|  ![Anguilla tissues comparison](doc/jVenn_Anguilla_tissues_comparison.png)  | ![Salvia tissues comparison](doc/jVenn_Salvia.png)|

Table: Intron retention comparison for different tissues from Anguilla a. (Animalae) and Salvia m. (Plantae).

Ovaries (SRR1532756), brain (SRR1532757) and muscles (SRR1532760) of Anguilla anguilla, according to the number of PASS candidates.
Replicate of the same tissu from Henan (SRR13996120), from Shaanxi (SRR15718796) and from Yunnan (SRR15718799) of Salvia miltiorrhiza,
according to the number of PASS candidates.


Species behave differently according to their genome complexity. 

The quality of the outputs depends on the complexity of the genome considered.

 | Species  | F1score mean  | std deviation|
 | :---        |      :----:     |          :--- |
 | Athal  |   99,12    |      0,48|
 | Cele   |   99,28    |      0,32|
 | Dmel   |   97,93    |      0,65|
 | Ggal   |   98,88    |      0,29|

Table: More or less complex genomes with repeats and noise (overlapping introns).


## Processing time

|Salvia miltiorrhiza  |  Anguilla anguilla   |   Gadus morhua     |      Drosophila miranda|
| :---        |      :----:      |      :----:     |          :--- |
|  Nb samples      |    3 fastq           |     5 fastq          |      3 fastq    |            7 fastq|
|  Nb tot sequences |   2 \* 60 964 336     |   2 \* 113 771 072    |   2 \* 94 374 017      |  2 \* 175 171 735|
|  Nb contigs      |    69 705 contigs    |     60 179 contigs     |    50 607 contigs      |   12 521 contigs|
|  Nb PASS candidats  | 26 357          |       6 299           |       3 396     |             137|
|  Elapsed time/mem   ||||                                                                      
|  HiSat2 alignment |   2:18:31 / 5.52 GB   |   5:57:51 / 6.86 GB     | 6:35:25 / 5.98 GB   |   4:01:27 / 1.89 GB|
|  Split Read Search |  1-10:28:52 / 2.22 GB |  1-16:07:04 / 3.22 GB  | 1-10:00:18 / 2.01 GB |  00:59:14 / 608.43 MB|
|  HTML report      |   00:56:14 / 16.72 GB  |  01:21:36 / 26.67 GB   | 01:11:45 / 29.23 GB  |  01:49:08 / 25.94 GB|

Table: IntronSeeker processing time and memory consumption for NCBI datasets.

Salvia miltiorrhiza BioProject PRJNA759933 TSA GJJN01
Anguilla anguilla BioProject PRJNA256923 TSA GFIC
Gadus morhua BioProject PRJNA256972 TSA GFIX
Drosophila miranda BioProject PRJNA208862 TSA GALP
Memory consumption is maximum RAM usage.


## Data availability

All data created or used during this study are openly available from DATA INRAE Archive at https://data.inrae.fr/ DOI\....

http://genoweb.toulouse.inra.fr/ smaman/intronSeeker/DATA-INRAE/\

|  Species    |     Samples' intronSeeker reports|
| :---        |         :--- |
|  Salvia m.      | [SRR13996120](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_salvia_SRR13996120.html)|
||                  [SRR15718796](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_salvia_SRR15718796.html)|
||                  [SRR15718799](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_salvia_SRR15718799.html)|
|  Anguilla a.  |   [SRR1532756](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Anguilla_CC_SRR1532756.html)|
||                  [SRR1532757](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Anguilla_CC_SRR1532757.html)|
||                  [SRR1532758](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Anguilla_CC_SRR1532758.html)|
||                  [SRR1532760](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Anguilla_CC_SRR1532760.html)|
||                  [SRR1532761](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Anguilla_CC_SRR1532761.html)|
|  Gadus m. |       [SRR2045415](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Gadus_SRR2045415.html)|
||                  [SRR2045416](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Gadus_SRR2045416.html)|
||                  [SRR2045417](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Gadus_SRR2045417.html)|
|  Drosophila m. |  [SRR364798](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/)|
||                  [SRR364799](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Droso_SRR364799.html)|
||                  [SRR364800](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Droso_SRR364800.html)|
||                  [SRR364801](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Droso_SRR364801.html)|
||                  [SRR364802](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Droso_SRR364802.html)|
||                  [SRR364803](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Droso_SRR364803.html)|
||                  [SRR364804](http://genoweb.toulouse.inra.fr/~smaman/intronSeeker/DATA-INRAE/REAL-reports/report_Droso_SRR364804.html)|

table: intronSeeker reports of NCBI samples used to optimize our strategy to find retained introns.

#  IntronSeeker SRA samples


|  Species (TSA)      |          Samples  |     Nb nucleotides/sequence |  Nb sequencies |  \% mapped   |
| :---        |      :----:      |      :----:     |     :----:     |          :--- |
|  Anguilla CC (GFIC)         |  SRR1532756  |  100 nt      |              19 688 393     | 82.11%    |  
||                               SRR1532757   | 100 nt              |      35 234 589    |  76.66 %     |
||                               SRR1532760   | 100 nt             |       21 815 169    |  83.55%  |    
|  Gadus morhua (GFIX)         | SRR2045415   | 100 nt              |      22 436 517    |  87.37%  |    
||                               SRR2045416  |  100 nt             |       36 457 070   |   75.63% |     
||                               SRR2045417  |  100 nt              |      35 480 430   |   87.28% |     
|  Drosophila miranda (GALP)   | SRR364798    | 101 nt               |     27 272 458    |  49.51%  |    
||                               SRR364800   | 101 nt               |     25 615 616  |    69.22% |     
||                               SRR364801   |  76 nt                |     32 719 471   |   41.49% |     
||                               SRR364802   |  76 nt                |     28 794 285   |   45.32% |     
|  Salvia miltiorrhiza (GJJN)  | SRR13996120  | 150 nt                |    22 503 120    |  79.09%  |    
||                               SRR15718796 |  150 nt               |     20 169 850   |   73.99% |     
||                               SRR15718799 |  150 nt               |     18 291 366  |    80.90%     | 

Table: SRA samples description.
%mapped= Nb reads mapped \* 100 / Nb tot reads

# More observations

## Different tissues give different introns

To compare tissues, all samples are reduced to 20 million sequences of 75 nt. In the main manuscript we should only have Venn diagrams without the bottom histogram or barplot and with the tissue names and not the dataset names. In the same plot we could have two or three examples.


| Samples             |      \% Retained introns |  
| :---                |            :---          |
|SRR1532756 - Ovaries |  **3.46%**               |       
|SRR1532757 - Brain   |  **3.36%**               |         
|SRR1532760 - Muscle  | **3.32%**                |         
| 3 tissues           | **0.90%**                |          

Table: introns found in different tissues of Anguilla anguilla data.


![The tissues are therefore specific. jVenn diagram shows PASS candidates common to all 3 tissues, to 2 tissues and those specific to
one tissue. ](doc/jVenn_Anguilla_tissues_comparison.png)

## Splice aligners enable to find introns from RNA-Seq de novo contigs

Mapping detection accuracy to find retained introns when only one unspliced form exists in the contigs, on simulated FRS data.

![STAR *versus* HISAT2 alignment results on no-mix state contig sets containing only contigs with an intron, and all sequencies contains a retained intron (r 100).](doc/mappers_accuracy_nb_detected_introns.png)

## Filter on mapping quality

The number of FPs drops from 138 (without filter) to 125 (MAPQ60) to 40 (MAPQ2). The more the quality is filtered, the more the number of FP drops.

Increase of the number of TP from 461 (without filter and MAPQ60) to 546 (MAPQ2).

![Without mapping filter.](doc/without_mapping_filter.png)

![MAPQ60 filter.](doc/MAPQ60.png)

![MAPQ2 filter.](doc/MAPQ2.png)

Test data: GBS Ggal, 1 000 contigs, mix state, tr 2000000, depth 1, HiSat2 mapping.

## Remove too complex introns because too long or overlapping

Detection is improved by removing too long introns.

![With and without too complex introns on GBS Ggal data, powerlaw 0.5, 1 000 contigs, 2 000 000 reads, mix states, depth 1 and d3 and HiSat2 mapping](doc/too_long_introns.png)

## The impact of the abundance bias on the detection for GBS data

The more complex the genome, as is the case for Gallus gallus, the more repeats, noise and overlapping introns there are. Thus, by adding the monotranscript parameter (-u), the results are optimized, especially for complex models, because the abundance bias is reduced. The -u option is not representative of reality, but it allows us to evaluate the impact of the abundance bias on the detection for GBS data, which are close to the real data. Let's study the impact of the suppression of monotranscript genes:

![Remove contigs with several transcripts improve detection accuracy.](doc/monotranscripts.png)

By keeping only the monotranscript genes, the sensitivity, and especially the specificity, is improved.

## Down sampling

![Mapping comparison between complete and reduced samples.](doc/DS_CompleteSample_pourcent_mapped_comparison.png)


##To complete

Randomly simulating data allows us to evaluate the impact of several dataset parameters on candidates detection. The aligner is important but also the data : the more uniform the read distribution the easier is the intron detection Impact of read abundance bias on the detection of retained introns. We therefore observe this impact between 0 and 2:

Powerlaws (pw) between 0.1 and 0.9 are similar to am (homogeneous). At more than 2, the abundance bias has no impact.

Cumulative coverage comparison between GBS simulation with powerlaw between 0.6 to 2, and real data (on 10 000 contigs): Specie : Atlantic cod - Gadus morhua - gadMor3.0.106 For real data : PRJNA256972 (50,607 contigs) - Sequence Read Archive: SRR2045415 - GFIX01

We checked these parameters on completely simulated and based on real information simulated data to define the best possible parameter set 

1/ Filter on mapping quality A non-modifiable filter built into the intronSeeker code is a filter on a minimum mapping quality of 2 for
HiSat2 (quality score: 0-1-60) and STAR (quality score: 1-2-3-255). Applying a filter on the mapping quality reduces the number of FPs and increases the number of TPs, thus improving detection: 

2/ The depth impacts the identification of the retained introns.

The aligner is wrong with a depth of 0 and optimal depth is between 1 and 3. intronSeeker parameters: HiSat, 1 000 contigs, 100% of sequences with retained intron, mix state (spliced state: find retained introns when both forms of the contig exist in the same fasta file), powerlaw 05.

4/ Filter on DP before/in/after candidate :

    DP_before = Mean DP for 10 bp before the candidate
    DP_in = Mean DP of candidate
    DP_after = Mean DP for 10 bp after the candidate

Samples Anguilla CC : SRR7820637 4h - SRR7820620 8h - SRR7820615 24h -
SRR7820616 24h ![jVenn Anguilla DP](doc/jVenn_Anguilla_DP.png)

 | Species        |       Samples     |  Nb candidats PASS  | Nb too complex introns|
 | :---        |      :----:      |   :----:      |          :--- |
 | Anguilla anguilla  |   SRR1532756   | 1 753      |         402 (68%+32%)|
 ||                       SRR1532757  |  1 015      |         564 (89%+11%)|
 ||                      SRR1532758  |  1 231      |         418 (88%+12%)|
 ||                       SRR1532760 |   1 261     |          480 (86%+14%)|
 ||                       SRR1532761 |   1 039     |          426 (86%+14%)|
  Gadus morhua       |   SRR2045415  |  1 720     |          477 (54%+46%)|
 ||                       SRR2045416 |   765     |            447 (89%+11%)|
 ||                       SRR2045417 |   911      |           696 (84%+16%)|
  Drosophila miranda |   SRR364798  |   10     |             4 319 (100%+0%)|
 ||                       SRR364799 |    14   |               2 646 (100%+0%)|
 ||                        SRR364800 |    10   |               2 375 (100%+0%)|
 ||                       SRR364801  |   32    |              5 282 ( 100%+ 0%)|
 ||                       SRR364802  |   22    |              395 (100%+0%)|
 ||                       SRR364803 |    17    |              2 441 (100%+0%)|
 ||                       SRR364804  |   32    |              435 (99%+1%)|
  Salvia miltiorrhiza |  SRR13996120  | 8 730  |             967 (14%+86%)|
 ||                       SRR15718796 |  8 181   |            847 (15%+85%)|
 ||                       SRR15718799  | 9 446   |            1005 (13%+87%)|

table: intronSeeker on NCBI datasets.


Nb too complex introns (%intron with len%+% Overlapping introns).
BioProject PRJNA256923 TSA GFIC 60 179 contigs
BioProject PRJNA256972 TSA GFIX 50 607 contigs
BioProject PRJNA208862 TSA GALP 12 521 contigs
BioProject PRJNA759933 TSA GJJN01 69 705 contigs


# Discussion 

-   complex contigs with overlapping splice alignments, complex contig rate depends on the genome complexity

-   The complexity of the reference genome impacts the detection of retained introns.

-   By comparing species on GBS data, we observe that the quality of the outputs depends on the complexity of the genome considered.
 
-   tissue related strategy : should we use all available tissues to check for introns or can we simplify the process

-   With different sequence sizes, different numbers of sequences, different depths per sample/tissue, it is difficult to compare     species and tissues within the same species.

-   Hence the need to reduce FASTQs with the same number of sequences and the same length of sequences.

-   To verify that reduce the dataset does not impact the rate of sequence alignment with the BAM flagstats, we compare the alignment between full samples and reduced samples (Down Sampling DS) at 20 million reads and 75 nucleotides.

-   what about long reads : reads to reads alignment, no more contigs

-   using IntronSeeker to find different contigs corresponding to the same transcript with a splice or intron retention event
