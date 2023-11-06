---
Title: Finding and removing introns from RNA-Seq de novo assemblies with IntronSeeker
Tags:
- RNA-Seq 
- de novo 
- assembly
- retained introns
Authors:
- Sarah Maman
- Philippe Bardou
- Emilien Lasguignes
- Faustine Oudin
- Floréal Cabanettes
- Christophe Klopp
Bibliography:
- IntronSeeker.bib
---

# Summary 

# Introduction

Short read RNA sequencing (RNA-Seq) is now routinely used to study gene expression. When a reference genome is available, RNA-Seq reads can be splice-aligned to the assembly and gene abundances can be measured by counting alignments found in each gene location. When no reference genome assembly is available, reads are usually assemble to build a reference transcriptome contig set. In this de novo approach reads are then aligned to the contigs without using a splice-aware aligner because they originate from mature transcripts.

In order not to sequence very abundant ribosomal RNAs, commonly used protocols include oligo(dT) transcript enrichment. Transcript poly-adhenylation and splicing take place in the nucleus before transfer in the cytoplasm. PolyA tail selection retrieves mainly mature spliced transcripts. Introns spanning reads are still found in most RNA-Seq sets as shown by  [@Ameur2011].

Intron retention is also a known biological gene regulation or alternative splicing mechanism. In plants it has been shown to increase the number of transcript splicing forms  [@Jia2020].  [@Braunschweig2014] have shown that, even in mammals, their is widespread intron retention linked with gene expression regulation.

Reads located on intron/exon boundaries harbor specific splice motifs. These splice motifs are di-nucleotides located at both intron ends. They can be canonical (GT/AG) or not. Canonical motifs are found in most site of most species.

Intron retention can be seen as an artefact or a biologically relevant mechanism and in both cases it is interesting to monitor. This is easy with a reference genome assembly because one can quantify reads aligned in introns. This is more complicated in the de novo approach in which the assembler can produce contigs with and without intron for the same transcript.

Contigs with introns are more difficult to annotated because introns are splitting coding sequences (CDS) and possibly generating several shorter protein alignments rather than a unique long match.

Canonical splice motifs presence, number of reads splicing at the same location, the minimum number of splice events and overlaps can be used to reduce faulty detection.

Hereafter we present IntronSeeker a software package enabling to search and remove introns from de novo assembled RNA-Seq contigs.

# Statement of need

RNA-Seq de novo assemblies are widely used to study transcription in species lacking a reference genome. The classical extraction protocol collects RNA fragments using their PolyA tails and there-fore should only gather mature RNAs. Unfortunately part of the extracted RNA molecules still com-prise retained introns which can therefore be present in some assembled contigs. These introns generate transcript frameshifts and thus render annotation more difficult. IntronSeeker searches introns by splice-realigning reads on contigs.  Introns candidates usually show canonical intron/exon boundaries supported by several sequences. IntronSeeker found between 0.02 and 11.96% of putative introns in twenty publicly available RNA-Seq de novo assembled samples. IntronSeeker produces an intron cleaned contig fasta file as well as a cleaning report. IntronSeeeker can be downloaded fromhttps://github.com/genotoul-bioinfo/intronseeker.

## Searching retained introns in public datasets 

The twenty public TSA contigs sets processed by IntronSeeker are classified in three Fungi, six Plantae, nine Animalia kingdoms and two Eukaryote superkingdoms. The number of contigs ranged from thirty thousand to three hundred thousand. The number of reads ranged from six millions to three hundred thirty millions. The proportion of retained intron candidates ranged from 0.02 to 11.96%. The figures are presented in Table 1.

| Species     |  (super)kingdom  | TSA id |Nb contigs | SRA and link to HTML intronSeeker report                                        |     Nb reads. |   \% cwi|
| :---       |  :----:   |   :----: |:----:  |    :----:  |      :--- |    :--- |
|  Salvia m.    |  Plantae   |  GJJN01 | 69 705 |  SRR15718805 |   23 086 599 | **11.96%**|
| Platichthys s.  |  Animalia  |  GAPK01 | 30 630 |  SRR1023744 | 516 791 904  | **10,71%**|
| Rigidoporus m.  |  Fungi    |  GDMN01 | 34 441 |  SRR2187438  |   75 600 628 |  **5,97%**|
| Isatis t.    |  Plantae   |  GARR01 | 33 238 | SRR1051997  |  113 134 348 |  **5,41%**|
| Goniomonas a.  |  Eukaryota  |  GGUN01 | 48 844 | SRR7601946   |  82 416 944 |  **5,27%**|
| Vriesea c.    |  Plantae   |  GHCB01 | 41 228 | SRR500874   |  85 726 288  | **3,94%**|
| Graminella n.  |  Animalia  |  GAQX01 | 37 537 | SRR857257   |  43 693 708  | **2,88%**|
| Carassius g.   |  Animalia  |  GJKR01 | 109 966 | SRR12596368  |  21 661 960 | **2.72%**|
| Rhizopus a.   |  Fungi    |  GDUK01 | 30 601 | SRR2104505   |  64 801 576  | **2,51%**|
| Diplonema p.   | Eukaryota  | GJNJ01  | 114 037 | SRR14933372  | 63 775 926   |  **2.35%**|
| Azolla f.    | Plantae   | GBTV01  | 36 091 | SRR1618559  |  122 059 452 |  **1,97%**|
| Tripidium r.   | Plantae   | GJDA01  | 106 494 | SRR14143372  | 15 357 748  | **1.50%**|
| Cimex l.     |  Animalia  | GBYH01  | 39 124 | SRR1660436   | 329 875 624  |  **1,24%**|
| Trichoplax sp. H2|  Animalia  | GFSF01   | 43 376 | SRR5819939   | 128 665 904  | **1,24%**|
| Phytophthora c. |  Plantae   | GBGX01   | 21 662 | SRR1206033   | 14 346 946  | **1.15%**|
| Piromyces sp.   |   Fungi   | GGXH01  | 124 096 | SRR7819335  |  15 615 535 |  **0.50%**|
| Sander l.    |  Animalia  | GJIW01  | 56 196 | SRR15372351  |  48 694 199  | **0.36%**|
| Rhodnius n.   |  Animalia  | GJJI01  | 67 217 | [SRR15602387](https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/blob/master/data/REAL_DATA/Rhodnius_n_SRR15602387/report_Rhod_SRR15602387.html)  | 6 775 534   | **0.11%**|
| Choromytilus c. | Animalia   | GJJD01  | 106 298 | [SRR15058678](https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/blob/master/data/REAL_DATA/Choromytilus_c_SRR15058678/report_Choromytilus_SRR15058678.html)  | 10 490 833   |  **0.10%**|
| Bombus t.    |  Animalia  | GHFS01  | 48 241 | SRR6148374   |  15 547 444  | **0.02%**|

Table: Rates of public contigs with introns, data from NCBI Transcript  Shotgun Archive (TSA) and Short Read Archive (SRA).

NCBI taxonomy browser, superkingdom used when kingdom not provided.

intronSeeker files are available in "master" branch : intronSeeker / data / REAL_DATA

Fraction of contigs with introns (cwi) = Nb of contigs PASS / Nb total of contigs


# Implementation

## General overview

IntronSeeker is a python script which includes three steps enabling to align read on contigs, find and remove retained introns and produce an html report. Figure 1. presents theses steps with input and output file formats.

IntronSeeker is open source (GNU General Public License) and can be download from <http://github.com/>

![IntronSeeker steps diagram](paper/Figures/IS_pipeline.png)


## Conda based installation procedure

To ease installation, intronSeeker includes an installation script (setup.sh) which run the installation of all the dependencies (1) but one (grinder) using conda and then installs grinder in the conda environment. Conda has to be installed beforehand. Installation can be checked using intronSeeker checkInstall program which will verify the presence and version of all dependencies. To run intronSeeker one has first to load the ISeeker_environment using conda 'source activate ISeeker_environment' command. An example data set named reduced_real_dataset is also provided in the data directory. It includes the result files which will enable manual comparison of the reference and produced intron files after complying the test. The installation procedure and the test command line can be found on the main page of the intronSeeker GITHUB WEB-page.

\(1\) IntronSeeker dependencies include seven software packages :
grinder [@angly2012grinder], gffread [@gffread], hisat2 [@kim2015hisat], STAR [@dobin2013star], samtools [@li2009sequence], TransDecoder [@haas2013novo], diamond [@buchfink2015fast].

Versions required before installation:
* Python version 3.6 or above.
* Miniconda 23.3.3 or above.

```
#Clone intronSeeker code from Git:
$ git clone https://forgemia.inra.fr/emilien.lasguignes/intronSeeker.git

#Load our miniconda environment and use libmamba
$ conda activate
$ conda update -n base conda
$ conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'

# Set up intronSeeker
$ cd intronSeeker/
$ CONDA_SOLVER="libmamba" /bin/bash ./setup.sh

#Activate ISeeker_environment and check installation
$ conda activate ISeeker_environment 
$ intronSeeker checkInstall 

# Command line help
$ intronSeeker -h
```

## Setting default detection parameters

Detected introns are filtered according to thresholds (number of reads overlapping an intron and maximum length), canonical junction (GT_AG or CT_AC), and complexity (too long or overlapping introns).

Different GBS simulations where performed to find the best possible thresholds for minimum number of splice events, maximum intron spanning size, presence of canonical splice sites and minimum overlap size to call a splice event an intron retention candidate.

![GBS parameters impacts](paper/Figures/GBS_params_impacts.png)

IntronSeeker parameters impacts. Graphics have been built on 10 samples of Arabidopsis thaliana data with powerlow abundance (1.2). Details are available in data/ directory.

**A** Coverage impact on detection: Increasing the coverage allows to quickly lose false candidates with a very limited impact on the number of true candidates (Arabidopsis thaliana data, powerlow abundance).

**B** Filter candidates on border improve detection. DPratio =DPin/(DPbefore+DPafter)

**C** Filter candidates on candidat length.

**D** Filter candidates on retained intron ratio (len contig / len candidat).

# Acknowledgement of financial support

The non permanent positions of Floréal Cabanette and Emilien Lasguignes were financed by projet France Génomique n° 31000523 and projet France Génomique n° 15000079.

