---
title: 'Finding and removing introns from RNA-Seq de novo assemblies with IntronSeeker'
tags:
  - RNA-Seq 
  - de novo 
  - assembly
  - retained introns
authors:
  - name: Sarah Maman
    corresponding: true
    orcid: 0000-0001-6453-902X
    affiliation: 1
    email: sarah.maman@inrae.fr
  - name: Philippe Bardou
    orcid: 0000-0002-0004-0251 
    affiliation: 1
  - name: Emilien Lasguignes
    affiliation: 1
  - name: Faustine Oudin
    affiliation: 1
  - name: Floréal Cabanettes
    affiliation: 1
  - name: Christophe Klopp
    orcid: 0000-0001-7126-5477
    affiliation: 2
affiliations:
  - name: SIGENAE, GenPhySE, Université de Toulouse, INRAE, ENVT, Castanet Tolosan, France
    index: 1
  - name: SIGENAE, Bioinfo Genotoul, UMIAT, INRAE, Castanet-Tolosan, France
    index: 2
date: 2023
bibliography: paper.bib
---

# Summary 

intronSeeker identifies potentially retained introns in *de novo* RNA-seq assembly in order to quantify and remove them.
To use it you have to provide a set of contigs resulting of a *de novo* transcriptome assembly and a set of RNA-Seq reads. The reads will be aligned on the contigs, splices sites will be searched and tested to check if they correspond to valid intron retention events.
It includes two types of RNA-seq data simulation strategies to validate the detection process and measure the false positive detection rate.
intronSeeker provides a list of potential intron candidates, it filters them and outputs a clean reference without introns in Fasta format.

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
| Salvia m.    |  Plantae   |  GJJN01 | 69 705 |  [SRR15718805](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Salvia_m_GJJN01.html) |   23 086 599 | **11.96%**|
| Platichthys s.  |  Animalia  |  GAPK01 | 30 630 |  [SRR1023744](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Platichthys_s_GAPK01.html) | 516 791 904  | **10,71%**|
| Rigidoporus m.  |  Fungi    |  GDMN01 | 34 441 |  [SRR2187438](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Rigidoporus_m_GDMN01.html)  |   75 600 628 |  **5,97%**|
| Isatis t.    |  Plantae   |  GARR01 | 33 238 | [SRR1051997](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Isatis_t_GARR01.html)  |  113 134 348 |  **5,41%**|
| Goniomonas a.  |  Eukaryota  |  GGUN01 | 48 844 | [SRR7601946](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Goniomonas_a_GGUN01.html)   |  82 416 944 |  **5,27%**|
| Vriesea c.    |  Plantae   |  GHCB01 | 41 228 | [SRR500874](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Vriesea_c_GHCB01.html)   |  85 726 288  | **3,94%**|
| Graminella n.  |  Animalia  |  GAQX01 | 37 537 | [SRR857257](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Graminella_n_GAQX01.html)   |  43 693 708  | **2,88%**|
| Carassius g.   |  Animalia  |  GJKR01 | 109 966 | [SRR12596368](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Carassius_g_GJKR01.html)  |  21 661 960 | **2.72%**|
| Rhizopus a.   |  Fungi    |  GDUK01 | 30 601 | [SRR2104505](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Rhizopus_a_GDUK01.html)   |  64 801 576  | **2,51%**|
| Diplonema p.   | Eukaryota  | GJNJ01  | 114 037 | [SRR14933372](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Diplonema_p_GJNJ01.html)  | 63 775 926   |  **2.35%**|
| Azolla f.    | Plantae   | GBTV01  | 36 091 | [SRR1618559]( https://bios4biol.pages.mia.inra.fr/intronseeker/report_Azolla_f_GBTV01.html)  |  122 059 452 |  **1,97%**|
| Tripidium r.   | Plantae   | GJDA01  | 106 494 | [SRR14143372]( https://bios4biol.pages.mia.inra.fr/intronseeker/report_Tripidium_r_GJDA01.html)  | 15 357 748  | **1.50%**|
| Cimex l.     |  Animalia  | GBYH01  | 39 124 | [SRR1660436](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Cimex_l_GBYH01.html)   | 329 875 624  |  **1,24%**|
| Trichoplax sp. H2|  Animalia  | GFSF01   | 43 376 | [SRR5819939](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Trichoplax_sp_H2_GFSF01.html)   | 128 665 904  | **1,24%**|
| Phytophthora c. |  Plantae   | GBGX01   | 21 662 | [SRR1206033](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Phytophtora_c_GBGX01.html)   | 14 346 946  | **1.15%**|
| Piromyces sp.   |   Fungi   | GGXH01  | 124 096 | [SRR7819335](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Piromyces_sp_GGXH01.html)  |  15 615 535 |  **0.50%**|
| Sander l.    |  Animalia  | GJIW01  | 56 196 | [SRR15372351](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Sander_l_GJIW01.html)  |  48 694 199  | **0.36%**|
| Rhodnius n.   |  Animalia  | GJJI01  | 67 217 | [SRR15602387](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Rhodnius_n__SRR15602387.html)  | 6 775 534   | **0.11%**|
| Choromytilus c. | Animalia   | GJJD01  | 106 298 | [SRR15058678](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Choromytilus_c_GJJD01.html)  | 10 490 833   |  **0.10%**|
| Bombus t.    |  Animalia  | GHFS01  | 48 241 | [SRR6148374](https://bios4biol.pages.mia.inra.fr/intronseeker/report_Bombus_t_GHFS01.html)   |  15 547 444  | **0.02%**|

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
- Python version 3.6 or above.
- Miniconda 23.3.3 or above.

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

**A** Impact of coverage on detection.
Increased coverage means that false candidates are quickly lost, with very limited impact on the number of true candidates (tested with Arabidopsis thaliana data and powerlow abundance).

**B** Filter candidates on border improves detection.
DPratio =DPin/(DPbefore+DPafter)

Filtering candidates on the border enhances detection. 
DPratio is calculated as follows :
$$\frac{DP_{in}}{(DP_{before} + DP_{after})}$$
With DPbefore corresponding to the mean DP for 10bp before the candidate, DPin to the mean DP of candidate and DPafter to the mean DP for 10bp after the candidate.

**C** Filter candidates on candidat length.

**D** Filter candidates on retained intron ratio (candidat length / contig length).

# Acknowledgement of financial support

The non permanent positions of Floréal Cabanette and Emilien Lasguignes were financed by projet France Génomique n° 31000523 and projet France Génomique n° 15000079.

# References

