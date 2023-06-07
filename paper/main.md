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

# Introduction

Short read RNA sequencing (RNA-Seq) is now routinely used to study gene expression. When reference genome and transcriptome are available, reads can be splice-aligned to the genome and gene abundances can be measured by counting reads matching each gene location or they can be quantified without alignment using the presence of gene specific kmers. When no reference genome or transcriptome is available, reads are usually assemble to build a reference contig set. In this de novo approach reads are then aligned to the reference contigs without using a splice-aware aligner because they originate from mature transcripts.

In order not to sequence very abundant ribosomal RNAs, commonly used protocols include oligo(dT) transcript enrichment. Transcript
poly-adhenylation and splicing take place in the nucleus before transfer in the cytoplasm. PolyA tail selection retrieves mainly mature spliced transcripts. Introns spanning reads are still found in most RNA-Seq samples as shown by  [@Ameur2011].

Intron retention is also a known biological gene regulation or alternative splicing mechanism. In plants it has been shown to increase
the number of transcript splicing forms  [@Jia2020].  [@Braunschweig2014] have shown that even in mammals their is widespread intron retention and its link to gene expression regulation.

Reads located on intron/exon boundaries harbor specific splice motifs. These splice motifs are di-nucleotides located at both intron ends. They can be canonical (GT/AG) or not. Canonical motifs are found in most site of most species.

Intron retention can be seen as an artefact or a biologically relevant mechanism and in both cases it is interesting to monitor it in RNA-Seq data sets. This is easy when having a reference genome because one can quantify reads aligned on introns. This is more complicated in a de novo approach in which the assembler can produce contigs with and without intron for the same transcript.

Contigs with introns are more difficult to annotated because introns are splitting coding sequences (CDS) and possibly generating several shorter protein alignments rather than a unique long match.

Canonical splice motifs presence, number of reads splicing at the same location, the minimum number of splice events and overlaps (XXXXXXX) can be used to reduce faulty detection.

Hereafter we present IntronSeeker a software package enabling to search and remove introns from de novo assembled RNA-Seq reads.

# Statement of need

RNA-Seq de novo assemblies are widely used to study transcription in species lacking a referencegenome.  The classical extraction protocol collects RNA fragments using their PolyA tails and there-fore should only gather mature RNAs.  Unfortunately part of the extracted RNA molecules still com-prise  retained  introns  which  can  therefore  be  present  in  some  assembled  contigs.   These  intronsgenerate  transcript  frameshits  and  thus  render  annotation  more  difficult.    IntronSeeker  searchesintrons  by  splice-realigning  reads  on  contigs.    Introns  candidates  showing  canonical  intron/exonboundaries supported by several sequences.  IntronSeeker found between 0.02 and 11.96% of pu-tative  introns  in  17  publicly  available  RNA-Seq  de  novo  samples.   IntronSeeker  produces  an  in-tron cleaned contig fasta file as well as a cleaning report.  IntronSeeeker can be downloaded fromhttps://github.com/genotoul-bioinfo/intronseeker.


# Implementation

## General overview

IntronSeeker is a python script which includes three steps enabling to align read on contigs, find and remove retained introns and produce an html report. Figure 1. presents theses steps with input and output file formats.

IntronSeeker is open source (GNU General Public License) and can be download from <http://github.com/>

![IntronSeeker steps diagram](paper/Figures/IS_pipeline.png)


## Conda based installation procedure

To ease installation, intronSeeker includes an installation script (setup.sh) which run the installation of all the dependencies (1) but one (grinder) using conda and then installs grinder in the conda environment. Conda has to be installed beforehand. Installation can be checked using intronSeeker checkInstall program which will verify the presence and version of all dependencies. To run  intronSeeker one has first to load the ISeeker_environment using conda 'source activate ISeeker_environment' command. An example data set named reduced_real_dataset is also provided in the data directory. It includes the result files which will enable manual comparison of the reference and produced intron files after complying the test. The installation procedure and the test command line can be found on the main page of the intronSeeker GITHUB WEB-page.

\(1\) IntronSeeker dependencies include seven software packages :
grinder [@angly2012grinder], gffread [@gffread], hisat2 [@kim2015hisat], STAR [@dobin2013star], samtools [@li2009sequence], TransDecoder [@haas2013novo], diamond [@buchfink2015fast].


Python version 3.6 or above is required before installation:

```
git clone https://forgemia.inra.fr/emilien.lasguignes/intronSeeker.git
# Genologin cluster:
module load system/Miniconda3-4.7.10  

# Genobioinfo cluster:
module load devel/Miniconda/Miniconda3

#For a faster installation, we advise you to use mamba,  available in the latest version of conda:
conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba 

# Run intronSeeker
cd intronSeeker/ 
./setup.sh 

# Check installation
source activate ISeeker_environment 
intronSeeker checkInstall 
\end{lstlisting}

# Command line help
$ intronSeeker -h
```

## Setting default detection parameters

Detected introns are filtered according to thresholds (number of reads overlapping an intron and maximum length), canonical junction (GT_AG or CT_AC), and complexity (too long or overlapping introns).

Different GBS simulations where performed to find the best possible thresholds for minimum number of splice events, maximum intron spanning size, presence of canonical splice sites and minimum overlap size to call a splice event an intron retention candidate.

![GBS parameters impacts](paper/Figures/GBS_params_impacts.png)

IntronSeeker parameters impacts. Graphics have been built on 10 samples of Arabidopsis thaliana data with powerlow abundance (1.2).

**A** Coverage impact on detection: Increasing the coverage allows to quickly lose false candidates with a very limited impact on the number of true candidates (Arabidopsis thaliana data, powerlow abundance).\

**B** Filter candidates on border improve detection. DPratio =DPin/(DPbefore+DPafter)\

**C** Filter candidates on candidat length.\

**D** Filter candidates on retained intron ratio (long contig / long candidat).\


