#!/bin/bash

#<IntronSeeker searches introns by splice-realigning reads on contigs.>
#Copyright (C) <2019-2024> INRAE
#<Sarah Maman, Philippe Bardou, Emilien Lasguignes, Faustine Oudin, Floréal Cabanettes, Christophe Klopp>
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Example of run:
#  ./complete_candidat_file.sh  srs_Ggal_d2_sample_1_HISAT2_candidates.txt  gbs_Ggal_d2_sample_1_transcripts-modified.fa    gbs_Ggal_d2_sample_1_transcripts-modified.gtf path/to/data path/to/outputDir


# Create outpir directory and path
dirPath=$4
outputDirName=$5

cd $4;
mkdir $5/;
cd $4/$5/;

# Get fasta length
~sigenae/bin/fasta_length.pl < $4/$2  > fasta.length

# Get correct intron information
cut -f 1,3,4,5 $4/$3  | awk '{print $1"|"$3"|"$4"\t"$2}' > introns.txt
# Copy candidate file
cp $4/$1  candidates.txt
# Merge candidate and length
awk 'NR==FNR{a[$1]=$2}NR!=FNR{if ($2 in a){print $0"\t"a[$2]}else{print $0"\tNone"}} ' fasta.length candidates.txt > contigs.length
# Merge candidate.length and introns
awk 'NR==FNR{a[$1]=$0}NR!=FNR{if ($1 in a){print $0"\t"a[$1]}else{print $0"\tFalse"}} ' introns.txt contigs.length  > candidates.Cglen.intron.txt

# Depth of false candidates
grep False candidates.Cglen.intron.txt | sed '1,3d' | cut -f 5 | sort -n | uniq -c > depth_False_candidats.out;
sed -i -e '1i\,False,Cov'  depth_False_candidats.out; sed 's/ \+/,/g' depth_False_candidats.out > depth_False_candidats.csv;
# Remove duplicated column 
cut -f 1,2,3,4,5,6,7,8,9,10,11,13 candidates.Cglen.intron.txt > candidates.Cglen.ri.txt

# Add candidat lenght
grep -v '^#' candidates.Cglen.intron.txt | awk -F '\t' -v col13=13 -v col4=4 -v col3=3 '{$col13=$col4-$col3;}1' - | awk '$(NF+1)="lCd"' -  > candidates.Cglen.intron.Cdlen.txt;

# Add DPin/(DPbefore+DPafter)
sed -i 's/ /\t/g' candidates.Cglen.intron.Cdlen.txt;
awk '{a=0;b=$9+$7;if(b!=0){a=$8/b;};print $0,a;}' candidates.Cglen.intron.Cdlen.txt | awk '$(NF+1)="DP"' - > candidates.Cglen.intron.Cdlen.DPratio.txt

# Replace "retained_intron" by "True"
sed -i 's/retained_intron/True/g'  candidates.Cglen.intron.Cdlen.DPratio.txt;
sed 's/ /,/g' candidates.Cglen.intron.Cdlen.DPratio.txt > candidates.Cglen.intron.Cdlen.DPratio.csv  #Convert in csv
sed -i -e '1i\contig,candidat,start,stop,depth,border,nb,nb,nb,filter,lcontig,Candidat,lCd,Critere_lCd,DPratio,Critere_DP'  candidates.Cglen.intron.Cdlen.DPratio.csv;   # Add title

# Add RI (= long contig / long candidat)
sed -i 's/ /\t/g' candidates.Cglen.intron.Cdlen.DPratio.txt ;  awk '$(NF+1)=$13/$11' candidates.Cglen.intron.Cdlen.DPratio.txt  | awk '$(NF+1)="RI"' -  > candidates.Cglen.intron.Cdlen.DPratio.RI.txt
sed 's/ /,/g'  candidates.Cglen.intron.Cdlen.DPratio.RI.txt > candidates.Cglen.intron.Cdlen.DPratio.RI.csv;
sed -i -e '1i\contig,candidat,start,stop,depth,border,nb,nb,nb,filter,lcontig,Candidat,lCd,Critere_lCd,DPratio,Critere_DP,RI,Critere_RI' candidates.Cglen.intron.Cdlen.DPratio.RI.csv;

# Add info about canonical junctions
egrep -v "(CT_AC|GT_AG)" candidates.Cglen.intron.Cdlen.DPratio.RI.csv  | awk '$(NF+1)="0"' - > tmp1
egrep "(CT_AC|GT_AG)"  candidates.Cglen.intron.Cdlen.DPratio.RI.csv | awk '$(NF+1)="1"' - > tmp2
cat tmp1 tmp2 > tmp3
awk '$(NF+1)="Junctions"' tmp3 | grep -v "contig,candidat,start,stop,depth,border,nb,nb,nb,filter" - > candidates.Cglen.intron.Cdlen.DPratio.RI.JC.csv
rm -rf tmp1 tmp2 tmp3;
sed -i 's/ /,/g' candidates.Cglen.intron.Cdlen.DPratio.RI.JC.csv;
sed -i -e '1i\contig,candidat,start,stop,depth,border,nb,nb,nb,filter,lcontig,Candidat,lCd,Critere_lCd,DPratio,Critere_DP,RI,Critere_RI,JC,Critere_JC' candidates.Cglen.intron.Cdlen.DPratio.RI.JC.csv;

## Supplementary analysis

# count lines
wc fasta.length introns.txt candidates.txt  contigs.length  candidates.Cglen.intron.txt  > count_lines.out

# Remove unuseful file(s)
rm -rf candidates.Cglen.ri.txt  candidates.Cglen.intron.txt candidates.Cglen.intron_min*.txt