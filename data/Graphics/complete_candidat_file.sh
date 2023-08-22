#!/bin/bash

# ./complete_candidat_file_V3.sh   srs_Athal_pw-02_d1_HISAT2_candidates.txt  gbs_Athal_pw-02_d1_transcripts-modified.fa  gbs_Athal_pw-02_d1_transcripts-modified.gtf
#  ./complete_candidat_file.sh  srs_Ggal_d2_sample_1_HISAT2_candidates_HEAD20.txt  gbs_Ggal_d2_sample_1_transcripts-modified_HEAD2000.fa    gbs_Ggal_d2_sample_1_transcripts-modified_HEAD100.gtf

# get fasta length
~sigenae/bin/fasta_length.pl < $2  > fasta.length

# get correct intron information
cut -f 1,3,4,5 $3  | awk '{print $1"|"$3"|"$4"\t"$2}' > introns.txt
# copy candidate file
cp $1  candidates.txt
# merge candidate and length
awk 'NR==FNR{a[$1]=$2}NR!=FNR{if ($2 in a){print $0"\t"a[$2]}else{print $0"\tNone"}} ' fasta.length candidates.txt > contigs.length
# merge candidate.length and introns
awk 'NR==FNR{a[$1]=$0}NR!=FNR{if ($1 in a){print $0"\t"a[$1]}else{print $0"\tFalse"}} ' introns.txt contigs.length  > candidates.Cglen.intron.txt
# Remove duplicated column 
cut -f 1,2,3,4,5,6,7,8,9,10,11,13 candidates.Cglen.intron.txt > candidates.Cglen.ri.txt
# Add candidat lenght
awk -F '\t' -v col13=13 -v col4=4 -v col3=3 '{$col13=$col4-$col3;}1' candidates.Cglen.ri.txt > candidates.Cglen.ri.Cdlen.txt;
sed -i 's/ /\t/g'  candidates.Cglen.ri.Cdlen.txt;
grep  -v '^#' candidates.Cglen.ri.Cdlen.txt | sed 's/\t/,/g' - > candidates.Cglen.ri.Cdlen_withoutHeader.csv;
 
# Pourcent ir/cg
# sed -i '/^##/d' candidates.Cglen.ri.Cdlen.txt   #Remove header comment
# awk -v col14=14 -v col13=13 -v col11=11 '{($col11 ? $col14=($col13/$col11)*100 : 0);}1'  candidates.Cglen.ri.Cdlen.txt > candidates.Cglen.ri.Cdlen_ircg.txt
# sed -i 's/ /\t/g' candidates.Cglen.ri.Cdlen_ircg.txt;

# supplementary analysis

# count lines
wc fasta.length introns.txt candidates.txt  contigs.length  candidates.Cglen.intron.txt  > count_lines.out

# Filter on depth d >=1
#for i in {1..5}
#do
#    awk '$5>="$i"' candidates.Cglen.ri.Cdlen.txt > candidates.Cglen.ri.Cdlen_min"$i".txt;
#    awk '$5>="$i"' candidates.Cglen.intron.txt > candidates.Cglen.intron_min"$i".txt
#done

awk '$5>=0' candidates.Cglen.ri.Cdlen.txt | grep -v '^#' - | sed 's/\t/,/g'  - > candidates.Cglen.ri.Cdlen_min0_withoutHeader.csv;
awk '$5>=1' candidates.Cglen.ri.Cdlen.txt | grep -v '^#' - | sed 's/\t/,/g'  - > candidates.Cglen.ri.Cdlen_min1_withoutHeader.csv;
awk '$5>=2' candidates.Cglen.ri.Cdlen.txt | grep -v '^#' - | sed 's/\t/,/g'  - > candidates.Cglen.ri.Cdlen_min2_withoutHeader.csv;
awk '$5>=3' candidates.Cglen.ri.Cdlen.txt | grep -v '^#' - | sed 's/\t/,/g'  - > candidates.Cglen.ri.Cdlen_min3_withoutHeader.csv;
awk '$5>=4' candidates.Cglen.ri.Cdlen.txt | grep -v "#"  - | sed 's/\t/,/g'  - > candidates.Cglen.ri.Cdlen_min4_withoutHeader.csv;
awk '$5>=5' candidates.Cglen.ri.Cdlen.txt | grep -v "#"  - | sed 's/\t/,/g'  - > candidates.Cglen.ri.Cdlen_min5_withoutHeader.csv;


awk '$5>=0' candidates.Cglen.intron.txt > candidates.Cglen.intron_min0.txt
awk '$5>=1' candidates.Cglen.intron.txt > candidates.Cglen.intron_min1.txt
awk '$5>=2' candidates.Cglen.intron.txt > candidates.Cglen.intron_min2.txt
awk '$5>=3' candidates.Cglen.intron.txt > candidates.Cglen.intron_min3.txt
awk '$5>=4' candidates.Cglen.intron.txt > candidates.Cglen.intron_min4.txt
awk '$5>=5' candidates.Cglen.intron.txt > candidates.Cglen.intron_min5.txt

# output file format
# depth of false candidates
grep False candidates.Cglen.intron.txt | sed '1,3d' | cut -f 5 | sort -n | uniq -c > depth_False_candidats.out;
sed -i -e '1i\,False,Cov'  depth_False_candidats.out; sed 's/ \+/,/g' depth_False_candidats.out > depth_False_candidats.csv;
# good intron depth
grep -v False candidates.Cglen.intron.txt | sed '1,3d' | cut -f 5 | sort -n | uniq -c > depth_true_candidats.out; 
sed -i -e '1i\,True,Cov' depth_true_candidats.out; sed 's/ \+/,/g' depth_true_candidats.out > depth_true_candidats.csv
# Nb junctions - With/without junctions
grep -vwE "(CT_AC|GT_AG)" candidates.Cglen.ri.Cdlen.txt  > candidates.Cglen.ri.Cdlen.WithoutJunctions.txt; sed -i 's/ /\t/g'  candidates.Cglen.ri.Cdlen.WithoutJunctions.txt;

for i in {0..5}
do
    grep False candidates.Cglen.intron_min"$i".txt | sed '1,3d' | cut -f 5 | sort -n | uniq -c  > depth_False_candidats_min"$i".out
    sed -i -e '1i\,False,Cov'  depth_False_candidats_min$i.out; sed 's/ \+/,/g' depth_False_candidats_min$i.out > depth_False_candidats_min$i.csv; rm -rf depth_False_candidats*.out
    grep -v False candidates.Cglen.intron_min"$i".txt | sed '1,3d' | cut -f 5 | sort -n | uniq -c > depth_true_candidats_min"$i".out
    sed -i -e '1i\,True,Cov' depth_true_candidats_min$i.out; sed 's/ \+/,/g' depth_true_candidats_min$i.out > depth_true_candidats_min$i.csv; rm -rf depth_true_candidats*.out
    grep -vwE "(CT_AC|GT_AG)" candidates.Cglen.ri.Cdlen_min"$i"_withoutHeader.csv  | sed  's/,/\t/g'  - > candidates.Cglen.ri.Cdlen.WithoutJunctions_min"$i".txt;
done


# Remove unuseful file(s)
rm -rf candidates.Cglen.ri.txt  candidates.Cglen.intron.txt candidates.Cglen.intron_min*.txt
