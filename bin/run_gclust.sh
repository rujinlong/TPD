#!/bin/bash

QRY=$1
OUT=$2
MEM=$3
ncpu=$4
sortgenome.pl --genomes-file $QRY --sortedgenomes-file sort_$QRY
gclust -both -nuc -threads $ncpu -memiden $MEM sort_$QRY > clust_${QRY%.*}.txt
pretty_cdhit.py -i clust_${QRY%.*}.txt -o ${QRY}.map
cut -f1 ${QRY}.map | sed '1d' | sort -u > ${QRY}.list
seqkit grep -f ${QRY}.list $QRY > $OUT

#rm sort_$QRY clust_$QRY ${QRY}.list
