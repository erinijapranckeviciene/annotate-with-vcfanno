#!/bin/bash
FILELIST=$1
# generate filelist of vcf.gz from which to compute variant frequency tables 
# ls /data/reports/[1,S]*/*gatk*-annotated-decomposed.vcf.gz > FILELIST 
for f in `cat $FILELIST`; do bcftools query -f '[%GT\t%CSQ=]\n' ${f} | tr "|" "\t" | cut -f1,46 | awk '{print $2 "\t" $1}' | grep chr | grep "/" | sed 's/0\/1/het/' | sed 's/1\/1/hom/' ; done > d.txt
# file d.txt contains variants from the listed variant vcf.gz files 
# count variant frequencies
cat d.txt | sort -k1,1 | datamash -g 1,2 count 1 |  datamash -g 1 sum 3 collapse 2 collapse 3 | awk '{print $1 "\t" $2 "\t" $2/150 "\t" $3 "\t" $4}' > dazniai.txt 
# create lua table
cat dazniai.txt | awk '{print $1"=Count:"$2"|Freq:"$3"|Zyg:"$4"|ZygN:"$5}' 
# > /projects/Annotations/vcfanno/tables/localvariants.luatable 
