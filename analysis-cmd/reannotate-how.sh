#!/bin/sh
# input will be list of directories
# cd into each dir and work with file
# ${name}_batch-gatk-haplotype-annotated-decomposed.vcf.gz
# IMPORTANT! Requires that vcfanno and vcf2db.py and gemini and sqlite3 should be in your system

vcfgzfile=$1
echo $vcfgzfile
currdir=`pwd`

workdir=`dirname $vcfgzfile`
echo $workdir
cd $workdir

name=`basename $vcfgzfile _batch-gatk-haplotype-joint-annotated-decomposed.vcf.gz`
# database will be name.db
database=${name}.db

echo "start" $name 
echo "start" $name >> $currdir/timelog
date >> $currdir/timelog

# reanotate vcf with vcfanoo
vcfanno -p 60 -lua reannot.lua reannot.conf ${vcfgzfile} | sed 's/Number=A/Number=1/g' > ${name}.vcf

#create  db with vcf2db.py 
# since it is only one sample at the moment, create ped file on the fly

echo "" >  ${name}.ped
sed -i '1i#Family_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\tPhenotype' ${name}.ped
sed -i '2iFAMILY\tSAMPLE\t-9\t-9\t0\t0' ${name}.ped
cat ${name}.ped | sed "s/FAMILY/${name}fam/" | sed "s/SAMPLE/${name}/"  > SAMPLE.ped
mv  SAMPLE.ped ${name}.ped
cat ${name}.ped

rm ${name}.db

echo vcf2db.py
vcf2db.py ${name}.vcf ${name}.ped ${name}.db

echo gemini query
# query genotypes and import the genotypes table to the databaset
gemini query --header -q "select variant_id, gts.${name},  gt_depths.${name}, gt_ref_depths.${name}, gt_alt_depths.${name}, gt_quals.${name}, gt_alt_freqs.${name} from variants" ${database} | sed "s/\.${name}//g"  > gts.csv

echo import annotations
sqlite3 ${database} ".mode tabs" ".headers on" "drop table if exists genotypes;" ".import gts.csv genotypes" 
sqlite3 ${database} < $currdir/report.sql > $name-report-reannot.csv.txt
# send rezult to some location
# cp $name-report-reannot.csv.txt /data/some-location/reannot/

echo "finished ${name}
echo "finished ${name} >> $currdir/timelog
date >> $currdir/timelog

cd $currdir
