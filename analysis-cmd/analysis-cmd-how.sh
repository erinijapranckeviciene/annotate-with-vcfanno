#!/bin/sh
fastqlist=$1

dline=`head -n1 ${fastqlist}`
dirname=`echo ${dline%/*}`
echo $dirname
currentdir=`pwd`
echo $currentdir
for f in `cat ${fastqlist}`; do basename $f ".fastq.gz"; done | sed 's/\..*//' | sed 's/_.*//' | uniq > projectlist
cat projectlist

for f in `cat projectlist`
 do 
  echo "samplename,description,batch,phenotype,sex,variant_regions" >${f}.csv
  echo "${f},${f},${f}_batch,,," >> ${f}.csv

  cat ${f}.csv 

  bcbio_nextgen.py -w template /projects/bcbio/TwistExome/exome_hotc.yaml ${f}.csv ${dirname}/${f}.*

  cat ${f}/config/${f}.yaml
  cd `pwd`/${f}/work
  ls ../config/* 
  bcbio_nextgen.py ../config/${f}.yaml -n32  2> log1>&1
  cd  ${currentdir}
 done
#!/bin/sh
projectlist=$1
wd=`pwd`

for f in `cat projectlist`
 do 
  echo "copy ${f} to reports"
  [ ! -d /data/reports/${f} ] && mkdir -p /data/reports/${f}
  cp -R ${wd}/${f}/final/20*/${f}* /data/reports/${f}/
  cp -R ${wd}/${f}/final/20*/multiqc /data/reports/${f}/

  cp -R ${wd}/${f}/final/${f}/* /data/reports/${f}/
  ls /data/reports/${f}
  echo "done"
  echo ""
 done
