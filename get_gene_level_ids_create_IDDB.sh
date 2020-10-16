#!/bin/bash
# This script REQUIRES convert2bed by A.Reynolds. To install: sudo apt install bedops.
# This script installs compbio toolkit written by O.Buske. 
#
# Reference genome Homo sapiens (human)assembly GRCh38.p13) 
# Download the reference genome annotation in gff format from NCBI
# To convert gff to bed use convert2bed .
# To have gene coordinates the NC_ seqID have to be recoded into the chromosomes.  

# Once we have coordinates, we need to match the GeneID in the annotation with the 
# GeneID in OMIM,ORPHANET and gene_info . Perhaps the most efficient way to do that 
# is to import the tables into the sqlite and  crossreference them 
# using the query language  and print the corresponding table to use in vcfanno . 

currentd=`pwd`
OMIMD=`pwd`/OMIM
echo "OMIMD $OMIMD"

[ ! -d `pwd`/tmp ] && mkdir -p `pwd`/tmp && cd tmp
currentdtmp=`pwd`


####################################################################################################
## Part1 
## get script for ORPHANET data prep if needed
## get bedops tools if not exist 
####################################################################################################

# clone comp-bio to use orphanet processing script

if [ ! -d `pwd`/scripts/compbio-toolkit ]; then 
  git clone https://github.com/buske/compbio-toolkit.git 
  mv compbio-toolkit `pwd`/scripts/
fi

# get bedops if not exist
if [ ! -d `pwd`/scripts/bedops ]; then 
   wget https://github.com/bedops/bedops/releases/download/v2.4.36/bedops_linux_x86_64-v2.4.36.tar.bz2
   tar jxvf bedops_linux_x86_64-v2.4.36.tar.bz2
   mv bin `pwd`/scripts/bedops
   rm bedops_linux_x86_64-v2.4.36.tar.bz2
      
fi
# export path to be able to use scripts this single time
export PATH=`pwd`/scripts/bedops:$PATH
echo $PATH

####################################################################################################
## Part2 - download or read sources into source-data
##         all spaces transform to _ , columns separator TAB
####################################################################################################
echo " Source data to source-data dir "

if [ ! -d `pwd`/source-data ]; then
  mkdir source-data; cd `pwd`/source-data
else
  cd `pwd`/source-data
  rm *
fi

## Working directory 
WD=`pwd`


#####################################################################################################
## RefSeq GRCh38.p13 RELEASE  annotations 
#####################################################################################################
## read in a GFF table
if [ ! -f GCF_000001405.39_GRCh38.p13_genomic.gff.gz ] ;
then
 wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
fi
echo "create NC_ id"
zcat GCF_000001405.39_GRCh38.p13_genomic.gff.gz | grep -v ^# | cut -f1 | grep "NC_00" | uniq  >NC_id
cat NC_id

echo "create tmp.gff leave only genes"
zcat GCF_000001405.39_GRCh38.p13_genomic.gff.gz | grep -v ^# | grep ^"NC_00" | awk '{if ($3=="gene") print $0}' > tmp.gff
head tmp.gff

echo "replace NC_ to chr"
for string in `cat NC_id` 
do 
 replacement=`echo $string | sed 's/\..*$//' | sed 's/NC_0\{1,\}/chr/'| sed 's/chr23/chrX/' | sed 's/chr24/chrY/' `
 echo "replacing $string to $replacement" 
 sed -i "s/^${string}/${replacement}/" tmp.gff
done 
echo "convert2bed  tmp.gff to GCF_000001405.39_GRCh38.p13_genomic.gff.bed"
convert2bed -i gff < tmp.gff > GCF_000001405.39_GRCh38.p13_genomic.gff.bed

echo "Create coordinate file with GeneID column"
cat GCF_000001405.39_GRCh38.p13_genomic.gff.bed | awk '{print $10}' | sed 's/^.*GeneID/GeneID/' | sed 's/[,;].*$//' | sed 's/GeneID://' | awk '{ print FNR "\t" $0 }' > GeneID
awk '{ print FNR "\t" $0 }'  GCF_000001405.39_GRCh38.p13_genomic.gff.bed | tr " " _  > gcf
join GeneID gcf | tr " " "\t" > tmp.gff

# Add header to coordinates file
echo 'num\tGeneID\tchrom\tstart\tend\tname\tscore\tstrand\tsource\tregion\treserved\tinfo'>tmp.header
cat tmp.header tmp.gff > hg38_genes_refseq_coordinates.tsv 

echo " ==========================================With coordinates GRCh38.p13 RefSeq  GCF_000001405.39_GRCh38.p13_genomic.gff"
echo ""
head -n3 hg38_genes_refseq_coordinates.tsv 
echo ""
tail -n3 hg38_genes_refseq_coordinates.tsv 
echo ""

#####################################################################################################
## Gene history changed and discontinued gene identifiers from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
#####################################################################################################
## 
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz
zcat gene_history.gz | grep ^9606 | cut -f2-5 > hs_gene_history.tsv 
sed -i '1i\GeneID\tDiscontinued_GeneID\tDiscontinued_Symbol\tDiscontinue_Date' hs_gene_history.tsv

echo " =======================================================================Discontinued genes gene_history"
echo ""
head -n10 hs_gene_history.tsv
echo ""

#####################################################################################################
## Full HUGO data set for gene names 
#####################################################################################################
## 

echo "HUGO gene names complete set https://www.genenames.org/download/statistics-and-files/"
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt

## make it into tsv
cat hgnc_complete_set.txt | tr " " _  >  hgnc.tsv

echo " =======================================================================HUGO gene names.tsv"
echo ""
head -n4 hgnc.tsv
echo ""

#####################################################################################################
## Obtain orphanet gene table 
#####################################################################################################
python3 ${WD}/../scripts/orphanet/get-orphanet-gene-xls.cgi | tail -n +5 | tr " " "_" | tr ";" ","  |   sed 's/P12081.*antigen/P12081,J_antigen/g'  > orpha.tsv

echo " =======================================================================orpha.tsv"
echo ""
head -n4 orpha.tsv
echo ""

#####################################################################################################
## OMIM data files
#####################################################################################################
#OMIMD=/projects/DB/OMIM
## OMIM files location
## OMIM 2019 has  mim2gene.txt, genemap2.txt,  morbidmap.txt , mimTitles.txt


# if OMIM folder exists clean files to import to sqlite3
if [ -d ${OMIMD} ]; then 

##############################################################################mim2gene#########
cat ${OMIMD}/mim2gene.txt | sed 's/^# MIM Number/MIM Number/' |  grep -v ^# | tr " " _  > mim2gene.tsv

echo " =======================================================================mim2gene.tsv"
echo ""
head -n4 mim2gene.tsv
echo ""

##############################################################################genemap2########
cat ${OMIMD}/genemap2.txt |  sed 's/^# Chromosome/Chromosome/' | grep -v ^# | tr " " _ | tr ";" "," | sed '$d' > genemap2.tsv

echo " =======================================================================genemap2.tsv"
echo ""
head -n4 genemap2.tsv
echo "tail"
tail -n4 genemap2.tsv
echo ""

##############################################################################morbidmap.txt###
cat ${OMIMD}/morbidmap.txt | sed 's/^# Phenotype/Phenotype/' | grep -v ^# | tr " " _ | tr ";" "," | sed '$d' | sed '$d' > morbidmap.tsv

echo " =======================================================================morbidmap.tsv"
echo ""
head -n4 morbidmap.tsv
echo "tail"
tail -n4 morbidmap.tsv
echo ""

##############################################################################mimTitles.txt###
cat ${OMIMD}/mimTitles.txt | sed 's/^# Prefix/Prefix/' | grep -v ^# | tr " " _ | tr ";" "," > mimTitles.tsv

echo " =======================================================================mimTitles.tsv"
echo ""
head -n4 mimTitles.tsv
echo ""

fi

#####################################################################################################
## Full HPO disease-gene-phenotype table from charite
#####################################################################################################
## read in HPO association table

echo "HPO all frequencies all sources HPO http://compbio.charite.de/   genes_to_phenotypes.txt"
if [ ! -f genes_to_phenotype.txt ]; then
  wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/genes_to_phenotype.txt
fi

if [ ! -f phenotype_to_genes.txt ]; then 
  wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/phenotype_to_genes.txt
fi

## change the header 
cat genes_to_phenotype.txt | tr " " _ | sed 's/#Format:_//' | sed 's/<tab>/\t/g' >hpo_genes_to_phenotype.tsv
cat phenotype_to_genes.txt | tr " " _ | sed 's/#Format:_//' | sed 's/<tab>/\t/g' >hpo_phenotype_to_genes.tsv

echo " =======================================================================hpo.tsv"
echo ""
head -n4 hpo_genes_to_phenotype.tsv
echo ""
echo ""
head -n4 hpo_phenotype_to_genes.tsv
echo ""

#####################################################################################################
## Obtain gene2ensemble mapping file
#####################################################################################################
echo "gene2ensembl.tsv"
# Process gene2ensembl 
    wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
    zcat gene2ensembl.gz | grep ^9606 | cut -f2-7 | sort -k1,1 | uniq > gene2ensembl.tsv
# Add header 
    sed -i '1 i\GeneID\tEnsembl_gene_identifier\tRNA_nucleotide_accession\tEnsembl_rna_identifier\tprotein_accession\tEnsembl_protein_identifier' gene2ensembl.tsv     

echo " =======================================================================gene2ensembl.tsv"
echo ""
head gene2ensembl.tsv
echo ""


#####################################################################################################
## Obtain gene_info file
#####################################################################################################
echo "gene_info.tsv"
# Process gene_info
  wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
  zcat gene_info.gz | grep ^9606 | tr " " _ | cut -f2,3,5,6,9,10,11,13,14,15  | sort -k1,1 | uniq | tr ";" ","  > gene_info.tsv
  sed -i '1 i\GeneID\tSymbol\tSynonyms\tdbXrefs\tdescription\ttype_of_gene\tSymbol_from_nomenclature_authority\tNomenclature_status\tOther_designations\tModification_date' gene_info.tsv
 
echo " =======================================================================gene_info.tsv"
echo ""
head gene_info.tsv
echo ""
echo ""

#####################################################################################################
## GENCODE  v33 data , add headers to files
#####################################################################################################
#####################################################################################################
## Obtain enst_to_geneid file
#####################################################################################################
echo "enst_to_geneid.tsv"
# Process gene_info
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.metadata.EntrezGene.gz
  zcat gencode.v33.metadata.EntrezGene.gz |sed 's/\.[0-9]\+//g' > enst_to_geneid.tsv
  sed -i '1 i\enst\tGeneID' enst_to_geneid.tsv
echo " =======================================================================enst_to_geneid.tsv"
echo ""
head enst_to_geneid.tsv
echo ""
echo ""
#####################################################################################################
## Obtain enst_to_symbol file
#####################################################################################################
echo "enst_to_symbol.tsv"
# Process HGNC
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.metadata.HGNC.gz
  zcat gencode.v33.metadata.HGNC.gz |sed 's/\.[0-9]\+//g'| sed 's/HGNC://' > enst_to_symbol.tsv
  sed -i '1 i\enst\tsymbol\thgnc' enst_to_symbol.tsv
echo " =======================================================================enst_to_symbol.tsv"
echo ""
head enst_to_symbol.tsv
echo ""
echo ""
#####################################################################################################
## Obtain enst_to_refseq file
#####################################################################################################
echo "enst_to_refseq.tsv"
# Process HGNC
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.metadata.RefSeq.gz
  zcat gencode.v33.metadata.RefSeq.gz |sed 's/\.[0-9]\+//g' > enst_to_refseq.tsv
  sed -i '1 i\enst\tNM\tNP' enst_to_refseq.tsv
echo " =======================================================================enst_to_refseq.tsv"
echo ""
head enst_to_refseq.tsv
echo ""
echo ""
#####################################################################################################
## Obtain enst_to_swissprot file
#####################################################################################################
echo "enst_to_refseq.tsv"
# Process HGNC
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.metadata.SwissProt.gz
  zcat gencode.v33.metadata.SwissProt.gz |sed 's/\.[0-9]\+//g' > enst_to_swissprot.tsv
  sed -i '1 i\enst\tuniprot1\tuniprot2' enst_to_swissprot.tsv
echo " =======================================================================enst_to_swissprot.tsv"
echo ""
head enst_to_swissprot.tsv
echo ""
echo ""

# GENCODE gene annotations, leaving out for now
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz
# zcat gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz | head


####################################################################################################
## Part3 - make SQLITE3 database of all identifiers
##         Import tsv using command line
##         Time stamp source data and IDDB.db
##         This database will ten be used to create links to gene level annotations
####################################################################################################

if [ -f IDDB.db ]; then
  rm IDDB.db
fi

sqlite3 IDDB.db ".mode tabs" ".import hgnc.tsv hgnc" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import hg38_genes_refseq_coordinates.tsv hg38" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import orpha.tsv orpha" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import mim2gene.tsv mim2gene" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import genemap2.tsv genemap2" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import morbidmap.tsv morbidmap" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import mimTitles.tsv mimTitles" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import hpo_genes_to_phenotype.tsv hpog2p" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import hpo_phenotype_to_genes.tsv hpop2g" ".exit"

sqlite3 IDDB.db ".mode tabs" ".import gene2ensembl.tsv gene2ensembl" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import gene_info.tsv geneinfo" ".exit"

sqlite3 IDDB.db ".mode tabs" ".import hs_gene_history.tsv genehistory" ".exit"

# GENCODE v33 identifiers
sqlite3 IDDB.db ".mode tabs" ".import enst_to_geneid.tsv enst2geneid" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import enst_to_refseq.tsv enst2refseq" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import enst_to_symbol.tsv enst2symbol" ".exit"
sqlite3 IDDB.db ".mode tabs" ".import enst_to_swissprot.tsv enst2swissprot" ".exit"


  echo ""
  echo "Tables in IDDB.db: "
  echo ""

sqlite3 IDDB.db ".tables" ".exit"

# move IDDB.db
mv IDDB.db  ${currentd}
cd ${currentd}

# proceed to create annotation tables
# use sql script on IDDB.db
sqlite3 IDDB.db < ./omim-orpha-ensg-annotations.sql
exit
