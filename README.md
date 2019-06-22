# annotate-with-vcfanno
Configuration and data preparation instructions for annotation with vcfanno using identifiers
------------------------------------------------------------------------------------------------

This is a just a tiny example how to add annotations associated with identifiers
such as MIM number, OMIM disease, ORPHANETdisorder, Gene description or Gnomad OE Lof annotations
It uses bcbio-nextgen data and configuration. Currently this is integrated into bcbio-nextgen. 

How to try it
------------------

git clone
cd annotate-with-vcfanno
vcfanno -p n -base-path /path/to/bcbio-nextgen/genome/HS/build -lua cheop.lua cheop.conf NA12878.chr22.example.vcf  

Annotation data sources
-----------------------------------

For usage with all data the identifier mapping is based on NCBI gene EntrezID and gene symbol.
MIM and OMIM annotations downloaded from https://www.omim.org/downloads/
The gene2ensembl.gz identifiers and  gene_info.gz are from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
The ORPHANET rare diseases with their associated genes are from http://www.orphadata.org/data/xml/en_product6.xml
The gnomad lof constraint data downloaded from 
https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

Data preprocessing steps
---------------------------

OMIM mim2gene:

$cat mim2gene-2019.txt| grep gene | grep -v ^# | awk '{print $3"="$1}' | grep -v ^=  > mimTable

OMIM genemap2.txt:

$cat genemap2.txt | tr " " _ | grep -v ^# | tr "\t" = | sed 's/==/=NA=/g' | sed 's/==/=NA=/g' | tr "=" "\t" | cut -f6,8,10,11,13 | tr -d " " | awk '{print $1"= \"Gene: "$2" | EntrezIDom: "$3" | Ensg: "$4" | Phenotype: "$5"\""}' > genemap2Table

NCBI GENE gene_info.gz:

$zcat gene_info.gz | grep ^9606 | tr " " _ | cut -f2,3,9 | awk '{print $1"=\""$2": "$3"\""}' > gene_infoTable

NCBI GENE gene2ensembl.gz
$zcat gene2ensembl.gz | grep ^9606 | tr " " _ | cut -f2,3 | awk '{print $1"="$2}' |  tr "\n" ";"  | sed 's/;$//' > gene2ensemblTable

ORPHANET ( remove irrelevant lines and convert to UTF-8 charset that XML parser accepts, and parse XML into table) :

$tail -n +10 en_product6.xml | head -n -1 | dos2unix -iso | iconv -f ISO-8859-15 -t UTF-8 > test_en_product6.xml
$perl orphanetXMLparse.pl test_en_product6.xml > orphanetTable          # rows: gene disease association pairs
$perl orphanetXMLparse_unique.pl test_en_product6.xml > orphanetTable2  # rows: unique gene name and its associated list of diseases

gnomAD lof by gene:

$mv gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz gnomad.v2.1.1.lof_metrics.by_gene.txt.gz
$gunzip gnomad.v2.1.1.lof_metrics.by_gene.txt.gz
$cut -f1,5,24,64 | grep -v gene | awk '{print $1"=\"gnomad_oe_mis: " $2" | gnomad_oe_lof: " $3 " | "$4"\""}' > gnomadlofTable
