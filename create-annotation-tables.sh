#!/usr/bin/bash
# Create  tables with annotations from  gene_info.gz and gene2ensemble.gz, omim, orphanet and gnomad constraint
# Tables crossrefernce gene symbol  ans ENSG identifiers with NCBI gene EntrezID
# Check requirements : datamash , dos2unix, iconv, perl

####################################################################################################
## Part1 - collect the data
####################################################################################################
echo " Source data prep"

if [ ! -d `pwd`/source-data ]; then
  mkdir source-data; cd source-data
else
  cd `pwd`/source-data
fi

WD=`pwd`

echo "gene_info.tsv"
#  Process gene_info 
if [ ! -f gene_info.tsv ]; then

  if [ ! -f gene_info.gz ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
  fi
  zcat gene_info.gz | grep ^9606 | tr " " "_" | cut -f2,3,5,6,9,10 | sort -k1,1 > gene_info.tsv

fi
echo ""
head gene_info.tsv
echo ""

echo "gene2ensembl-uniq.tsv"
# Process gene2ensembl 
if [ ! -f gene2ensembl-uniq.tsv ]; then

  if [ ! -f gene2ensembl.gz ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
  fi
  zcat gene2ensembl.gz | grep ^9606 | cut -f2,3,4,5,7 > gene2ensembl.tsv
  cut -f1,2 gene2ensembl.tsv | sort -k1,1 | uniq  > gene2ensembl-uniq.tsv
fi 
echo ""
head gene2ensembl-uniq.tsv
echo ""

# The gene_info_ensembl table is used to create OMIM , Orphanet and Gnomad constraint tables

echo "gnomad_constraint.tsv"
# Process gnomad constraint file 
if [ ! -f gnomad_constraint.tsv ]; then

  if [ ! -f gnomad.v2.1.1.lof_metrics.by_gene.txt.gz ]; then
    wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
    mv gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz gnomad.v2.1.1.lof_metrics.by_gene.txt.gz
  fi
  zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.gz | cut -f1,5,21,22,23,24,64 > gnomad_constraint.tsv

fi
echo ""
head gnomad_constraint.tsv
echo ""

echo "orphanet_unique_gene_disorder.tsv"
# Process Orphanet file 
if [ ! -f orphanet_unique_gene_disorder.tsv ]; then
  if [ ! -f en_product6.xml ]; then
    wget http://www.orphadata.org/data/xml/en_product6.xml
  fi
  tail -n +10 en_product6.xml | head -n -1 | dos2unix -iso | iconv -f ISO-8859-15 -t UTF-8 > modif_en_product6.xml
  perl ${WD}/../orphanetXMLparse.pl modif_en_product6.xml > orphanet_gene_disorder_pairs.tsv
  perl ${WD}/../orphanetXMLparse_unique.pl modif_en_product6.xml | sed 's/;/,/g' | sed 's/=/\t/' > orphanet_unique_gene_disorder.tsv

fi
echo ""
head "orphanet_unique_gene_disorder.tsv"
echo ""

# Process OMIM
# With OMIM restriction- if you have genemap2.txt file, copy it to the source-data directory
# Mim to gene file is freely downloadable
echo "mim2gene.tsv"
if [ ! -f mim2gene.tsv ]; then
  if [ ! -f mim2gene.txt ]; then
    wget https://www.omim.org/static/omim/data/mim2gene.txt
  fi
## Clean and make 4 column file: MIM number \t EntrezID \t Symbol \ ensg_om and sort by MIM number
  cat mim2gene.txt| grep -v ^# | grep gene | cut -f1,3,4,5 | sort -k2,2 | tr "\t" = | sed 's/=$/=NA/' | sed 's/==/=NA=/' | grep -v "==NA" | tr = "\t" | sort -k1,1 > mim2gene.tsv

fi
echo ""
head "mim2gene.tsv"
echo ""

echo "genemap2.tsv"
## in genemap2.txt there are unexpected = , but symbol &  is not used, we use it to transform and clean the data
if [ ! -f genemap2.tsv ]; then
  if [ -f genemap2.txt ]; then
    cat genemap2.txt| tr " " "_" | grep -v ^# | cut -f6-8,10,11,13 | tr "\t" "&" | sed 's/\&$/\&NA/' | sed 's/\&\&/\&NA\&/' | sed 's/\&\&/\&NA\&/' | tr ";" ","  | tr "&" "\t" |  sort -k1,1 > genemap2.tsv
  fi
# posible error echo "NO genemap2.txt, OMIM phenotypes will be limited to only MIM number, or copy genemep2.txt to source-data"
fi
echo ""
head "genemap2.tsv"
echo ""

# To keep identifiers clean and tidy - the Orphanet, Omim and Gnomad constraint will be queried by  gene Entrez ID  only
# join the files with gene_info, keep all information from file 1, so as to know what identifiers are missing

echo "Finished source data prep"

####################################################################################################
## Part2 - create harmonized identifiers and annotations - annotations referenced by gene EntrezID
##  In annotations the ; is used to separate annotations, annotation can be 
##   multiple parts separated by | and  ,. 
####################################################################################################
echo ""
echo "Creating helper tables for annotation tables"
echo ""
## Create  gene_info that cotains ENSG identifiers from gene2ensembl joining by  entrezID
## In  joining  make the data into two columns , the symbol that apparently is not used in gene_info is '&'
cat gene_info.tsv | sed 's/\t/\&/4g' > a.tsv
join -a 1 -1 1 -2 1 a.tsv gene2ensembl-uniq.tsv | tr " " "\t" | awk '{if (length($5)==0) print $0 "\t" "NA"; else print $0}' | datamash --group 1 unique 2-5 | tr ";" "," > a_gene_info_ensembl  #  change ";" into ","
echo ""
echo " joined gene_info.tsv and ensembl-uniq.tsv into a_gene_info_ensembl"
echo ""
head -n5 a_gene_info_ensembl

## After joining a_gene_info_ensembl has entrezID \t symbol \t synonyms \t gene_info \t ensg; this file contains not only protein coding genes
## but all genes that are in gene_info  


## Create a table that contains all gene symbols and their synomyms each per line sort by symbol
perl ${WD}/../inflate_gs_synonyms_in_gene_info.pl a_gene_info_ensembl | sort -k2,2 > a_gene_info_inflated
echo ""
echo " created a_gene_info_inflated from the  a_gene_info_ensembl"
echo ""
head -n5 a_gene_info_inflated


#######################################################################################################################
##  OMIM Phenotypes
#######################################################################################################################
## Not all mim2gene records can be crossreferenced  with the genemap2.txt records. Joining these two files as of June, 2019
## The 1333 genemap2.txt(of 2017)  records do not have matching mim2gene MIM numbers. 
## It seems that  mim2gene and genemap2 annotations should be  separated. 
## No processsing is required for mim2gene.tsv, since it has MIM number, EntrezID , Symbol and ENSG 
## 
## In genemap2.tsv a part of records(16249) have EntrezID and another part(262 records) do not have EntrezID, but have gene symbol
## The second part we will try to crossreference with a_gene_info_inflated
cat genemap2.tsv | awk '{print $4 "\t" $0}' | sort -k1,1 | grep -v ^NA  > b_genemap2_has_entrezid.tsv
cat genemap2.tsv | awk '{print $4 "\t" $0}' | sort -k1,1 | grep  ^NA  > b_to_curate_genemap2_no_entrezid.tsv
echo ""
echo " created b_genemap2_has_entrezid from genemap2.tsv"
echo ""
head -n5 b_genemap2_has_entrezid.tsv

#######################################################################################################################
##  ORPHANET disorders
#######################################################################################################################
## Join Orphanet records and gene_info_inflated by gene symbol
## we have to transform all spaces into _ on orphanet_unique_gene_disorder
cat orphanet_unique_gene_disorder.tsv  | tr " " "_" > c
## join orphanet disorders list with inflated gene_info by symbol to assign 
## EntrezID to the gene symbol associated with orpha
## There might be different synonyms gene symbols associated with different disorders, but  having same entrezID
## As in this case:
## 9968	MED12	"Blepharophimosis-intellectual_disability_syndrome,_MKB_type,_ORPHAID:_293707_|_FG_syndrome_type_1,_ORPHAID:_93932"
## 9968	OPA1	"Autosomal_dominant_optic_atrophy_plus_syndrome,_ORPHAID:_1215"
## Such annotations  have to be collapse , we will use datamash 

join -a 1 -1 1 -2 2 c  a_gene_info_inflated | tr " " "\t" |  awk '{print $3 "\t" $1 "\t" $2}' | sort -k1,1 | uniq > c1
cat c1 | grep -v ^[[:space:]] | grep -v ^[a-z,A-Z] | datamash --group 1 collapse 2 collapse 3 >  c_entrezid_symbol_orphanet.tsv
echo ""
echo " created c_entrezid_symbol_orphanet.tsv"
echo ""
head -n5 c_entrezid_symbol_orphanet.tsv

## 22 gene symbol IDs in ORPHA do not have associated EntrezID they are removed from annotation table now, to be curated later
##     1		C11ORF80	"Complete_hydatidiform_mole,_ORPHAID:_254688"
##     2		C11ORF95	"Ependymoma,_ORPHAID:_251636"
##     3		C12ORF57	"Temtamy_syndrome,_ORPHAID:_1777"
##     4		C12ORF65	"Combined_oxidative_phosphorylation_defect_type_7,_ORPHAID:_254930_|_Autosomal_recessive_spastic_paraplegia_type_55,_ORPHAID:_320375"
##     5		C15ORF41	"Congenital_dyserythropoietic_anemia_type_I,_ORPHAID:_98869"
##     6		C19ORF12	"Mitochondrial_membrane_protein-associated_neurodegeneration,_ORPHAID:_289560_|_Autosomal_recessive_spastic_paraplegia_type_43,_ORPHAID:_320370"
##     7		C8ORF37	"Bardet-Biedl_syndrome,_ORPHAID:_110"
##     8		C9ORF72	"Behavioral_variant_of_frontotemporal_dementia,_ORPHAID:_275864_|_Huntington_disease-like_syndrome_due_to_C9ORF72_expansions,_ORPHAID:_401901_|_Progressive_non-fluent_aphasia,_ORPHAID:_100070_|_Semantic_dementia,_ORPHAID:_100069"
##     9		MT-ATP6	"NARP_syndrome,_ORPHAID:_644_|_MT-ATP6-related_mitochondrial_spastic_paraplegia,_ORPHAID:_320360"
##    10		MT-ATP8	"Periodic_paralysis_with_later-onset_distal_motor_neuropathy,_ORPHAID:_397750"
##    11		MT-CYB	"Histiocytoid_cardiomyopathy,_ORPHAID:_137675"
##    12		MT-ND3	"Leber_plus_disease,_ORPHAID:_99718"
##    13		MT-ND6	"Leber_hereditary_optic_neuropathy,_ORPHAID:_104"
##    14		MT-TE	"Myopathy_and_diabetes_mellitus,_ORPHAID:_2596_|_Mitochondrial_myopathy_with_reversible_cytochrome_C_oxidase_deficiency,_ORPHAID:_254864"
##    15		MT-TF	"MELAS,_ORPHAID:_550_|_MERRF,_ORPHAID:_551"
##    16		MT-TK	"Mitochondrial_DNA-related_cardiomyopathy_and_hearing_loss,_ORPHAID:_1349_|_Maternally-inherited_diabetes_and_deafness,_ORPHAID:_225"
##    17		MT-TL1	"Kearns-Sayre_syndrome,_ORPHAID:_480_|_Mitochondrial_DNA-related_progressive_external_ophthalmoplegia,_ORPHAID:_663_|_Hypertrophic_cardiomyopathy_and_renal_tubular_disease_due_to_mitochondrial_DNA_mutation,_ORPHAID:_324525"
##    18		MT-TS1	"Palmoplantar_keratoderma-deafness_syndrome,_ORPHAID:_2202_|_Mitochondrial_non-syndromic_sensorineural_deafness,_ORPHAID:_90641_|_Mitochondrial_non-syndromic_sensorineural_deafness_with_susceptibility_to_aminoglycoside_exposure,_ORPHAID:_168609"
##    19		MT-TT	"Lethal_infantile_mitochondrial_myopathy,_ORPHAID:_254857"
##    20		MT-TV	"Mitochondrial_DNA-associated_Leigh_syndrome,_ORPHAID:_255210"
##    21		SCA32	"Spinocerebellar_ataxia_type_32,_ORPHAID:_276183"
##    22		SLC7A2-IT1	"Ravine_syndrome,_ORPHAID:_99852"
## 

#######################################################################################################################
##  GNOMAD constraint
#######################################################################################################################
## Join Gnomad constraint and gene_info_ensembl by ensg to map the GnomadLof to EntrezID
## joining by symbol leaves 632  records unjoined out of 20949 
## joining by ENSG leaves 1060  records unjoined 
## The  Symbol is better 
sort -k1,1 gnomad_constraint.tsv | grep -v  pLI > d
sort -k2,2 a_gene_info_inflated > d1
## after joining  we have to collapse EntrezID . Column 8 - ENSG in GnomadConstraint, Column 9 is ENSG in Gene info
join -a 1 -1 1 -2 2 d d1 | tr " " "\t" | awk '{print $8 "\t" $0}' | sort -k1,1 | grep -v ^[[:space:]] | cut -f1-8,10 > d2
cat d2 | datamash --group 1 collapse 2-9 > d_gnomad_constraint
## After collapsing - situations like following needs investigation
##9962	SLC23A1,SLC23A2	7.0539e-01,4.8400e-01	2.3504e-02,6.2629e-03	1.1039e-04,3.1505e-05	9.7639e-01,9.9371e-01	2.9783e-01,3.1156e-01	ENSG00000170482,ENSG00000089057	ENSG00000089057,ENSG00000089057
##9963	SLC23A1,SLC23A2	7.0539e-01,4.8400e-01	2.3504e-02,6.2629e-03	1.1039e-04,3.1505e-05	9.7639e-01,9.9371e-01	2.9783e-01,3.1156e-01	ENSG00000170482,ENSG00000089057	ENSG00000170482,ENSG00000170482
echo ""
echo " creted d_gnomad_constraint"
echo ""
head -n5 d_gnomad_constraint


####################################################################################################
## Part3 - create  Tables as EntrezID = " annotation list "
##  The ; is used to separate  annotations in vcfanno 
##  Multiple fields in OMIM, ORPHANET , GNOMAD constraint , GENE description  are separated 
##  by | ,  the general format for multiple fields and values is name: value. 
##
##  Sources : c_entrez_symbol_orphanet.tsv to OrphanetTable
##            d_gnomad_constraint     to GnomadConstraintTable
##            a_gene_info_ensembl      to GeneDescriptionTable
##            b_genemap2_has_entrezid.tsv  to Genemap2Table
##            mim2gene.tsv             to Mim2geneTable
####################################################################################################

tables="../annotation-tables"
if [ ! -d ${tables} ]; then
  mkdir -p ${tables}
fi

##  Annotation tables will be save in annotation-tables
## OrphanetTable
cat c_entrezid_symbol_orphanet.tsv | tr  \" \'  | awk '{print $1"=\"Sym:"$2"|Orpha:"$3"\""}' > ${tables}/OrphanetTable
cat d_gnomad_constraint | awk '{print $1"=\"Sym:"$2"|oe_mis:"$3"|pLI:"$4"|pNull:"$5"|pRec:"$6"|oe_lof:"$7"|ensg_gn:"$8"|ensg_gi:"$9"\""}' > ${tables}/GnomadConstraintTable
cat a_gene_info_ensembl | awk '{if (length($5)==0) print $0 "\t" "NA"; else print $0}' | tr "&" "\t" | awk '{print $1"=\"Sym:"$2"|Syn:"$3"-|Desc:"$5"|ensg_gi:"$7"\""}' > ${tables}/GeneDescriptionTable
cat b_genemap2_has_entrezid.tsv | awk '{print $1"=\"Mim:"$2"|Sym:"$3"|OmimG:"$4"|OmimEnsg:"$6"|OmimPhen:"$7"\""}' > ${tables}/Genemap2Table
cat mim2gene.tsv | awk '{print $2"=\"Mim:"$1"|Sym:"$3"|OmimEnsg:"$4"\""}' > ${tables}/Mim2geneTable

## clean 
## comment the following line if you want to keep intermediate data
rm a* b* c* d*
cd ../

