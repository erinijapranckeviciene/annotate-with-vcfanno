# annotate-with-vcfanno

The [vcfanno](https://github.com/brentp/vcfanno) annotates vcf files with another bed and vcf files. Sometimes there is a need to annotate with information that is associated with gene or other identifier such as MIM number, OMIM disease, ORPHANETdisorder, Gene description or Gnomad OE Lof annotations.
The lua and conf files is this example attempts to define annotations that are associted with gene EntrezID and gene symbol. The annotation source table examples in *annotation-tables* contain some examples that work with Lua functions defined in cheop.lua. 

Here is a tiny self-sufficient example that demonstrates annotation by identifiers only with three variants in vcf file. 

## How to try example
Install [vcfanno](https://github.com/brentp/vcfanno). Then

```

git clone https://github.com/erinijapranckeviciene/annotate-with-vcfanno.git
cd annotate-with-vcfanno
vcfanno -lua cheop.lua cheop.conf NA12878.example.vcf  

```
You can add -p n option to the call to vcfanno to use n cores. The output will be vcf annotated with OMIM, ORPHANET and Gnomad OE Lof annotations. The configuration cheop.conf also contains instructions to compute gnomaAD_POPMAX annotation.
Instructions below show how to prepare and use full annotation data tables. 

##Instructions and data sources to create annotation data tables

###Annotation data sources

The identifier mapping is based on NCBI gene EntrezID and gene symbol. MIM and OMIM annotations downloaded from https://www.omim.org/downloads/. The gene2ensembl.gz identifiers and  gene_info.gz are from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/. The ORPHANET rare diseases with their associated genes are from http://www.orphadata.org/data/xml/en_product6.xml. The gnomad lof constraint data downloaded from 
https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

###Annotation data preprocessing steps

OMIM mim2gene:
---------------
`$cat mim2gene-2019.txt| grep gene | grep -v ^# | awk '{print $3"="$1}' | grep -v ^=  > mimTable`

OMIM genemap2.txt:
-------------------
`$cat genemap2.txt | tr " " _ | grep -v ^# | tr "\t" = | sed 's/==/=NA=/g' | sed 's/==/=NA=/g' | tr "=" "\t" | cut -f6,8,10,11,13 | tr -d " " | awk '{print $1"= \"Gene: "$2" | EntrezIDom: "$3" | Ensg: "$4" | Phenotype: "$5"\""}' > genemap2Table`

NCBI GENE gene_info.gz:
------------------------
`$zcat gene_info.gz | grep ^9606 | tr " " _ | cut -f2,3,9 | awk '{print $1"=\""$2": "$3"\""}' > gene_infoTable`

NCBI GENE gene2ensembl.gz:
---------------------------
`$zcat gene2ensembl.gz | grep ^9606 | tr " " _ | cut -f2,3 | awk '{print $1"="$2}' |  tr "\n" ";"  | sed 's/;$//' > gene2ensemblTable`

ORPHANET: 
----------
Since ORPHANET data is in xml format we need to remove irrelevant lines and convert to UTF-8 charset that XML parser accepts, and parse XML file into the table. The script orphanetXMLparse.pl parses into rows that contain gene disease association pairs for each gene. The orphanetXMLparse_unique.pl parses into lines each containing a unique gene and a list of associated diseases with their ORPHA number.   
```

$tail -n +10 en_product6.xml | head -n -1 | dos2unix -iso | iconv -f ISO-8859-15 -t UTF-8 > test_en_product6.xml
$perl orphanetXMLparse.pl test_en_product6.xml > orphanetTable          
$perl orphanetXMLparse_unique.pl test_en_product6.xml > orphanetTable2

```

gnomAD lof by gene:
--------------------
```

$mv gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz gnomad.v2.1.1.lof_metrics.by_gene.txt.gz
$gunzip gnomad.v2.1.1.lof_metrics.by_gene.txt.gz
$cut -f1,5,24,64 | grep -v gene | awk '{print $1"=\"gnomad_oe_mis: " $2" | gnomad_oe_lof: " $3 " | "$4"\""}' > gnomadlofTable

```
