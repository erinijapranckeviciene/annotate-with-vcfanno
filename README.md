# annotate-with-vcfanno

The [vcfanno](https://github.com/brentp/vcfanno) annotates vcf files with another bed and vcf files. Sometimes there is a need to annotate with information that is associated with gene or other identifier such as MIM number, OMIM disease, ORPHANETdisorder, Gene description or Gnomad OE Lof annotations.
The lua and conf files is this example attempts to define annotations that are associated with gene EntrezID and gene symbol. The annotation source table examples in *annotation-tables* contain some examples that work with Lua functions defined in cheop.lua. 

This self-sufficient example demonstrates annotation by identifiers. 

## How to try this example
Install [vcfanno](https://github.com/brentp/vcfanno). Then execute:

```
git clone https://github.com/erinijapranckeviciene/annotate-with-vcfanno.git
cd annotate-with-vcfanno
vcfanno -lua cheop.lua cheop.conf example.vcf  
```

Add -p n option to use n cores. The  vcf annotated with EnrezID, OMIM, ORPHANET and Gnomad OE Lof annotations will print to the stdout. 
This example uses only Entrez ID 85358. 

To try full annotation tables use script *create-annotation-tables.sh*. The OMIM genemap2.txt should be present in annotation-tables folder for that script to create genemap2Table. 

## Query gemini db created from annotated vcf files

Annotated vcf files in [bcbio-nextgen](https://bcbio-nextgen.readthedocs.io/en/latest/contents/introduction.html) can be uploaded to 
[gemini](https://gemini.readthedocs.io/en/latest/index.html) sqlite3 database, in which INFO fields from vcf header become columns in variant table.

To create geminin db offline with [vcf2db.py](https://github.com/quinlan-lab/vcf2db) you will need vcf and ped files and follow the provided instructions. 
The created gemini db can be queried with sqlite3 SQL statements to extract required fields. 
The hypothetical family.db on which some queries are shown here was created as folows:

```
vcf2db.py family.vcf family.ped family.db
```

The *gemini-query-with-grouping.sh* and *gemini-query-with-grouping-platypus-first.sh* examples show some queries for family.db created from family.vcf file
annotated by identifiers using cond and lua files from this repository. The *gemini-query-answer-examples* folder contains several answered queries derived as:

```
sh gemini-query-with-grouping.sh family.db | head -n5 | datamash transpose | cat -n > gemini-query-answer-examples/query-answer1-example.csv
```


