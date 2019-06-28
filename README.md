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

Suggestions and input are welcome!
