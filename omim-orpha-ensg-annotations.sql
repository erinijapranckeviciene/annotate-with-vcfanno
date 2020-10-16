DROP TABLE IF EXISTS ginfoomim;
CREATE TABLE ginfoomim AS SELECT g.GeneID, 
  g.Symbol, 
  g.description, 
  o.MIM_number,
  o.Ensembl_Gene_ID, 
  o.Gene_Symbols, 
  o.Approved_Symbol, 
  o.Phenotypes  
FROM geneinfo g LEFT JOIN genemap2 o ON g.GeneID=o.Entrez_Gene_ID;

DROP TABLE IF EXISTS ginfoomimorpha;
CREATE TABLE ginfoomimorpha AS SELECT g.*, o.* 
 FROM ginfoomim g LEFT JOIN orpha o ON g.Symbol=o.Gene_symbol;

.mode tabs
.headers on
.out orpha_to_curate.csv
SELECT * FROM orpha WHERE Gene_symbol NOT IN (SELECT Symbol FROM geneinfo);
.out

DROP TABLE IF EXISTS ginfoomimorphaens;
CREATE TABLE ginfoomimorphaens AS SELECT 
  g.*, e.* FROM ginfoomimorpha g LEFT JOIN gene2ensembl e ON g.GeneID=e.GeneID;

DROP TABLE IF EXISTS annosource;
CREATE TABLE annosource AS SELECT
  rowid AS ID,
  GeneID, 
  Symbol, 
  description, 
  MIM_number, 
  Phenotypes, 
  Disorder_name as Orpha_disorder, 
  Disorder_OrphaCode as Orpha_number, 
  Association_type as Orpha_association, 
  Ensembl_gene_identifier as ENSG, 
  RNA_nucleotide_accession as NM, 
  SUBSTR(Ensembl_rna_identifier,1,15) as ENST 
FROM ginfoomimorphaens;

DROP TABLE IF EXISTS annosourceentrez;
CREATE TABLE annosourceentrez AS SELECT 
 ID, 
 GeneID as gid, 
 Symbol, 
 description, 
 MIM_number, 
 Phenotypes, 
 Orpha_disorder, 
 Orpha_number, 
 Orpha_association, 
 ENSG as ENSG,
 GeneID as ENTREZ,
 NM, 
 ENST 
FROM annosource;

DROP TABLE IF EXISTS annosourceensg;
CREATE TABLE annosourceensg AS SELECT 
 ID, 
 ENSG as gid, 
 Symbol, 
 description, 
 MIM_number, 
 Phenotypes, 
 Orpha_disorder, 
 Orpha_number, 
 Orpha_association, 
 ENSG as ENSG,
 GeneID as ENTREZ,
 NM, 
 ENST 
FROM annosource;

DROP TABLE IF EXISTS annottmp;
.mode tabs
.headers on
.out tmp/annottmp
SELECT * FROM annosourceentrez;
SELECT * FROM annosourceensg;

.out
.import tmp/annottmp annottmp
 
.mode tabs
.headers on
.out omim-orpha-ensg-annotations.csv
SELECT rowid as NUM, 
       gid as GeneID, 
       Symbol, 
       description as Description, 
       MIM_Number as MIM, 
       Phenotypes as OMIM_Phenotypes, 
       Orpha_disorder, 
       Orpha_number, 
       Orpha_association, 
       ENSG,
       ENTREZ,
       NM, 
       ENST 
FROM annottmp WHERE NOT INSTR(Symbol,"Sym") ORDER BY ID;

.exit
