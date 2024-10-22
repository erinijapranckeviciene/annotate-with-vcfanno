-- Lua functions to return MIM, OMIM, ORPHANET, Gene description, Gnomad OE Lof  annotations
-- Identifier mappings are based on the gene NCBI EntrezID
-- Annotation tables generated by create-annotation-tables.sh

mimFile  = "annotation-tables/Mim2geneTable"
genemap2File  = "annotation-tables/Genemap2Table"
geneinfoFile      = "annotation-tables/GeneDescriptionTable"
gnomadlofFile     = "annotation-tables/GnomadConstraintTable"
orphanetFile      = "annotation-tables/OrphanetTable"

-- ##################### FUNCTION TO READ FILE TO  TABLES ##########################################
-- Function to read external data as a string of type "key1=value1;...; keyX=valueX"
-- in which key ia alphanumeric type but value can be anything text " 00001=any %^&*: type of symbol; 88786=Gene : somegene_and_*(@#$%^) and anything ; ... "
function TableFromFile(File, newT)
for line in io.lines(File) do
  for key, value in string.gmatch(line,"(%w+)=(.+)") do
--    print(key)
--    print(value)
    newT[key]=value
  end
end  
end

-- ##################### INITIALIZE GLOBAL TABLES ####################################################
-- global table to keep entrezid to mim number mapping
mimTable={}
--print(mimFile)
TableFromFile(mimFile, mimTable)

-- global table for OMIM annotations
genemap2Table={}
--print(genemap2File)
TableFromFile(genemap2File,genemap2Table)

-- global table for gene name annotations from gene_info
geneinfoTable={}
--print(geneinfoFile)
TableFromFile(geneinfoFile,geneinfoTable)

-- global table for ORPHANET annotations
orphanetTable={}
--print(orphanetFile)
TableFromFile(orphanetFile,orphanetTable)

-- global table for ORPHANET annotations
gnomadlofTable={}
--print(gnomadlofFile)
TableFromFile(gnomadlofFile,gnomadlofTable)

-- test the contents of the table
-- TableToTest=mimTable
--for k,v in pairs(TableToTest) do print(k,v) end

-- ##################### DEFINE FUNCTIONS ################################################################ 
-- Function returns annotation associated with gene EntrezID from the table 
function EntrezIDToAnnotation (entrezid, annotationTable)
  return annotationTable[tostring(entrezid)]
end

-- Test script 
--id = {10, 100, 1000, 9991, 9992, 9994, 9997}
--print(" Test annotations")
--for i=1,7 do
--  entrez=id[i]
--  print("For EntrezID")
--  print(entrez) 
--  print("mim to gene")
--  print( EntrezIDToAnnotation(entrez, mimTable) )
--  print("gememap2")
--  print( EntrezIDToAnnotation(entrez, genemap2Table) )
--  print("gene description")
--  print( EntrezIDToAnnotation(entrez, geneinfoTable) )
--  print("orphanet")
--  print( EntrezIDToAnnotation(entrez, orphanetTable) )
--  print("gnomad lof")
--  print( EntrezIDToAnnotation(entrez, gnomadlofTable) )
--print("\n\n")
--end

