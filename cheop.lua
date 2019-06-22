-- External files
mimFile  = "annotation-tables/mimTable"
omimgenemap2File  = "annotation-tables/genemap2Table"
geneinfoFile      = "annotation-tables/gene_infoTable"
gene2ensemblFile  = "annotation-tables/gene2ensemblTable"
gnomadlofFile     = "annotation-tables/gnomadlofTable"
orphanetFile      = "annotation-tables/orphanetTable2"


-- ##################### FUNCTIONS TO READ FILE TO  TABLES ##########################################
-- Function to read external data as a string of type "key1=value1;...; keyX=valueX"
-- in which key and value are only alpahnumeric types" 00001=990000; 88786=ENSG000001345; ... "
function ReadDataAsString(File, newT)
  local f=assert(io.open(File,"r"))
  local s=f:read("*all")
  for key, value in string.gmatch(s,"(%w+)=(%w+)") do
    newT[key]=value
  end
  f:close()
  --return newT
end

-- Function to read external data as a string of type "key1=value1;...; keyX=valueX"
-- in which key ia alphanumeric type but value can be anything text " 00001=any %^&*: type of symbol; 88786=Gene : somegene_and_*(@#$%^) and anything ; ... "
function ReadDataFromFile(File, newT)
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
entrezid2mimT={}
ReadDataFromFile(mimFile,entrezid2mimT)
--for k,v in pairs(entrezid2mimT) do print(k,v) end

-- global table to keep entrezid to ensembl ENSG mapping
gene2ensemblT={}
ReadDataAsString(gene2ensemblFile,gene2ensemblT)
--for k,v in pairs(gene2ensemblT) do print(k,v) end


-- global table for OMIM annotations
genemap2T={}
--print(omimgenemap2File)
ReadDataFromFile(omimgenemap2File,genemap2T)
--for k,v in pairs(genemap2T) do print(k,v) end

-- global table for gene name annotations from gene_info
geneinfoT={}
--print(geneinfoFile)
ReadDataFromFile(geneinfoFile,geneinfoT)
--for k,v in pairs(geneinfoT) do print(k,v) end

-- global table for ORPHANET annotations
orphanetT={}
--print(geneinfoFile)
ReadDataFromFile(orphanetFile,orphanetT)
--for k,v in pairs(orphanetT) do print(k,v) end

-- global table for ORPHANET annotations
gnomadlofT={}
--print(geneinfoFile)
ReadDataFromFile(gnomadlofFile,gnomadlofT)
--for k,v in pairs(gnomadlofT) do print(k,v) end

-- ##################### DEFINE FUNCTIONS ################################################################ 
-- Function returns MIM number associated with gene EntrezID
function EntrezIDToMim (entrezid)
  return entrezid2mimT[tostring(entrezid)]
end

-- Function returns OMIM gene and phenotype accociated with mim number that is associated with gene EntrezID
function OMIMp (entrezid)
  mim  = entrezid2mimT[tostring(entrezid)]
  -- consider changing the : into the = and the | into the ; 
  -- omimannot = string.gsub(string.gsub(genemap2T[tostring(mim)],': ','='), ' | ',';')
  omimannot = genemap2T[tostring(mim)]
  return omimannot
end

-- Function returns gene symbol and description associated with gene EntrezID
function EntrezIDToGeneDescription (entrezid)
  return geneinfoT[tostring(entrezid)]
end

-- Function returns Ensembl ID  associated with gene EntrezID
function EntrezIDToEnsembl (entrezid)
  return gene2ensemblT[tostring(entrezid)]
end

-- Need to return gene symbol and query ORPHANET and GNOMAD

-- Function returns ORPHANET query
function EntrezIDToOrpha (entrezid)
-- extract genesymbol from line
  line = geneinfoT[tostring(entrezid)]
  genesymbol = string.gsub(string.gsub(line,':.*',''),'"','')
  return orphanetT[tostring(genesymbol)]
end

-- Function returns Gnomad OE Lof  query
function EntrezIDToGnomadLof (entrezid)
-- extract genesymbol from line
  line = geneinfoT[tostring(entrezid)]
  genesymbol = string.gsub(string.gsub(line,':.*',''),'"','')
  return gnomadlofT[tostring(genesymbol)]
end

-- test query in lua
-- uncomment the following lines
-- id={1, 859}
--for i=1,2 do
--  entrez=id[i]
--  print("entrez")
--  print(entrez)
--  print("MIM")
--  print(EntrezIDToMim(entrez))
--  print("OMIM")
--  print(OMIMp(entrez))
--  print("Gene Description")
--  print(EntrezIDToGeneDescription(entrez))
--  print("Ensembl ID")
--  print(EntrezIDToEnsembl(entrez))
--  print("orpha")
--  print(EntrezIDToOrpha(entrez))
--  print("Gnomad OE Lof")
--  print(EntrezIDToGnomadLof(entrez))
-- end
