#!/usr/bin/perl
# Synopsis: print gene_info line for each gene symbol and its synonyms
# Use file with 
use strict; 
use warnings;

my $filename=shift;

open(DAT,"<$filename") or die "err $! \n";
while(<DAT>)
{ chop;
  my( $entrezid, $symbol, $synonyms, $geneinfo, $ensg) = split("\t",$_);
  #print $synonyms,"\n";
  print $entrezid, "\t", $symbol , "\t", defined $ensg ? $ensg : "NA", "\t", $geneinfo, "\n";

  foreach my $synonym( split("\\|",$synonyms) ){
    if ($synonym ne '-'){
      print $entrezid, "\t", $synonym, "\t" , defined $ensg ? $ensg : "NA", "\t",  $geneinfo, "\n";}
  } 
}
close(DAT);

