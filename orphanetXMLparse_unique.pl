#!/usr/bin/perl
# parse into gene = list of disorders

use strict; 
use warnings;

use XML::LibXML;

my $filename=shift;

my $disorderline;
my $genesymbol;

my %DisorderToGene;
my %GeneToDisorder;

## Parse disorders and their associated gene names into DisorderToGene
my $dom = XML::LibXML->load_xml(location => $filename);

foreach my $disorder ($dom->findnodes('/DisorderList/Disorder')){
    $disorderline=$disorder->findvalue("./Name")."; ORPHAID: ".$disorder->findvalue("./OrphaNumber");

    foreach my $gene ($disorder->findnodes("./DisorderGeneAssociationList/DisorderGeneAssociation/Gene")){
        $genesymbol=$gene->findvalue("./Symbol");        

    $DisorderToGene{$disorderline}=$genesymbol;
    }
}

## create the  gene = list of disorders hash
foreach my $disorder( keys(%DisorderToGene)){
    my $gene=$DisorderToGene{$disorder};
    if (exists($GeneToDisorder{$gene})){ 
        $GeneToDisorder{$gene}=$GeneToDisorder{$gene}." | ".$disorder;}
    else { $GeneToDisorder{$gene}=$disorder;}
}

## print the hash uniquely for genes
foreach my $gene (sort(keys(%GeneToDisorder) ) ){
    print $gene."=\"".$GeneToDisorder{$gene}."\"\n" }
