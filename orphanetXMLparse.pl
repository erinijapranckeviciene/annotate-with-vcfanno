#!/usr/bin/perl
use strict; 
use warnings;

use XML::LibXML;

my $filename=shift;

my $disorderline;
my $geneline;
my $genesymbol;

my $dom = XML::LibXML->load_xml(location => $filename);

foreach my $disorder ($dom->findnodes('/DisorderList/Disorder')){
    $disorderline="Disorder: ".$disorder->findvalue("./Name").";orphaID: ".$disorder->findvalue("./OrphaNumber").";";

    foreach my $gene ($disorder->findnodes("./DisorderGeneAssociationList/DisorderGeneAssociation/Gene")){
        $genesymbol=$gene->findvalue("./Symbol");        
        $geneline=$genesymbol.": ".$gene->findvalue("./Name")."|";

        foreach my $externalreference ($gene->findnodes("./ExternalReferenceList/ExternalReference")){
            $geneline=$geneline.$externalreference->findvalue("./Source").": ".$externalreference->findvalue("./Reference")."|"; 
        }
    print $genesymbol."=\"".$disorderline.$geneline,";\"\n"
    }
}
