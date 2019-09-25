#!/usr/bin/perl
use strict;

# script to parse biotype and IDs info from fa files

open (inFile, $ARGV[0]) || die "Can't open in $!\n";
open (oFile_ID, ">$ARGV[1]") || die "Can't open IDs $!\n";
open (oFile_Bio, ">$ARGV[2]") || die "Can't open biotypes $!\n";

print oFile_ID "TransID\tGeneID\n";
print oFile_Bio "TransID\tBiotype\tGeneName\n";



while (my $line = <inFile>){
	next unless ($line =~ /^\>([\w-.]*).*gene:(\w+)/);

	print oFile_ID "$1\t$2\n";

	my $transcript = $1;

	$line =~ /gene_biotype:(\w+)\s/;
	my $biotype = $1;

	$line =~ /gene_symbol:(\w+)\s/;
	my $gene = $1;

	print oFile_Bio "$transcript\t$biotype\t$gene\n";

}

close inFile;
close oFile_ID;
close oFile_Bio;

__DATA__

>ENSMUST00000179664.1 cdna chromosome:GRCm38:14:54113468:54113478:1 gene:ENSMUSG00000096749.2 gene_biotype:TR_D_gene transcript_biotype:lncRNA gene_symbol:Trdd1 description:T cell receptor delta diversity 1 [Source:MGI Symbol;Acc:MG
