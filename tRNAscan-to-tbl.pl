#!/usr/local/bin/perl
use warnings;
use strict;

# usage: perl tRNAscan-to-tbl.pl tRNAscan-SE_prediction.txt anticodon-database.tsv tRNA_result.tbl
# fengleiluck@gmail.com 2018-06-07

####  example of tRNAscan-SE_prediction.txt ####
#Please note (due to the simulation of a circular input sequence by GeSeq)
#tRNA predictions outside the range of the submitted sequence might be displayed.
#For details please refer to our Documentation.
#
#Sequence                tRNA    Bounds  tRNA    Anti    Intron Bounds   Inf           
#Name            tRNA #  Begin   End     Type    Codon   Begin   End     Score   Note
#--------        ------  -----   ------  ----    -----   -----   ----    ------  ------
#DNA           1       11407   11478   Arg     TCT     0       0       63.8
#DNA           2       30019   30089   Cys     GCA     0       0       64.7
#DNA           3       34347   34418   Thr     GGT     0       0       65.9

#### example of anticodon-database.tsv ####
#DNA     Codon   Anticodon       AminoAcid       Symbol  Fullname
#TTT     UUU     AAA     Phe     F       Phenylalanine
#TTC     UUC     GAA     Phe     F       Phenylalanine
#TTA     UUA     UAA     Leu     L       Leucine
#TTG     UUG     CAA     Leu     L       Leucine
#CTT     CUU     AAG     Leu     L       Leucine
#CTC     CUC     GAG     Leu     L       Leucine
#CTA     CUA     UAG     Leu     L       Leucine
#CTG     CUG     CAG     Leu     L       Leucine
#ATT     AUU     AAU     Ile     I       Isoleucine
#ATC     AUC     GAU     Ile     I       Isoleucine
#ATA     AUA     UAU     Ile     I       Isoleucine
#ATG     AUG     CAU     Met     M       Methionine


open TRNASCAN, $ARGV[0] or die $!;
open PRODUCTS, $ARGV[1] or die $!;
open TBL, ">$ARGV[2]" or die $!;

select TBL;

my $line = 0;
while(<TRNASCAN>){
	chomp $_;
	$line++; 
	if($line<8){next;}
	my @array=split(/\s+/);
#	print "$array[4]\t$array[2]\t$array[3]\n";

	my $AntiCodon = $array[5];
	$AntiCodon =~ tr/ACGTacgt/ACGUacgu/;
	my $Codon = reverse( $array[5] );
	$Codon =~ tr/ACGTacgt/UGCAugca/;
	my $AAsymbol;
#	print "$AntiCodon\n";
	while(<PRODUCTS>){
		chomp $_;
		my @codonarray=split(/\t/);
#		print "$codonarray[2]\n";
		if($codonarray[2] eq $AntiCodon){ $AAsymbol = $codonarray[4]; last; }
	}
	close PRODUCTS;
	open PRODUCTS, $ARGV[1] or die $!;
#	my $AAsymbol = &anticodonToSymbol($AntiCodon);
	my @name=split(//, $array[4]);
#	print "trn$name[0]-$AntiCodon\ttRNA-$array[4]\n";

	if($array[6] == 0){
		print "$array[2]\t$array[3]\tgene\n";
		print "\t\t\tgene\ttrn$AAsymbol-$AntiCodon\n";
		print "$array[2]\t$array[3]\ttRNA\n";
		print "\t\t\tnote\tanticodon:$AntiCodon\n";
#		print "\t\t\tprotein_id\tlcl| G-trn$AAsymbol\n";
		print "\t\t\tproduct\ttRNA-$array[4]\n";
	}

	if($array[6] > 0 ){
		my $exon1_end; my $exon2_start;
		if($array[6] < $array[7]){ $exon1_end=$array[6]-1;  $exon2_start=$array[7]+1;}
		if($array[6] > $array[7]){ $exon1_end=$array[6]+1;  $exon2_start=$array[7]-1;}
		print "$array[2]\t$array[3]\tgene\n";
		print "\t\t\tgene\ttrn$AAsymbol-$AntiCodon\n";
		print "$array[2]\t$exon1_end\ttRNA\n";
		print "$exon2_start\t$array[3]\n";
		print "\t\t\tnote\tanticodon:$AntiCodon\n";
#		print "\t\t\tprotein_id\tlcl| G-trn$AAsymbol\n";
		print "\t\t\tproduct\ttRNA-$array[4]\n";
		print "$array[2]\t$exon1_end\texon\n";
		print "\t\t\tgene\ttrn$AAsymbol-$AntiCodon\n";
		print "\t\t\tnumber 1\n";
		print "$exon2_start\t$array[3]\texon\n";
		print "\t\t\tgene\ttrn$AAsymbol-$AntiCodon\n";
		print "\t\t\tnumber 2\n";
		print "$array[6]\t$array[7]\tintron\n";
		print "\t\t\tgene\ttrn$AAsymbol-$AntiCodon\n";
		print "\t\t\tnumber 1\n";
	}


}


sub anticodonToSymbol {	
	my $AntiCodon = $_;
	while(<PRODUCTS>){
		chomp $_;
		my @arr=split(/\t/);
		if($arr[2] eq $AntiCodon){ return $arr[4]; }
	}
	close PRODUCTS;
	open PRODUCTS, $ARGV[1] or die $!;
}

close PRODUCTS;
close TRNASCAN;
close TBL;
