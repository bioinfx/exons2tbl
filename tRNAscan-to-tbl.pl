#!/home/fenglei/local/bin/perl
use warnings;
use strict;

# usage: perl tRNAscan-to-tbl.pl tRNAscan-SE_prediction.txt tRNA-product.tsv tRNA_result.tbl
# fengleiluck@gmail.com 2018-06-07

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
