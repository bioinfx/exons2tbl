#!/usr/local/bin/perl
use strict;
use warnings;

open GENEPOS, $ARGV[0] || die $!;
open PRODUCT, $ARGV[1] or die $!;

my %hashproduct;
while(<PRODUCT>){
        chomp $_;
        my @product=split("\t");
        $hashproduct{$product[0]}=$product[1];
        }

while(<GENEPOS>){
        chomp $_;
                my ($name, $pos) = split("\t");
                my @posarray = split(/,/,$pos);
                my $exonnum=@posarray/2;
#               print "$name\t$exonnum\n";
####    gene position
                print "$posarray[0]\t$posarray[@posarray-1]\tgene\n";
                print "\t\t\tgene\t$name\n";
####    CDS
                my $exon;
                for($exon=1; $exon<=$exonnum; $exon++){
                        if($exon == 1){print "$posarray[0]\t$posarray[1]\tCDS\n";}
                        else{print "$posarray[$exon*2-2]\t$posarray[$exon*2-1]\n";}
                }
####    product info
                print "\t\t\tproduct\t$hashproduct{$name}\n";
                print "\t\t\tprotein_id\tlcl| G-$name\n";
####    exon info
                for($exon=1; $exon<=$exonnum; $exon++){
                        print "$posarray[$exon*2-2]\t$posarray[$exon*2-1]\texon\n";
                        print "\t\t\tnumber\t$exon\n";
                }
####    intron info
                my $intron; my $intronstart; my $intronend;
                for($intron=1; $intron<$exonnum; $intron++){
                        if($posarray[0] < $posarray[1]){
                                $intronstart=$posarray[$intron*2-1]+1; # sense strand gene
                                $intronend=$posarray[$intron*2]-1;
                        }
                        else{
                                $intronstart=$posarray[$intron*2-1]-1;  # antisense strand gene
                                $intronend=$posarray[$intron*2]+1;
                        }
                        print "$intronstart\t$intronend\tintron\n";
                        print "\t\t\tnumber\t$intron\n";
                }
#               print "$name\n"; print "$posarray[0]\t$posarray[@posarray-1]\n"; exit;
        }

close GENEPOS;
close PRODUCT;
