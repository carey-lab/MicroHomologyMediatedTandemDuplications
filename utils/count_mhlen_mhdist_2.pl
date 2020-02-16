#!/usr/bin/perl
# reads in multiple files created by count_mhlen_mhdist.pl
# sums them and outputs a sigle file
#  for the human genome (multiple chromosomes) on the cluster (one chr per node)
#
#  January 2020 ; LBC
#
use warnings;
use strict; 

# input file format:
# MHlen	InterMHDistance	N_MHPairs
# 10	100	204
# 9	93	670

my %h ; 
while(<>){
	chomp;
	next if m/^MHlen/ ; 
	my @l = split ;
	$h{"$l[0]\t$l[1]"} += $l[2] ;
}

print "MHlen\tInterMHDistance\tN_MHPairs\n" ;
foreach my $ID (sort keys %h){
	print "$ID\t$h{$ID}\n";
}

