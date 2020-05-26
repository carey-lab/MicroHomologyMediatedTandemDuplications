#!/usr/bin/env perl
# given two or more .counts.tsv files (last two columns are MTDs and MTCs) 
#    count the number of times each MTD and MTC occurs across all files
#    LBC May 2020

use warnings; 
use strict; 
use Data::Dumper ; 

my $EMSG = "$0 catchsig1.counts.tsv catchsig2.counts.tsv catchsig3.counts.tsv ...\n";
die $EMSG if ($#ARGV < 1); # not enough files

my @FH ; 
my $fI = 0;
foreach my $fn (@ARGV) {
	die "file $fn not found\n" unless (-e $fn);
	open( $FH[$fI] , $fn );
	$fI++;
}

while ( ! eof($FH[0]) ) {
	my $l = readline($FH[0]);
	my @l = split($l);
	print "@l\n";
	#	print "$l[0]\t$l[1]\t$l[2]\t$l[3]\t";
	#foreach my $fn (@ARGV) {
}

