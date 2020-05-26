#!/usr/bin/env perl
# given two or more .counts.tsv files (last two columns are MTDs and MTCs) 
#    count the number of times each MTD and MTC occurs across all files
#    LBC May 2020

use warnings; 
use strict; 
use Data::Dumper ; 

my $EMSG = "$0 catchsig1.counts.tsv catchsig2.counts.tsv catchsig3.counts.tsv ...\n";
die $EMSG if ($#ARGV < 1); # not enough files

# open all files, storing array of filehandles
my @FH ; 
my $fI = 0;
foreach my $fn (@ARGV) {
        die "file $fn not found\n" unless (-e $fn);
        open( $FH[$fI] , $fn );
        $fI++;
}

# read through the each file, counting the number of files with an MTD or MTC, and print per line
while (my @l = split( /\t/ ,  readline $FH[0] ) ){
	my $MTDs = 0 ;
	my $MTCs = 0 ; 
	print "$l[0]\t$l[1]\t$l[2]\t$l[3]\t$l[4]\t" ;
	$MTDs = 1 if ($l[5]>0) ;
	$MTCs = 1 if ($l[6]>0) ;
	for (my $fI = 1 ; $fI <= $#FH ; $fI++){
		my $line = readline $FH[$fI];
		my @line = split( /\t/ , $line );
		$MTDs++ if ($line[5]>0) ;
		$MTCs++ if ($line[6]>0) ;
	}
	print "$MTDs\t$MTCs\n";
}
