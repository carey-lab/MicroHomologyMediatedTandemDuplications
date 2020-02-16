#!/usr/bin/perl
# give an .sign.count.tsv file
#  count the number of MHpairs for each MHlength and inter-MH distance
#  distance includes MHpairs
#
# chr6	60002	60005	60024	60027	0	0
use warnings; 
use strict; 

my %h  ;
while(<>){
	my @l = split;
	my $MHLen = $l[2] - $l[1] + 1 ; 
	my $InterMHDistance = $l[4] - $l[1] + 1 ; 
	my $ID = "$MHLen\t$InterMHDistance" ; 
	$h{$ID}++;
}


print "MHlen\tInterMHDistance\tN_MHPairs\n" ;
foreach my $ID (sort keys %h){
	print "$ID\t$h{$ID}\n";
}
