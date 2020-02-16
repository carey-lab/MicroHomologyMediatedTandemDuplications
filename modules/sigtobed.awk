#! /usr/local/bin/gawk -f

{
	print $1, $2 | "sort -k1,1 -k2n,2 | uniq"
	print $1, $3 | "sort -k1,1 -k2n,2 | uniq"
	print $1, $4 | "sort -k1,1 -k2n,2 | uniq"
	print $1, $5 | "sort -k1,1 -k2n,2 | uniq"
}