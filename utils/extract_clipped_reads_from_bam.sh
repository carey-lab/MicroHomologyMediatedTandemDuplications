#!/bin/bash
# extract soft-clipped reads
# LBC May 2020
# 
# example: 
#  extract_clipped_reads_from_bam.sh test.cram  > test.bam
{
samtools view -H $1
samtools view $1 | cut -f 1-10  \
	| grep  '[[:space:]|[:digit:]M][[:digit:]]S' \
	| awk '{print NR "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t"  $9 "\t" $10 "\t" "*" }' 
} | samtools view -h -O BAM - -o - 
