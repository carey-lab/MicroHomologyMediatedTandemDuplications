#!/bin/bash
# extract soft-clipped reads
# LBC May 2020
# 
# example: 
#  extract_clipped_reads_from_bam.sh test.cram  > test.bam
#     grep  '[[:space:]|[:digit:]M][[:digit:]]S' 
{
samtools view -H $1
samtools view $1 | cut -f 2-10  | \
	awk '{ if ($5 ~ /[SID]/) print NR "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t"  $8 "\t" $9 "\t" "*" } '
} | samtools view -h -O BAM - -o - 
