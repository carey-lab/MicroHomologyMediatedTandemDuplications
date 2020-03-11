#!/usr/bin/env python3
#
# given as input a .vcf file and an fasta(.gz) file
# output fasta sequence with all REF sequences and all ALT seqs, plus flanking seq on either side
import pyfaidx
import vcf
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument( '-f' , '--fasta' , required=True , help='fasta file')
parser.add_argument( '-v' , '--vcf' , required=True ,  help='.vcf file')
parser.add_argument( '-nt' , '--nt_flank' , default=250 , type=int , help='get NT on either side of the insertion')
parser.add_argument( '--MIN_LENGTH_INSERTION' , default=5 , type=int , help='min insertion length to process')

args = parser.parse_args()

# ftp://ftp.ensemblgenomes.org/pub/fungi/release-46/variation/vcf/schizosaccharomyces_pombe/schizosaccharomyces_pombe_incl_consequences.vcf.gz.csi
# ftp://ftp.ensemblgenomes.org/pub/fungi/release-46/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz
WORKDIR = '/Users/lcarey/Downloads/test/'
#args.fasta = WORKDIR + 'Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz'
#args.vcf   = WORKDIR + 'schizosaccharomyces_pombe_indels.vcf'
#nt_flank = 50
#cons = pyfaidx.FastaVariant( fasta_file_name , vcf_file_name , het=True , hom=True )


fa = pyfaidx.Faidx(args.fasta)

vcf_reader = vcf.Reader(open(args.vcf, 'r'))


for record in vcf_reader:
    print_me_flag = False
    for r in record.ALT:
        if len(r.sequence) >= args.MIN_LENGTH_INSERTION :
            print_me_flag = True
    if print_me_flag:
        left_seq  = fa.fetch(record.CHROM , record.start - args.nt_flank , record.start)
        right_seq = fa.fetch(record.CHROM , record.end + 1 , record.end + args.nt_flank + 1)
        print( '>' + record.CHROM + '_' + str(record.POS)  + '_' + record.REF + '_REF' )
        print( left_seq.seq  + record.REF +  right_seq.seq )
        for r in record.ALT:
            if len(r.sequence) >= args.MIN_LENGTH_INSERTION :
                print( '>' + record.CHROM + '_' + str(record.POS)  + '_' + record.REF + '_ALT_' + r.sequence[0:args.MIN_LENGTH_INSERTION] + '_' + str(len(r.sequence)) )
                print( left_seq.seq  + r.sequence +  right_seq.seq )
        


