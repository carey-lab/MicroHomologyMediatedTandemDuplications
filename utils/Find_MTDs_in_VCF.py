#!/usr/bin/env python3
#
# given as input a .vcf file and an fasta(.gz) file
# output fasta sequence with all REF sequences and all ALT seqs, plus flanking seq on either side
import pyfaidx
import vcf
import argparse
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


parser = argparse.ArgumentParser()
parser.add_argument( '-f' , '--fasta' , required=True , help='fasta file')
parser.add_argument( '-v' , '--vcf' , required=True ,  help='.vcf file')
parser.add_argument( '-nt' , '--nt_flank' , default=250 , type=int , help='get NT on either side of the insertion')
parser.add_argument( '--MIN_LENGTH_INSERTION' , default=5 , type=int , help='min insertion length to process')
parser.add_argument( '-o' , '--output_basename' , required=True , help='save output with <OUTPUT_BASENAME>.fasta')
parser.add_argument( '--no-run_trf' , action='store_false' , dest='run_trf', help='optional: run TRF on output (default = TRUE)')
parser.add_argument( '--run_trf' , default=True  , dest='run_trf', action='store_true' , help='optional: run TRF on output (default = TRUE)')
parser.add_argument( '--trf_args' , default='  1 5 1000 80 10 20 100  -f -d -l 2 -ngs -h  '  , help='arguments to pass to TRF')
parser.add_argument( '--extract_awk_cmd' , default='gawk -f ~/Develop/MicroHomologyMediatedTandemDuplications/reference_genome_TRs/results/extract.awk'  , help='parse TRF output into tab file for MHP finding')
parser.add_argument( '--find_MHPs_cmd' , default='Rscript ~/Develop/MicroHomologyMediatedTandemDuplications/reference_genome_TRs/code/TR_MH_build_results_table.R'  , help='parse processed TRF tab file & find MHPs')
args = parser.parse_args()

# ftp://ftp.ensemblgenomes.org/pub/fungi/release-46/variation/vcf/schizosaccharomyces_pombe/schizosaccharomyces_pombe_incl_consequences.vcf.gz.csi
# ftp://ftp.ensemblgenomes.org/pub/fungi/release-46/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz
#  WORKDIR = '/Users/lcarey/Downloads/test/'
#  args.fasta = WORKDIR + 'Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz'
#  args.vcf   = WORKDIR + 'schizosaccharomyces_pombe_indels.vcf'
#nt_flank = 50
#cons = pyfaidx.FastaVariant( fasta_file_name , vcf_file_name , het=True , hom=True )

# intput and output files
fasta_output_filename = args.output_basename + '.fasta'
seqs_to_write = [] # empty list of seqs we'll save
fa = pyfaidx.Faidx(args.fasta)
vcf_reader = vcf.Reader(open(args.vcf, 'r'))

# for each record in the VCF file with an insertion length > MIN_LENGTH_INSERTION
#   generate the REF and each ALT seq
#   save these to a fasta file
for record in vcf_reader:
    print_me_flag = False
    for r in record.ALT:
        if len(r.sequence) >= args.MIN_LENGTH_INSERTION :
            print_me_flag = True
    if print_me_flag:
        try:
            left_seq  = fa.fetch(record.CHROM , record.start - args.nt_flank , record.start)
            right_seq = fa.fetch(record.CHROM , record.end + 1 , record.end + args.nt_flank + 1)
            seqrec = SeqRecord( Seq(left_seq.seq  + record.REF +  right_seq.seq) , 
                id=record.CHROM + '_' + str(record.POS)  + '_' + record.REF + '_REF' , description='')
            seqs_to_write.append(seqrec)
            for r in record.ALT:
                if len(r.sequence) >= args.MIN_LENGTH_INSERTION :
                    seqname = record.CHROM + '_' + str(record.POS)  + '_' + record.REF + '_ALT_' + r.sequence[0:args.MIN_LENGTH_INSERTION] + '_' + str(len(r.sequence)) 
                    seqrec = SeqRecord( Seq(left_seq.seq  + r.sequence +  right_seq.seq) , id=seqname , description='')
                    seqs_to_write.append(seqrec)
        except:
            pass

SeqIO.write( seqs_to_write , fasta_output_filename , 'fasta' ) # write all sequences at once to the fasta file

# optionally, run trf on the output
if (args.run_trf):
    
    # run TRF
    cmd = 'trf ' + fasta_output_filename + args.trf_args + ' > ' + args.output_basename + '.trf.txt'
    print( cmd )
    subprocess.run(  cmd  ,  shell=True , check=True)

    # process the output of TRF
    cmd = args.extract_awk_cmd + ' ' + args.output_basename + '.trf.txt' + ' > ' + args.output_basename + '.trf.out'
    print( '\n' , cmd )
    subprocess.run(  cmd  ,  shell=True , check=True)
 
    # R script to find MHPs
    cmd = args.find_MHPs_cmd + ' --fasta_file ' + fasta_output_filename \
        + '  --trf_output_file ' + args.output_basename + '.trf.out' \
            + ' -o ' + args.output_basename + '.txt'
    print( '\n' , cmd )
    subprocess.run(  cmd  ,  shell=True , check=True)