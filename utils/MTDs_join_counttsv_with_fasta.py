#!/usr/bin/env python3

# read in a counts.tsv file from catch_signatures.awk
# get MHseq, and GC content of MHseq and inter-MH distance
# also transform into bed-like format for the first five columns
#  chr MHPstart MHPend name(MHLen) score(NDupReads)
#
#
# (1) read in counts.tsv file
# (2) read fasta file
#
# LBC March 2020


import pyfastx
import argparse
import os
import pandas
import re
import string
import numpy as np
import time


def GC(x):
    return round( len(re.findall("[GCgc]",x)) / len(x)  * 100 )

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', help='FASTA input filename')
parser.add_argument('--countstsv', help='counts.tsv file produced by catch_signatures.awk')

args = parser.parse_args()

# for development
os.chdir('/Users/lcarey/Downloads')
args.fasta = '/Users/lcarey/CareyLab/ExternalData/PomBase/pombe_n.fasta'
args.countstsv = 'test.tsv'
args.countstsv = '/Users/lcarey/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/SRR7817502_n.sign.count_500bp.tsv'

# load input files
print(args.fasta, ' is ', os.path.getsize(args.fasta), ' bytes')
fa = pyfastx.Fasta(args.fasta)

print(args.countstsv, ' is ', os.path.getsize(args.countstsv), ' bytes')
tic = time.time()
cnts_tsv_df = pandas.read_csv(args.countstsv, sep='\t', names=["chr", "s1", "e1", "s2", "e2", "NDups", "NCol"])
toc = time.time()
print(f"... loaded. Adding Sequence data took {toc - tic:0.2f} seconds")

cnts_tsv_df["s01"] = cnts_tsv_df.s1 - 1
cnts_tsv_df["e01"] = cnts_tsv_df.e1 - 1
cnts_tsv_df["s02"] = cnts_tsv_df.s2 - 1
cnts_tsv_df["e02"] = cnts_tsv_df.e2 - 1

# Assign MHseq, MHlen, & GC
cnts_tsv_df["MH_length"] = cnts_tsv_df.e1 - cnts_tsv_df.s1 + 1
MH = ["" for x in range(len(cnts_tsv_df))]
interMH = ["" for x in range(len(cnts_tsv_df))]

tic = time.time()
for I in range(len(cnts_tsv_df)):
    MH[I] = fa.fetch( cnts_tsv_df.chr[I] , (cnts_tsv_df.s1[I].item() , cnts_tsv_df.e1[I].item()) )
    interMH[I] = fa.fetch( cnts_tsv_df.chr[I] , (cnts_tsv_df.e1[I].item()+1 , cnts_tsv_df.s02[I].item()) )

toc = time.time()
cnts_tsv_df["MH_Sequence"] = MH
cnts_tsv_df["InterMH_Sequence"] = interMH
print(f"... done. Extracting sequences took {(toc - tic)/60:0.2f} minutes. Adding GC & modulo 3")

cnts_tsv_df["MH_Sequence_GC"] = cnts_tsv_df["MH_Sequence"].apply(GC)
cnts_tsv_df["InterMH_Sequence_GC"] = cnts_tsv_df["InterMH_Sequence"].apply(GC)
cnts_tsv_df["DuplicatedRegion_GC"] = cnts_tsv_df["MH_Sequence"].add(cnts_tsv_df["InterMH_Sequence"]).apply(GC)
cnts_tsv_df["CompleteMHP_GC"] = cnts_tsv_df["MH_Sequence"].add( cnts_tsv_df["MH_Sequence"].add(cnts_tsv_df["InterMH_Sequence"]) ).apply(GC)
cnts_tsv_df["Duplication_Creates_modulo3_insertion"] = (cnts_tsv_df["MH_Sequence"].apply(len)+cnts_tsv_df["InterMH_Sequence"].apply(len)) % 3==0


# save
tic = time.time()
cnts_tsv_df.to_pickle( "/Users/lcarey/df.pickle" )
toc = time.time()
print(f"Saving as pickle took {(toc - tic)/60:0.2f} minutes.")

tic = time.time()
cnts_tsv_df.to_hdf( "/Users/lcarey/df.h5" , key='df', mode='w')
toc = time.time()
print(f"Saving as h5 took {(toc - tic)/60:0.2f} minutes.")

tic = time.time()
cnts_tsv_df.to_feather( "/Users/lcarey/df.feather" , key='df', mode='w')
toc = time.time()
print(f"Saving as feather took {(toc - tic)/60:0.2f} minutes.")