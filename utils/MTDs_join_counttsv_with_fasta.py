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
parser.add_argument('--output_filename', help='<file name to save dataframe>.feather')

args = parser.parse_args()

## for development
#os.chdir('/Users/lcarey/Downloads')
#args.fasta = '/Users/lcarey/CareyLab/ExternalData/PomBase/pombe_n.fasta'
#args.countstsv = '/Users/lcarey/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/SRR7817502_n.sign.count_500bp.tsv'

# load input files
print(args.fasta, ' is ', os.path.getsize(args.fasta), ' bytes')
fa = pyfastx.Fasta(args.fasta)

print(args.countstsv, ' is ', os.path.getsize(args.countstsv), ' bytes')
tic = time.time()
cnts_tsv_df = pandas.read_csv(args.countstsv, sep='\t', names=["chr", "s1", "e1", "s2", "e2", "NDups", "NCol"])
toc = time.time()
print(f"... loaded. Adding Sequence data took {toc - tic:0.2f} seconds")

## smaller df for debugging 
#cnts_tsv_df = cnts_tsv_df.loc[ (cnts_tsv_df["chr"]=="I") ]
cnts_tsv_df = cnts_tsv_df.loc[ (cnts_tsv_df["chr"]=="I") | (cnts_tsv_df["chr"]=="II")  | (cnts_tsv_df["chr"]=="III")  ] 
#cnts_tsv_df = cnts_tsv_df.loc[ (cnts_tsv_df["s1"] < 1000000) | (cnts_tsv_df["NDups"] > 0) ] 
cnts_tsv_df = cnts_tsv_df.reset_index(drop=True) # rebuild row numbers in the shrunken dataframe

cnts_tsv_df["s01"] = cnts_tsv_df.s1 - 1
cnts_tsv_df["e01"] = cnts_tsv_df.e1 - 1
cnts_tsv_df["s02"] = cnts_tsv_df.s2 - 1
cnts_tsv_df["e02"] = cnts_tsv_df.e2 - 1

# put the position columns into the order needed to make a .bed file  : s01 & e2
cnts_tsv_df["MH_length"] = cnts_tsv_df.e1 - cnts_tsv_df.s1 + 1
cnts_tsv_df = cnts_tsv_df[['chr','s01','e2','MH_length','NDups','s1','s2','e1','e01','s02','e02']]




# Assign MHseq, MHlen, & GC
tic = time.time()

# way way way way faster to not work with dataframes
#   extract all info into lists, save into a list, and put the result into the df
MHseq = ["" for x in range(len(cnts_tsv_df))]
interMHseq = ["" for x in range(len(cnts_tsv_df))]

c=cnts_tsv_df['chr'].to_list()
s1=cnts_tsv_df['s1'].to_list()
e1=cnts_tsv_df['e1'].to_list()
s02=cnts_tsv_df['s02'].to_list()

for I in range(len(cnts_tsv_df)):
    MHseq[I] = fa.fetch( c[I] , (s1[I] , e1[I] ) )
    interMHseq[I] = fa.fetch( c[I] , (e1[I]+1 , s02[I]) )

toc = time.time()
cnts_tsv_df["MH_Sequence"] = MHseq
# cnts_tsv_df["InterMH_Sequence"] = interMHseq     # makes the dataframe huge, and we don't need it
print(f"... done. Extracting sequences took {(toc - tic)/60:0.2f} minutes. Adding GC & modulo 3  ...")
del s1,e1,s02,c,MHseq


cnts_tsv_df["MH_Sequence_GC"] = cnts_tsv_df["MH_Sequence"].apply(GC)
cnts_tsv_df["InterMH_Sequence_GC"] = list(map(GC,interMHseq))
cnts_tsv_df["DuplicatedRegion_GC"] = cnts_tsv_df["MH_Sequence"].add(interMHseq).apply(GC)
cnts_tsv_df["CompleteMHP_GC"] = cnts_tsv_df["MH_Sequence"].add( cnts_tsv_df["MH_Sequence"].add(interMHseq) ).apply(GC)
cnts_tsv_df["Duplication_Creates_modulo3_insertion"] = (list(map(len,interMHseq))+cnts_tsv_df["MH_Sequence"].apply(len)) % 3==0
#cnts_tsv_df["InterMH_Sequence_GC"] = cnts_tsv_df["InterMH_Sequence"].apply(GC)
#cnts_tsv_df["DuplicatedRegion_GC"] = cnts_tsv_df["MH_Sequence"].add(cnts_tsv_df["InterMH_Sequence"]).apply(GC)
#cnts_tsv_df["CompleteMHP_GC"] = cnts_tsv_df["MH_Sequence"].add( cnts_tsv_df["MH_Sequence"].add(cnts_tsv_df["InterMH_Sequence"]) ).apply(GC)
#cnts_tsv_df["Duplication_Creates_modulo3_insertion"] = (cnts_tsv_df["MH_Sequence"].apply(len)+cnts_tsv_df["InterMH_Sequence"].apply(len)) % 3==0

# this  make files smaller 
cnts_tsv_df["MH_Sequence_GC"] = cnts_tsv_df["MH_Sequence_GC"].astype('uint8')
cnts_tsv_df["InterMH_Sequence_GC"] = cnts_tsv_df["InterMH_Sequence_GC"].astype('uint8')
cnts_tsv_df["DuplicatedRegion_GC"] = cnts_tsv_df["DuplicatedRegion_GC"].astype('uint8')
cnts_tsv_df["CompleteMHP_GC"] = cnts_tsv_df["CompleteMHP_GC"].astype('uint8')
cnts_tsv_df["MH_length"] = cnts_tsv_df["MH_length"].astype('uint8')

# other variables needed for predicting MTD frequency
cnts_tsv_df['InterMH_Distance'] = list(map(len,interMHseq))
cnts_tsv_df['InterMH_Distance'] = cnts_tsv_df['InterMH_Distance'].astype('uint16')

# we only need to save the position columns needed to make a .bed file  : s01 & e2
cnts_tsv_df = cnts_tsv_df.drop(columns=['s1', 'e1', 's2', 'e01', 's02', 'e02'])

#give start and end column headers good names
cnts_tsv_df = cnts_tsv_df.rename(columns={"s01": "MHP_start", "e2": "MHP_end"})

# save
tic = time.time()
cnts_tsv_df.to_pickle( args.output_filename )
toc = time.time()
print(f"Saving as pickle took {(toc - tic)/60:0.2f} minutes to file : {args.output_filename :s} ")

tic = time.time()
cnts_tsv_df.to_hdf( "/Users/lcarey/Downloads/df.h5" , key='df', mode='w')
toc = time.time()
print(f"Saving as h5 took {(toc - tic)/60:0.2f} minutes.")

tic = time.time()
cnts_tsv_df.to_feather( "/Users/lcarey/Downloads/df.feather" )
toc = time.time()
print(f"Saving as feather took {(toc - tic)/60:0.2f} minutes.")

