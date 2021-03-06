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
parser.add_argument('--output_basename', help='will save results into two files: <filename>.feather & <filename>.bed' , required=True)
parser.add_argument('--head', help='only process the first N MHPs. Default: 0, process all.' , default=None , type=int)
parser.add_argument('--calc-inter_MHP_GC', help='calculate GC content for inter-MHP sequences' , default=False  , action='store_true')
parser.add_argument('--testing', help='set paramters for quick debugging' , default=False  , action='store_true')


args = parser.parse_args()

## for development
#os.chdir('/Users/lcarey/Downloads')
if (args.testing==True):
    args.fasta = '/Users/lcarey/CareyLab/ExternalData/PomBase/pombe_n.fasta'
    args.countstsv = '/Users/lcarey/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/SRR7817502_n.sign.count_500bp.tsv'
    args.head = 1000 


print(args.fasta, ' is ', os.path.getsize(args.fasta), ' bytes')
fa = pyfastx.Fasta(args.fasta)
for s in fa:
    print( "sequences in the FASTA file : ", s.id , " , name=" , s.name  , " desc=" , s.description , 'l=' , len(s.seq) )

print(args.countstsv, ' is ', os.path.getsize(args.countstsv), ' bytes')
tic = time.time()
cnts_tsv_df = pandas.read_csv(args.countstsv, sep='\t', names=["chr", "s1", "e1", "s2", "e2", "NDups", "NCol"] , nrows=args.head)
toc = time.time()
print(f"... loaded. Adding Sequence data took {toc - tic:0.2f} seconds")

## smaller df for debugging 
#cnts_tsv_df = cnts_tsv_df.loc[ (cnts_tsv_df["chr"]=="I") ]
#cnts_tsv_df = cnts_tsv_df.loc[ (cnts_tsv_df["chr"]=="I") | (cnts_tsv_df["chr"]=="II")  | (cnts_tsv_df["chr"]=="III")  ] 
#cnts_tsv_df = cnts_tsv_df.loc[ (cnts_tsv_df["s1"] < 1000000) | (cnts_tsv_df["NDups"] > 0) ] 
#cnts_tsv_df = cnts_tsv_df.reset_index(drop=True) # rebuild row numbers in the shrunken dataframe

cnts_tsv_df["s01"] = cnts_tsv_df.s1 - 1
cnts_tsv_df["e01"] = cnts_tsv_df.e1 - 1
cnts_tsv_df["s02"] = cnts_tsv_df.s2 - 1
cnts_tsv_df["e02"] = cnts_tsv_df.e2 - 1

# put the position columns into the order needed to make a .bed file  : s01 & e2
cnts_tsv_df["MH_length"] = cnts_tsv_df.e1 - cnts_tsv_df.s1 + 1
cnts_tsv_df = cnts_tsv_df[['chr','s01','e2','MH_length','NDups','s1','s2','e1','e01','s02','e02']]
print( cnts_tsv_df.head() )




# Assign MHseq, MHlen, & GC
# way way way way faster to not work with dataframes
#   extract all info into lists, save into a list, and put the result into the df

tic = time.time()
MHseq = ["" for x in range(len(cnts_tsv_df))]
interMHseq = ["" for x in range(len(cnts_tsv_df))]

c=cnts_tsv_df['chr'].to_list()
s1=cnts_tsv_df['s1'].to_list()
e1=cnts_tsv_df['e1'].to_list()
s02=cnts_tsv_df['s02'].to_list()
print('Extracting sequences...')
for I in range(len(cnts_tsv_df)):
    MHseq[I] = fa.fetch( c[I] , (s1[I] , e1[I] ) )
    if (args.calc_inter_MHP_GC):
        interMHseq[I] = fa.fetch( c[I] , (e1[I]+1 , s02[I]) )

toc = time.time()
cnts_tsv_df["MH_Sequence"] = MHseq
# cnts_tsv_df["InterMH_Sequence"] = interMHseq     # makes the dataframe huge, and we don't need it
print(f"... done. Extracting sequences took {(toc - tic)/60:0.2f} minutes. Adding GC & modulo 3  ...")
del s1,e1,s02,c,MHseq


# addd GC content to df
cnts_tsv_df["MH_Sequence_GC"] = cnts_tsv_df["MH_Sequence"].apply(GC)
cnts_tsv_df["MH_Sequence_GC"] = cnts_tsv_df["MH_Sequence_GC"].astype('uint8')

cnts_tsv_df["MH_length"] = cnts_tsv_df["MH_length"].astype('uint8')


if (args.calc_inter_MHP_GC):
    cnts_tsv_df["InterMH_Sequence_GC"] = list(map(GC,interMHseq))
    cnts_tsv_df["DuplicatedRegion_GC"] = cnts_tsv_df["MH_Sequence"].add(interMHseq).apply(GC)
    cnts_tsv_df["CompleteMHP_GC"] = cnts_tsv_df["MH_Sequence"].add( cnts_tsv_df["MH_Sequence"].add(interMHseq) ).apply(GC)
    cnts_tsv_df["InterMH_Sequence_GC"] = cnts_tsv_df["InterMH_Sequence_GC"].astype('uint8')
    cnts_tsv_df["DuplicatedRegion_GC"] = cnts_tsv_df["DuplicatedRegion_GC"].astype('uint8')
    cnts_tsv_df["CompleteMHP_GC"] = cnts_tsv_df["CompleteMHP_GC"].astype('uint8')


# we only need to save the position columns needed to make a .bed file  : s01 & e2
cnts_tsv_df = cnts_tsv_df.drop(columns=['s1', 'e1', 's2', 'e01', 's02', 'e02'])

#give start and end column headers good names
cnts_tsv_df = cnts_tsv_df.rename(columns={"s01": "MHP_start", "e2": "MHP_end"})


# other variables needed for predicting MTD frequency
cnts_tsv_df['InterMH_Distance'] = (cnts_tsv_df['MHP_end'] - cnts_tsv_df['MHP_start']  - 2 * cnts_tsv_df['MH_length']).astype('uint16')

cnts_tsv_df["Duplication_Creates_modulo3_insertion"] = (cnts_tsv_df['InterMH_Distance'] + cnts_tsv_df["MH_length"]) % 3==0


#sort 
cnts_tsv_df = cnts_tsv_df.sort_values( by=['chr','MHP_start','MHP_end'] , ascending=True )
cnts_tsv_df = cnts_tsv_df.reset_index(drop=True)

# save
tic = time.time()
cnts_tsv_df.to_feather( args.output_basename + '.feather' )
toc = time.time()
print(f"Saving as feather took {(toc - tic)/60:0.2f} minutes to file : {args.output_basename :s} ")

# output a .bed file with column header
cnts_tsv_df = cnts_tsv_df.rename(columns={'chr' : '#chr'})
cnts_tsv_df.to_csv( args.output_basename + '.bed' , sep='\t' , header=True , index=False) 
