#!/usr/bin/env python3

#%% import & set paths
import pysam
import numpy as np 
import tempfile, argparse, os, datetime
from pprint import pprint
from pybedtools import BedTool

parser = argparse.ArgumentParser()
parser.add_argument('--input_bam',type=str,required=True)
parser.add_argument('--output_bam',type=str,required=True)
parser.add_argument('--ORIs_anno_file',type=str,required=True)

parser.add_argument('--watson',action='store_true',default=False)
parser.add_argument('--crick',action='store_true',default=False)
parser.add_argument('--leading',action='store_true',default=False)
parser.add_argument('--lagging',action='store_true',default=False)

parser.add_argument('--kb_to_keep',type=int,default=15)
parser.add_argument('--min_mapping_qual',type=int,default=20)
parser.add_argument('--min_avg_base_qual_score',type=int,default=20)
opts = parser.parse_args()

nt_to_take = opts.kb_to_keep*1000

def decide_if_to_keep_read(this_read):
    if ( not this_read.cigarstring=='150M' and 
        this_read.is_proper_pair and
        this_read.mapping_quality >= opts.min_mapping_qual and 
        not this_read.is_qcfail and 
        np.mean(this_read.query_qualities) >= opts.min_avg_base_qual_score 
    ) : return True

if os.path.getsize(opts.input_bam)<10000:
    raise Exception(f"{opts.input_bam} does not exist or is too small.")
ORIs = BedTool(opts.ORIs_anno_file)


bam_fh = pysam.AlignmentFile(opts.input_bam, "rb") # read in from this file
tmpfn = tempfile.NamedTemporaryFile(suffix='.bam')
tmp_fh = pysam.AlignmentFile(tmpfn.name,'wb',header=bam_fh.header) # write out to this file

# extract leading:  read 1 maps to the :  
# the  Crick  strand, 15kb 5’ of the ORI
# the Watson strand, 15kb 3’ of the ORI

# get leading
for interval in ORIs:
    read_counter = 0
    if (opts.watson and opts.leading) or (opts.crick and opts.lagging) :
        iter = bam_fh.fetch(interval.chrom,max([interval.start-nt_to_take,0]),interval.start)
        for x in iter:
            if decide_if_to_keep_read(x): 
                tmp_fh.write(x)
                read_counter += 1
    elif (opts.crick and opts.leading) or (opts.watson and opts.lagging) :
        iter = bam_fh.fetch(interval.chrom,interval.stop,interval.stop+nt_to_take)
        for x in iter:
            if decide_if_to_keep_read(x): 
                tmp_fh.write(x)
                read_counter += 1
    else:
        pprint(opts)
        raise Exception(f'invalid combination of W/C and Lead/Lag')
    print(f"{read_counter}\t{datetime.datetime.now()}\t{interval}",end='')

tmp_fh.close()
pysam.sort("-m", "8G","-o", opts.output_bam, tmpfn.name)
pysam.index(opts.output_bam)

# get lagging
#iter = C.fetch('chrXIV',ars_loc-nt_to_take,ars_loc)
#for x in iter:
#    if decide_if_to_keep_read(x): lagging.write(x)
#iter = W.fetch('chrXIV',ars_loc,ars_loc+nt_to_take)
#for x in iter:
#    if decide_if_to_keep_read(x): lagging.write(x)
#pysam.sort("-o", lagging_reads_bam_fn, tmpfn2.name)


# %%
