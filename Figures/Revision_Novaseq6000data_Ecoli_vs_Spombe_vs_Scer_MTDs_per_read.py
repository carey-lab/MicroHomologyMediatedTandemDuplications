#%% load the .txt files that summarize data from the cluster
#  reviewers point out that lots of evidence says that MTDs happen in E. coli
#  our data show that MTDs are much less common in E coli
#
#   calculate the number of observed MTDs per million reads for each species
#
# LBC November 2020
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, glob, re

PROJDIR = '/Volumes/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/'
DATADIR = PROJDIR + '/DataFromCluster/'


# %% load metadata (species, lane, library #) by parsing filenames
df = pd.DataFrame()
df['txt_files_full_path'] = glob.glob(DATADIR + '[ES][cp]*/*txt')
df['txt_files_last_part'] = [x.replace(DATADIR,'') for x in df['txt_files_full_path']]
df['species'] = df.txt_files_last_part.str.replace(
    '/.*','',regex=True).astype('category')
df['file_name'] = df.txt_files_last_part.str.replace(
    r'.*/(.*).txt',r'\1',regex=True)
df['SampleN'] = df.file_name.str.replace('_.*','').astype('category')
df['Lane'] = df.file_name.str.replace(r'.*_L(.)\..*',r'\1',regex=True).astype('category')
df['data_type'] = df.file_name.str.replace(r'.*L[12]\.','')

# we no longer need these columns
df.drop(columns='txt_files_last_part',inplace=True)
df.drop(columns='file_name',inplace=True)

# we don't need data from these files
df = df[np.logical_not(df.data_type=='sites_found')].reset_index(drop=True)
df = df[np.logical_not(df.data_type=='good_read_counts')].reset_index(drop=True)

# %% load data -- the first column of first line in each txt file
data = list()
for filename in df.txt_files_full_path:
    fid = open(filename,'r')
    try:
        l = fid.readline().split('\t')
        data.append(int(l[0]))
        print(f'{filename}\t{data[-1]}')
    except:
        print(f'{filename}\tERROR')
        data.append(-99)

df['reads'] = data
df.drop(columns='txt_files_full_path',inplace=True)
print(df)

# %% rename columns & re-index
df.data_type[df.data_type=='q30_f0x2.counts.tsv.good_read_counts'] = 'MTD_counts'
df.data_type[df.data_type=='q30_f0x2.n_mapped_reads'] = 'total_reads'

df.set_index(['species','SampleN','Lane','data_type'],drop=True,inplace=True)

# calculate the ratio (observed MTDs per million mapped reads)
df = df.unstack()
df['MTD_freq'] = df.reads.MTD_counts / df.reads.total_reads * 1e6
print(df)

## %% B11 L2 is missing data. remove (not necessary, as the ratio will be NaN)
#df = df[np.logical_not(df.index.isin([('Spombe','B11','2',slice(None))]))]

#%% boxplots of MTD frequency for each of the three species
fig, ax= plt.subplots(figsize=(4,4))
df.boxplot(by='species',column='MTD_freq',ax=ax,sym='',widths=0.7)
ax.set_ylabel('MTD frequency\n(MTD sites per million reads)')
ax.set_yscale('log')
yt = [0.5,1,5,10,50]
ax.set_yticks(yt)
ax.set_yticklabels(yt)
ax.grid()
plt.savefig('Revision__ComparSpeices__Boxplots__MTD_sites_per_million_reads.png',dpi=600)

print( df.groupby('species')['MTD_freq'].median() )
# %%
