#%% set constants 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats 
import os


DATADIR = os.path.dirname(__file__) # folder with this script; also contains the data files
image_file_basename = DATADIR + '/Scer_HapVSDip__MetaData_Ploidy_and_MAlines___'
data_file_name = DATADIR+'/Scer_HapVSDip__dup_sites_found____genes_with_essentiality.txt.xz'
anno_file_name = DATADIR+'/Scer_HapVSDip__MetaData_Ploidy_and_MAlines.csv.xz'

#%% Load data
# I	2707	7235	intergenic	0	0	I	3725	3743	SRR6997021	0	3732
col_names = ['chr','ORFstart','ORFnd','type','gene','orf','chr2','MHPleftStart','MHPright','SRA','Nreads','MHPleftend']
df = pd.read_csv(data_file_name,sep='\t',names=col_names)
df['MHlen'] = df.MHPleftend - df.MHPleftStart
df['Essential'] = df.type.str.match('Essential')
df['Non_essen'] = df.type.str.match('Non-essent')
df['Intergenic'] = df.type.str.match('intergenic')
df['SizeDivByThree'] = (df.MHPright-df.MHPleftStart)%3==0 ; 
df.drop(columns='chr2',inplace=True)
df = df[df.chr!='Mito']

## ### ### #### ### ### ### ### 
# remove Nreads==0 
#  must check why these exist!
df = df[df.Nreads>0]
df.reset_index(inplace=True,drop=True)
## ### ### #### ### ### ### ### 


A = pd.read_csv(anno_file_name)
A = A[['SRA','Ploidy']]

df = df.merge(A,on='SRA')

df.head()
# %% count in how many SRAs we observe each MTD
#  and only keep MTDs observed in only a single SRA
#
# command line equivalents: 
#  total number:
#      xz -dkc Scer_HapVSDip__dup_sites_found____genes_with_essentiality.txt.xz | wc -l 
#         4853
#
#  number of unique MTDs:
#      xz -dkc Scer_HapVSDip__dup_sites_found____genes_with_essentiality.txt.xz | cut -f 1-9 |sort|uniq -c|sort -nk 1|wc -l
#           993
#
cols_to_use = ['ORFstart','ORFnd','type','orf','gene','Essential','Non_essen','Intergenic']
cols_to_use += ['chr','MHPleftStart','MHPright','MHPleftend','MHlen','SizeDivByThree']

hap = df[df.Ploidy == 'hap']
hap = hap.groupby(cols_to_use).agg('count')
hap = hap[hap.Nreads<2]
hap.reset_index(drop=False,inplace=True)

dip = df[df.Ploidy == 'dip']
dip = dip.groupby(cols_to_use).agg('count')
dip = dip[dip.Nreads<2]
dip.reset_index(drop=False,inplace=True)

print(len(dip))
print(len(hap))

# %%
print('')
print(f"# Essential & diploid = {np.sum(dip.Essential)}")
print(f"# Essential & haploid = {np.sum(hap.Essential)}")
print('')
print(f"# Non-Essential & diploid = {np.sum(dip.Non_essen)}")
print(f"# Non-Essential  & haploid = {np.sum(hap.Non_essen)}")
print('')
print(f"# intergenic & diploid = {np.sum(dip.Intergenic)}")
print(f"# intergenic & haploid = {np.sum(hap.Intergenic)}")
print('')
print(f"total # dip = {len(dip)}")
print(f"total # hap = {len(hap)}")
print('')


print(f"Dip mean Div3 Essential = {np.mean(dip['SizeDivByThree'][dip.Essential])}")
print(f"Hap mean Div3 Essential = {np.mean(hap['SizeDivByThree'][hap.Essential])}")
print(f"Dip mean Div3 Non-Essen = {np.mean(hap['SizeDivByThree'][dip.Non_essen])}")
print(f"Hap mean Div3 Non-Essen = {np.mean(hap['SizeDivByThree'][hap.Non_essen])}")
print(f"Dip mean Div3 Intergenic = {np.mean(hap['SizeDivByThree'][dip.Intergenic])}")
print(f"Hap mean Div3 Intergenic = {np.mean(hap['SizeDivByThree'][hap.Intergenic])}")


def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:d}".format(absolute)
    #return "{:.1f}%\n({:d})".format(pct, absolute)

fig, ax = plt.subplots(1,2)
labels = ['Essential', 'Non-essen', 'Intergenic']
sizesDip = [14,32,65]
sizesHap = [1,10,25]
explode = (0.1,0,0)  # only "explode" the 2nd slice (i.e. 'Hogs')
# autopct='%1.1f%%'
ax[0].pie(sizesDip, explode=explode, labels=labels, autopct=lambda pct: func(pct, sizesDip),
        shadow=True, startangle=90)
ax[0].set_title('Diploid')

ax[1].pie(sizesHap, explode=explode, labels=labels, autopct=lambda pct: func(pct, sizesHap),
        shadow=True, startangle=90)
ax[1].set_title('Haploid')
                
# %%
x = [[sizesDip[0], sizesHap[0]],[sum(sizesDip[-2:]),sum(sizesHap[-2:])] ]
oddsratio, pvalue = stats.fisher_exact(x)
print(f"FE test: OR={oddsratio}\tp={pvalue:0.05f}")
# %%
