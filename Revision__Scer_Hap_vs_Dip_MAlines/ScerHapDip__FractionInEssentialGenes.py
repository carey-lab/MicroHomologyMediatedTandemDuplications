#%% set constants 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats 
import os

try:
        DATADIR = os.path.dirname(__file__) # folder with this script; also contains the data files
except:
        DATADIR = os.getenv('HOME') + '/Develop/MicroHomologyMediatedTandemDuplications/Revision__Scer_Hap_vs_Dip_MAlines/'
image_file_basename = DATADIR + '/Scer_HapVSDip__MetaData_Ploidy_and_MAlines___'
data_file_name = DATADIR+'/Scer_HapVSDip__dup_sites_found____genes_with_essentiality.txt.xz'
anno_file_name = DATADIR+'/Scer_HapVSDip__MetaData_Ploidy_and_MAlines.csv.xz'

del_fitness_file_name = DATADIR + '/Baryshnikova10.tab'
#%% Load data
# I	2707	7235	intergenic	0	0	I	3725	3743	SRR6997021	0	3732
col_names = ['chr','ORFstart','ORFnd','type','gene','orf','chr2','MHPleftStart','MHPright','SRA','Nreads','MHPleftend']
df = pd.read_csv(data_file_name,sep='\t',names=col_names,keep_default_na=False
        ,dtype={'gene':str,'orf':str})
df['MHlen'] = df.MHPleftend - df.MHPleftStart
df['Essential'] = df.type.str.match('Essential')
df['Non_essen'] = np.logical_or(df.type.str.match('Non-essent'),df.type.str.len()==0)
df['Intergenic'] = df.type.str.match('intergenic')
df['SizeDivByThree'] = (df.MHPright-df.MHPleftStart)%3==0
df['DuplicationLength'] = (df.MHPright-df.MHPleftStart)
df.drop(columns='chr2',inplace=True)
df = df[df.chr!='Mito']
print(df.head())
## ### ### #### ### ### ### ### 
# remove Nreads==0 
#  must check why these exist!
df = df[df.Nreads>0]
df.reset_index(inplace=True,drop=True)
## ### ### #### ### ### ### ### 


A = pd.read_csv(anno_file_name)
A = A[['SRA','Ploidy']]

df = df.merge(A,on='SRA')
fitness_df = pd.read_csv(del_fitness_file_name,sep='\t')[['orf','fitness']]
#df = df.merge(fitness_df,on='orf',how='left')

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
cols_to_use += ['chr','MHPleftStart','MHPright','MHPleftend','MHlen','SizeDivByThree','DuplicationLength'
        ]

hap = df[df.Ploidy == 'hap']
hap = hap.groupby(cols_to_use).agg('count')
hap = hap[hap.Nreads<2]
hap.reset_index(drop=False,inplace=True)

dip = df[df.Ploidy == 'dip']
dip = dip.groupby(cols_to_use).agg('count')
dip = dip[dip.Nreads<2]
dip.reset_index(drop=False,inplace=True)

hap = hap.merge(fitness_df,on='orf',how='left')
dip = dip.merge(fitness_df,on='orf',how='left')
print(len(dip))
print(len(hap))

#%% function to return dict of useful values for bootstrapping
def useful_values(df,id):
        d = {'ID':id,
        'NumEssential':np.sum(df.Essential),
        'PctEssential':np.mean(df.Essential)*100,
        'NumNonEssential':np.sum(df.Non_essen),
        'PctNonEssential':np.mean(df.Non_essen)*100,
        'NumIntergenic':np.sum(df.Intergenic),
        'PctIntergenic':np.mean(df.Intergenic)*100,
        'PctDeleterious':np.sum(df.fitness<0.9) / np.sum(df.fitness>=0)*100,
        'Total':len(df),
        'Div3_pct_Essential':np.mean(df['SizeDivByThree'][df.Essential])*100,
        'Div3_pct_NonEssential':np.mean(df['SizeDivByThree'][df.Non_essen])*100,
        'Div3_pct_Intergenic':np.mean(df['SizeDivByThree'][df.Intergenic])*100,
        'Div3_pct_All':np.mean(df['SizeDivByThree'])*100,
        'Div3_pct_Unfit':np.mean(df['SizeDivByThree'][
                np.logical_and(df.fitness>0,df.fitness<0.9)])*100,
        }
        return d
# %% Bootstrap
BS = pd.DataFrame()
for I in range(0,1000,2):
        BS = BS.append(
                pd.DataFrame(useful_values(dip.sample(frac=1,replace=True),'Diploid'),index=[I])
        )
        BS = BS.append(
                pd.DataFrame(useful_values(hap.sample(frac=1,replace=True),'Haploid'),index=[I+1])
        )

#%% plot based on bootstrap
clrs = ('lightgrey','dimgrey')
g = BS.groupby('ID').agg(['mean','std'])
fig,ax = plt.subplots(2,4)
fig.set_figwidth(12)
ax[0,0].bar([1,2],g.PctEssential['mean'],yerr=g.PctEssential['std'],color=clrs)
ax[0,0].set_ylabel(r'% of MTDs in')
ax[0,0].set_title('Essential genes')
ax[0,2].bar([1,2],g.PctNonEssential['mean'],yerr=g.PctNonEssential['std'],color=clrs)
ax[0,2].set_title('Non-essential genes')
ax[0,3].bar([1,2],g.PctIntergenic['mean'],yerr=g.PctIntergenic['std'],color=clrs)
ax[0,3].set_title('Intergenic regions')
ax[0,1].bar([1,2],g.PctDeleterious['mean'],yerr=g.PctDeleterious['std'],color=clrs)
ax[0,1].set_title('Deleterious mutants')

ax[1,0].bar([1,2],g.Div3_pct_Essential['mean'],yerr=g.Div3_pct_Essential['std'],color=clrs)
ax[1,0].set_ylabel(r'% in-frame')
ax[1,2].bar([1,2],g.Div3_pct_NonEssential['mean'],yerr=g.Div3_pct_NonEssential['std'],color=clrs)
#ax[1,1].set_ylabel(r'% in-frame')
ax[1,3].bar([1,2],g.Div3_pct_Intergenic['mean'],yerr=g.Div3_pct_Intergenic['std'],color=clrs)
ax[1,1].bar([1,2],g.Div3_pct_Unfit['mean'],yerr=g.Div3_pct_Unfit['std'],color=clrs)

#ax[1,2].set_ylabel(r'% in-frame')

for I in range(4):
        ax[1,I].set_yticks([0,33,66,100])
        ax[1,I].axhline(y=33,linestyle=':',color='red')
        ax[1,I].set_ylim(0,100)

for I in range(4):
        for J in range(2):
                try:
                        ax[J,I].set_xticks([1,2])
                        ax[J,I].set_xticklabels(list(g.index.values))
                except:
                        pass

plt.savefig(image_file_basename + 'bar_plots_bootstrap_errorbars.png',dpi=600)
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
labels = ['Essential', 'Non-essen', 'Intergenic', 'Deleterious']
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
