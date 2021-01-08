#!/usr/bin/env python3
#  Reviewer writes:   I'd like to get some sense of how many reads are supporting each of the MTDs. 
#%% import
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

FIGDIR = os.getenv('HOME') + '/Downloads/'
FIGNAME = FIGDIR + 'MTD_revisions__how_many_reads_support_each_MTD_in_the_MA_lines.png'

DATADIR = '/Volumes/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/'
FN = 'all_MTDs_found_in_MA_lines.txt'
df = pd.read_csv(DATADIR + FN, sep='\t',names=['chr','s','e','type','NReads','SRA'])

df = df[df.type.str.match('DUP_')]
df.sort_values(['chr','s','e'],inplace=True)
dfI = df.set_index(['chr','s','e','type'])

g = dfI.groupby(['chr','s','e','type'])['SRA'].count().rename('N_SRAs')

dfI = dfI.merge(pd.DataFrame(g), left_index=True,right_index=True,how='outer')
dfI.sort_values('NReads',ascending=False,inplace=True)

# haploid vs diploid
DATADIR = os.getenv('HOME') + '/Develop/MicroHomologyMediatedTandemDuplications/Revision__Scer_Hap_vs_Dip_MAlines/'
anno_file_name = DATADIR+'/Scer_HapVSDip__MetaData_Ploidy_and_MAlines.csv.xz'
A = pd.read_csv(anno_file_name, usecols=['SRA','Ploidy'])
dfI = dfI.merge(A,how='left',on='SRA')

# for each MTD, count the number of SRAs where it's found
q = dfI.groupby(['NReads','N_SRAs']).count().rename(columns={'SRA':'N_MTDs'}).reset_index()
q.sort_values('N_MTDs',inplace=True)
q.drop(columns='Ploidy',inplace=True)
q['N_MTDs_log10']  = np.log10(q.N_MTDs)
#%% plot only for haploid, or only for diploid
# for each MTD, count the number of SRAs where it's found
# q = dfI[dfI.Ploidy=='hap'].groupby(['NReads','N_SRAs']).count().rename(columns={'SRA':'N_MTDs'}).reset_index()
# q.sort_values('N_MTDs',inplace=True)
# q.drop(columns='Ploidy',inplace=True)
# q['N_MTDs_log10']   = np.log10(q.N_MTDs)
# %%
plt.tight_layout()
sns.set_style("white")
ax = sns.scatterplot(x="N_SRAs", y="NReads", hue="N_MTDs_log10", palette='viridis', data=q)
norm = plt.Normalize(q['N_MTDs_log10'].min(), q['N_MTDs_log10'].max())
sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
sm.set_array([])
ax.get_legend().remove()
cbh = ax.figure.colorbar(sm,ticks=[0,1,2])
cbh.ax.set_yticklabels(['1','10','100'])
cbh.ax.set_title('# of MTDs')
ax.set_ylabel("# of reads supporting the MTD")
ax.set_xlabel("# of M.A. lines with the MTD")
ax.set_xticks([1,43,175,q.N_SRAs.max()])
ax.figure.set_size_inches((4,4))
ax.figure.savefig(FIGNAME,dpi=600)
plt.show()

# %%
#%% 
fig,ax = plt.subplots(1,1,figsize=(4,4))
ax.hist(dfI.NReads,bins=np.linspace(0.5,17.5,18), label='all MTDs')
ax.hist( dfI.loc[g.index[g>10]]['NReads'] , bins=np.linspace(0.5,17.5,18),histtype='step',linewidth=2, label='MTDs found in >1 M.A. line')

ax.set_xticks([1,5,10,15])
ax.set_xlabel('# of reads supporting the MTD')
ax.set_ylabel('# of MTDs')
ax.legend(loc='nort east')
#ax.set_yscale('log')