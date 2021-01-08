
#%% import & set constants
# Reviewers wonderd why MTDs hadn't been picked up in 5-FOA screens using URA4
#   Likely they have been, but have not been noticed
#  one reviewer proposed that URA4 is depleted for microhomology pairs (MHPs)
# lets compare the total number of MHPs in URA4 vs SSP1
#
# LBC:  I do NOT think we should normalize by gene length. The probability of a mutation occuring  in a gene
#   is higher for longer genes. 
#
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

ura3 = 'YEL021W'
ura4 = 'SPCC330.05c'
ssp1 = 'SPCC297.03'
MAX_MHP_LENGTH = 8

figure_size  = (7, 2)

figure_output_path = '/Users/lcarey/Downloads/'
boxplot_figure_name  = figure_output_path + 'Revision__URA4_depleted_MHPS_boxplot.png'

#  load model predictions for all genes
PROJDIR = '/Volumes/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/'
DATADIR = PROJDIR + '/Sarah/MH_project/Manuscript-todo/processeddata/'
PREDFILE = DATADIR + 'MHRSumPreinGene.txt'
model_prediction_figure_name = figure_output_path + 'Revision_SSP1_vs_URA4_ModelPredictionScores_CDFs.png'
df = pd.read_csv(PREDFILE,sep='\t')
df.set_index('systematic_name',drop=True,inplace=True)


#  Load data & truncate all MHP homology sequence lengths to MAX_MHP_LENGTH
#
#  ids_to_genes = pd.read_csv('ids_2_genes.tab',sep='\t') # not needed unless we want to map systematic ORF names to common gene names
#  Spombe__cds_introns_utrs.mhps is from running find_mh on the fasta file containing UTRs, introns, and exons, from PomBase
#  'ftp://ftp.pombase.org/pombe/genome_sequence_and_features/feature_sequences/cds%2Bintrons%2Butrs.fa.gz'
#  
#   MHPS_FILE = Spombe__cds_introns_utrs.mhps
#   count_mhps_per_gene:
#     rm -f $(MHPS_FILE)
#     find_mh $(MHPS_FILE) < Spombe__cds_introns_utrs.fasta
#     find ./ -type f -name \*out -exec cat {} >> $(MHPS_FILE) \;
#     find ./ -type f -name \*out -exec rm -f {} \;
#     xz -e -T 10 Spombe__cds_introns_utrs.mhps
DDIR2 = '/Volumes/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Revision_Ura4_MHP_Frequency/'
filename_MHPs_for_all_genes = DDIR2 + 'Spombe__cds_introns_utrs.mhps.xz'
mhps = pd.read_csv(filename_MHPs_for_all_genes,sep='\t',names=['MH1_start','MH2_start','MHPlen','gene'],index_col='gene')

mhps_scer = pd.read_csv( DDIR2 + 'Scer_genomic.mhps.xz',sep='\t',names=['MH1_start','MH2_start','MHPlen','gene'],index_col='gene')

# set max homology sequence length
mhps_scer.MHPlen[mhps_scer.MHPlen>MAX_MHP_LENGTH] = MAX_MHP_LENGTH
mhps.MHPlen[mhps.MHPlen>MAX_MHP_LENGTH] = MAX_MHP_LENGTH
print(mhps.head())

# Count the number of MHPs, with each MHP homology sequence length, in each gene
g = pd.DataFrame(mhps.groupby(['MHPlen','gene'])['MH1_start'].count()).reset_index()
g.MH1_start = g.MH1_start.astype(float)
print(g)

gScer = pd.DataFrame(mhps_scer.groupby(['MHPlen','gene'])['MH1_start'].count()).reset_index()
gScer.MH1_start = gScer.MH1_start.astype(float)
#%% MPHs per KB to normalize by gene length

gene_lengths = mhps.reset_index().groupby('gene').max()
gene_lengths.MH2_start += gene_lengths.MHPlen
gene_lengths.drop(columns=['MH1_start','MHPlen'],inplace=True)
gene_lengths.rename(columns={'MH2_start':'GeneLength'},inplace=True)
gene_lengths.reset_index(inplace=True)
gene_lengths = q.merge(gene_lengths,on='gene')
gene_lengths['MHPs_per_kb'] = gene_lengths.MH1_start / gene_lengths.GeneLength / 1000
gene_lengths.set_index(['MHPlen','gene'],inplace=True)
gene_lengths.drop(columns=['MH1_start','GeneLength'],inplace=True)
# %% FIG : boxplots of the # of MHPs with each homology sequence length
#  highlighting SSP1 and URA4

fig, ax= plt.subplots(figsize=figure_size)
g.boxplot(by='MH1_start',ax=ax,sym='',whis=(10,90))
ax.set_xlabel('MHP length (nt)') 
ax.set_ylabel('# of MHPs in gene')
ax.set_yscale('log')
print(type(ax))


ssp1_data = g.loc[(slice(None),ssp1),:]
x = np.array([x for x,n in ssp1_data.index])
y = ssp1_data.MH1_start.to_list()
ax.plot(x-3,y,'o',label='SSP1',color='purple')

ura4_data = g.loc[(slice(None),ura4),:]
x = np.array([x for x,n in ura4_data.index])
y = ura4_data.MH1_start.to_list()
ax.plot(x-3,y,'o',label='URA4',color='red')

ax.set_title('# of MHPs in S. pombe genes')
ax.set_xticklabels([4,5,6,7,'>=8'])
ax.legend(edgecolor='black')
ax.grid()
plt.savefig(boxplot_figure_name,dpi=600)

#%% S. cerevisiae boxplot

fig, ax= plt.subplots(figsize=figure_size)
_ = gScer[['MHPlen','MH1_start']].boxplot(by='MHPlen',ax=ax,sym='',whis=(10,90))
ax.set_xlabel('MHP length (nt)') 
ax.set_ylabel('# of MHPs in gene')
ax.set_yscale('log')

ura3_data = gScer[gScer.gene==ura3]
x = ura3_data.MHPlen.values
y = ura3_data.MH1_start.values
ax.plot(x-3,y,'o',label='URA3',color='orange')
ax.legend(edgecolor='black')
ax.grid()
ax.set_xticklabels([4,5,6,7,'>=8'])

plt.savefig(boxplot_figure_name+'URA3.png',dpi=600)
#%%  boxplots of total # of MHPs in each gene
fig, ax= plt.subplots(figsize=figure_size)
g.boxplot(by='MHPlen',ax=ax,sym='',whis=(10,90))
ax.set_xlabel('MHP length (nt)') 
ax.set_ylabel('# of MHPs in gene')
ax.set_yscale('log')
print(type(ax))


ssp1_data = g.loc[(slice(None),ssp1),:]
x = np.array([x for x,n in ssp1_data.index])
y = ssp1_data.MH1_start.to_list()
ax.plot(x-3,y,'o',label='SSP1',color='purple')

ura4_data = g.loc[(slice(None),ura4),:]
x = np.array([x for x,n in ura4_data.index])
y = ura4_data.MH1_start.to_list()
ax.plot(x-3,y,'o',label='URA4',color='red')

ax.set_title('# of MHPs in S. pombe genes')
ax.set_xticklabels([4,5,6,7,'>=8'])
ax.legend(edgecolor='black')
ax.grid()
plt.savefig(boxplot_figure_name,dpi=600)


#%%  boxplots of MHPs in each gene normalized by gene length
fig, ax= plt.subplots(figsize=(7,6))
gene_lengths.boxplot(by='MHPlen',ax=ax,sym='',whis=(10,90))
ax.set_xlabel('MHP length (nt)') 
ax.set_ylabel('# of MHPs per kb')
ax.set_yscale('log')

ssp1_data = gene_lengths.loc[(slice(None),ssp1),:]
x = np.array([x for x,n in ssp1_data.index])
y = ssp1_data.MHPs_per_kb.to_list()
ax.plot(x-3,y,'o',label='SSP1',color='purple')

ura4_data = gene_lengths.loc[(slice(None),ura4),:]
x = np.array([x for x,n in ura4_data.index])
y = ura4_data.MHPs_per_kb.to_list()
ax.plot(x-3,y,'o',label='URA4',color='red')

ax.set_title('# of MHPs per kb')
ax.set_xticklabels([4,5,6,7,'>=8'])
ax.legend(edgecolor='black')
ax.grid()
#plt.savefig(boxplot_figure_name,dpi=600)




# %% FIG : Model predictions for all genes, highlighting SSP1 and URA4
fig, ax= plt.subplots(1,1,figsize=(4,4))
plt.sca(ax)
plt.hist(df.sumFloat,bins=1000,label='all genes',color='black',cumulative=True,histtype='step',density=True)
plt.vlines(df.sumFloat.loc[ura4],0,1,label='URA4',color='red',linewidth=3)
plt.vlines(df.sumFloat.loc[ssp1],0,1,label='SSP1',color='purple',linewidth=3)
ax.set_xlabel('Prediction Score')
ax.set_ylabel('Fraction of genes')
ax.set_yticks(np.linspace(0,1,11))
ax.grid()
ax.set_xlim(0,220)
ax.set_ylim(0,1)
ax.legend(loc='lower right')

plt.savefig(model_prediction_figure_name,dpi=600)



# %%  also look at other paramters, eg: gene length, GC content
fasta_file = DDIR2 + 'Spombe__cds_introns_utrs.fasta'
name = list()
seq = list()
GC = list()
gene_len = list()

for cur_gene in SeqIO.parse(fasta_file,'fasta'):
    name.append(cur_gene.name)
    seq.append(cur_gene.seq)
    gc = float((cur_gene.seq.count('C') + cur_gene.seq.count('G')) / len(cur_gene.seq))
    gene_len.append(len(cur_gene.seq))
    GC.append(gc)

df = pd.DataFrame(zip(name,GC,gene_len),columns=['name','GC','gene_len'])
df.set_index('name',drop=True,inplace=True)
print(df)

fig, ax= plt.subplots(2,1,figsize=(10, 4))
plt.sca(ax[0])
plt.hist(df.GC,bins=100,label='all genes',color='grey')
plt.vlines(df.GC.loc[ura4],0,500,label='URA4',color='red',linewidth=3)
plt.vlines(df.GC.loc[ssp1],0,500,label='SSP1',color='purple',linewidth=3)
ax[0].set_xlabel('% GC')
ax[0].set_ylabel('# of genes')

plt.sca(ax[1])
plt.hist(df.gene_len,bins=100,label='all genes',color='grey')
plt.vlines(df.gene_len.loc[ura4],0,500,label='URA4',color='red',linewidth=3)
plt.vlines(df.gene_len.loc[ssp1],0,500,label='SSP1',color='purple',linewidth=3)
ax[1].set_xlabel('gene length (nt)')
ax[1].set_ylabel('# of genes')
# %%
fig, ax= plt.subplots(2,1,figsize=(10, 6))
plt.sca(ax[0])
plt.hist(df.GC,bins=100,label='all genes',color='black',cumulative=True,histtype='step',density=True)
plt.vlines(df.GC.loc[ura4],0,1,label='URA4',color='red',linewidth=3)
plt.vlines(df.GC.loc[ssp1],0,1,label='SSP1',color='purple',linewidth=3)
ax[0].set_xlabel('% GC')
ax[0].set_title('% GC')
ax[0].set_ylabel('# of genes')
ax[0].set_yticks(np.linspace(0,1,11))
ax[0].grid()
ax[0].set_xlim(0.25,0.6)
ax[0].legend()

plt.sca(ax[1])
plt.hist(df.gene_len,bins=100,label='all genes',color='black',cumulative=True,histtype='step',density=True)
plt.vlines(df.gene_len.loc[ura4],0,1,label='URA4',color='red',linewidth=3)
plt.vlines(df.gene_len.loc[ssp1],0,1,label='SSP1',color='purple',linewidth=3)
ax[1].set_xlabel('gene length (nt)')
ax[1].set_ylabel('# of genes')
ax[1].set_yticks(np.linspace(0,1,11))
ax[1].grid()
#ax[1].set_xscale('log')
ax[1].set_xlim(250,5000)
ax[1].legend()

plt.savefig('Revision_SSP1_vs_URA4_GC_GeneLength__CDFs.png',dpi=600)


# %%
