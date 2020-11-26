
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


#%%  Load data & truncate all MHP homology sequence lengths to MAX_MHP_LENGTH
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

filename_MHPs_for_all_genes = 'Spombe__cds_introns_utrs.mhps.xz'
mhps = pd.read_csv(filename_MHPs_for_all_genes,sep='\t',names=['MH1_start','MH2_start','MHPlen','gene'],index_col='gene')

# set max homology sequence length
mhps.MHPlen[mhps.MHPlen>MAX_MHP_LENGTH] = MAX_MHP_LENGTH
print(mhps.head())

#%% Count the number of MHPs, with each MHP homology sequence length, in each gene
g = pd.DataFrame(mhps.groupby(['MHPlen','gene'])['MH1_start'].count())
print(g)

# %% FIG : boxplots of the # of MHPs with each homology sequence length
#  highlighting SSP1 and URA4

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

# %% FIG : Model predictions for all genes, highlighting SSP1 and URA4
fig, ax= plt.subplots(1,1,figsize=figure_size)
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
fasta_file = 'Spombe__cds_introns_utrs.fasta'
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

# %%
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
