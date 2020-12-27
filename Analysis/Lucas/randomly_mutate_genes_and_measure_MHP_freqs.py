#%% import
import gzip, argparse, sys, random, os
from Bio import SeqIO, Seq
from Levenshtein import hamming
from Bio.SeqUtils import GC, CodonUsage, CodonUsageIndices
from subprocess import run, PIPE
from os import path
from io import StringIO
import pandas as pd 
import numpy as np
from scipy.stats import gaussian_kde
from pprint import pprint

sys.argv = [''] #if running in interactive window
parser = argparse.ArgumentParser()
parser.add_argument('--number_of_shuffled_versions',default=20,type=int)
parser.add_argument('--cds_fasta_filename',default='cds.fasta.gz',type=str)
parser.add_argument('--min_gene_length',default=100,type=int)
parser.add_argument('--MAX_MHPLEN',default=7,type=int)

opts = parser.parse_args()

#cai = CodonUsage.CodonAdaptationIndex()
#cai.set_cai_index(CodonUsageIndices.SharpEcoliIndex)

#%% CodonShuffler (taken from Wilke's codon_tools)
class CodonShuffler:
    """
    Class to shuffle codons within an ORF, without changing the overall
    codon frequency.
    """
    def __init__(self):
        pass

    def make_lookup_table(self, seq):
        """
        Make a reverse codon lookup table that lists codons used in seq for
        each amino acid. The same codon may be listed more than once for each
        amino acid.
        """
        codons = {}

        for i in range(int(len(seq)/3)):
            codon = seq[3*i:3*i+3]
            amino_acid = str(Seq.translate(codon))
            if amino_acid in codons:
                codons[amino_acid].append(str(codon))
            else:
                codons[amino_acid] = [str(codon)]

        return codons

    def shuffle_codons(self, seq, start_window = 0, end_window = 0):
        """
        Shuffle codons in a sequence without changing codon frequencies.
        """
        # Trim sequence based on start and end windows
        start = start_window*3
        end = len(seq) - end_window*3
        seq_trim = seq[start:end]

        # Make sure sequence doesn't have any partial codons
        assert len(seq_trim) % 3 == 0

        # Construct dictionary of all codons present in the sequence
        codons = self.make_lookup_table(seq_trim)
        # Shuffle codons
        for amino_acid in codons:
            random.shuffle(codons[amino_acid])
        # Translate original sequence
        seq_aa = seq_trim.translate()
        shuffled_seq = seq[0:start]
        # Reconstruct sequence with codons shuffled
        for amino_acid in str(seq_aa):
            shuffled_seq += codons[amino_acid].pop()
        shuffled_seq += seq[end:]
        return shuffled_seq

# %% shuffle all genes
shuffler = CodonShuffler()
c=0
l = pd.DataFrame()
results_df = pd.DataFrame()
with gzip.open(opts.cds_fasta_filename,'rt') as fasta_fh:
    for gene in SeqIO.parse(fasta_fh,'fasta'):
        c+=1
        if gene.seq.startswith('ATG') and gene.seq.translate().endswith('*') \
            and (gene.seq.translate().count('*')==1) and (len(gene.seq)>opts.min_gene_length):
            d = dict()
            d['GC'] = GC(gene.seq)
            d['ID'] = gene.id
            d['len'] = len(gene.seq)
            #and (not gene.id in EssentialGenes):
            #print(f"HD={hamming(str(shuffled_seq),str(gene.seq))}")
            #print(f"CAI orig={cai.cai_for_gene(str(gene.seq)):0.3f}\tnew={cai.cai_for_gene(str(shuffled_seq)):0.3f}")
            input_seq_orig = f">{gene.id}\n{gene.seq}"
            p = run(['find_mh','STDOUT'], stdout=PIPE,input=input_seq_orig, encoding='ascii')
            df = pd.read_csv(StringIO(p.stdout),sep='\t',names=['s','e','MHPlen','gene'])
            df.MHPlen[df.MHPlen>opts.MAX_MHPLEN] = opts.MAX_MHPLEN
            d['N_MHPs_Orig'] = df.groupby('MHPlen')['gene'].count().rename('# of MHPs')

            t = pd.DataFrame()
            for i in range(opts.number_of_shuffled_versions):
                shuffled_seq = shuffler.shuffle_codons(gene.seq)
                input_seq_recoded = f">{gene.id + '_shuffled'}\n{shuffled_seq}"
                p = run(['find_mh','STDOUT'], stdout=PIPE,input=input_seq_recoded, encoding='ascii')
                df = pd.read_csv(StringIO(p.stdout),sep='\t',names=['s','e','MHPlen','gene'])
                df.MHPlen[df.MHPlen>opts.MAX_MHPLEN] = opts.MAX_MHPLEN
                n_mhps = df.groupby('MHPlen')['gene'].count().rename('# of MHPs') 
                t = t.append(pd.DataFrame(n_mhps))
            t = t.reset_index().groupby('MHPlen').mean()
            d['N_MHPs_Shuff'] = t.squeeze()
            results_df = results_df.append(d,ignore_index=True)
            if c&100==0: print('.',end='')
            if c&1000==0: print(' ')

            # create ratio series with genes as index; can be unstacked into df
            r = (d['N_MHPs_Orig'] / d['N_MHPs_Shuff'])
            r.index = (pd.MultiIndex.from_tuples([(gene.id,a) for a in r.index]))
            l = l.append(pd.DataFrame(r))

results_df = l.unstack()
results_df['median'] = results_df.median(axis=1)
results_df['mean'] = results_df.mean(axis=1)
results_df.sort_values('median',inplace=True,ascending=False)


#%% join with metadata
meta = pd.read_csv('Spombe_gene_names_map.tsv.xz',sep='\t',index_col='Systematic ID')
q = results_df.merge(meta,how='left',left_index=True,right_index=True)
q.to_csv('MHP_enrichment.tsv',sep='\t')
# %% histograms for paper
xl = np.linspace(-0.5,0.5,1000)
bw = 0.1
lw=3
plt.figure(figsize=(4,4))
x = [x.loc[4] for i,x in enumerate(l) if i%2==0]
y = [x.loc[4] for i,x in enumerate(l) if i%2==1]
density = gaussian_kde(np.log2(np.array(x)/np.array(y)))
density.covariance_factor = lambda : bw
density._compute_covariance()
plt.plot(xl,density(xl),label='4',linewidth=lw)

plt.vlines(0,0,density(xl).max(),linestyles='dashed',color='grey')

x = [x.loc[5] for i,x in enumerate(l) if i%2==0]
y = [x.loc[5] for i,x in enumerate(l) if i%2==1]
density = gaussian_kde(np.log2(np.array(x)/np.array(y)))
density.covariance_factor = lambda : bw
density._compute_covariance()
plt.plot(xl,density(xl),label='5',linewidth=lw)

x = [x.loc[6] for i,x in enumerate(l) if i%2==0]
y = [x.loc[6] for i,x in enumerate(l) if i%2==1]
density = gaussian_kde(np.log2(np.array(x)/np.array(y)))
density.covariance_factor = lambda : bw
density._compute_covariance()
plt.plot(xl,density(xl),label='6',linewidth=lw)

x = [x.loc[7] for i,x in enumerate(l) if i%2==0]
y = [x.loc[7] for i,x in enumerate(l) if i%2==1]
density = gaussian_kde(np.log2(np.array(x)/np.array(y)))
density.covariance_factor = lambda : bw
density._compute_covariance()
plt.plot(xl,density(xl),label='>=7',linewidth=lw)

plt.xlim([-0.5,0.5])
plt.xlabel('log2( native / shuffled )')
plt.ylabel('Density of genes')
plt.legend(loc='upper left')
plt.savefig(FIGPATH + 'MHP enrichment in native coding sequences - all_genes.png', dpi=600)
plt.show()

# %%
# %% ECDFs
xl = np.linspace(-0.5,0.5,1000)
x = [x.loc[4] for i,x in enumerate(l) if i%2==0]
y = [x.loc[4] for i,x in enumerate(l) if i%2==1]
plt.hist(np.log2(np.array(x)/np.array(y)),density=True,bins=xl,cumulative=True,histtype='step',label='4')

x = [x.loc[5] for i,x in enumerate(l) if i%2==0]
y = [x.loc[5] for i,x in enumerate(l) if i%2==1]
plt.hist(np.log2(np.array(x)/np.array(y)),density=True,bins=xl,cumulative=True,histtype='step',label='5')

x = [x.loc[6] for i,x in enumerate(l) if i%2==0]
y = [x.loc[6] for i,x in enumerate(l) if i%2==1]
plt.hist(np.log2(np.array(x)/np.array(y)),density=True,bins=xl,cumulative=True,histtype='step',label='6')

x = [x.loc[7] for i,x in enumerate(l) if i%2==0]
y = [x.loc[7] for i,x in enumerate(l) if i%2==1]
plt.hist(np.log2(np.array(x)/np.array(y)),density=True,bins=xl,cumulative=True,histtype='step',label='>=7')


plt.grid()
plt.legend(loc='upper left')
plt.xlim([-0.5,0.45])
plt.yticks(np.linspace(0,1,11))
plt.xlabel('log2( native / shuffled )')
plt.ylabel('Fraction of genes')


# %%