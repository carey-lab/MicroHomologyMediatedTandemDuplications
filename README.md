## Contains all the scripts for making figures and downstream data-analysis for the paper: 
A rapidly reversible mutation generates subclonal genetic diversity and unstable drug resistance

Lufeng Dan, Yuze Li, Shuhua Chen, Jingbo Liu, Fangting Li, Yu Wang, Xiangwei He*, Lucas B. Carey*

Most genetic changes have negligible reversion rates. As most mutations that confer resistance to an adversary condition (e.g., drug treatment) also confer a growth defect in its absence, it is challenging for cells to genetically adapt to transient environmental changes. Here we identify a set of rapidly reversible drug resistance mutations in S. pombe that are caused by Microhomology mediated Tandem Duplication (MTD), and reversion back to the wild-type sequence. Using 10,000x coverage whole-genome sequencing we identify near 6000 subclonal MTDs in a single clonal population, and determine using machine learning how MTD frequency is encoded in genome. We find that sequences with the highest predicted MTD rates tend to generate insertions that maintain the correct reading frame suggesting that MTD formation has shaped the evolution of coding sequences. Our study reveals a common mechanism of reversible genetic variation that is beneficial for adaptation to environmental fluctuations and facilitates evolutionary divergence.

https://doi.org/10.1101/2020.03.03.972455


## Run on Cluster

Scripts for one-command (actually two-command) running on cluster is presented in the folder `run_on_cluster`.

1. Up load `fastq` data and reference sequences `fasta`
2. Make working directory and creating sub-directories

```
mkdir ref  # Here you put / link your reference sequences
mkdir data # Here you put / link your fasta files
mkdir probes
mkdir mapped
mkdir results
```

3. Do sequence alignment

```
for i in 1 2 3 4 ; do
	alignForMTDCatching "lib"$i "ref/lib"$i".fa" "data/Library_"$i"_clean_R1.fastq.gz" "data/Library_"$i"_clean_R2.fastq.gz"
done
```

Change `1 2 3 4` in above code to the identifier of your own file. Or run the command one-by-one

4. Generate MTD pattern reference and catch MTD reads

```
for i in 1 2 3 4 ; do
    countForSignatures "lib"$i "ref/lib"$i".fa" "mapped/lib"$i".sort.bam"
done
```

## Separated Modules

### 1. Find MH pairs across given sequence(s)

#### Description

Micro-homology pairs (MH pairs) are short identical subsequences flanking a interval of sequnce. For example, in sequence `**ATCG**AGCGCACTTAG**ATCG**`, MH pair with identical sequence `ATCG` are flanking with a interval.

The indel size means the length from the first base of the first MH to the first base of the second MH. This indel size is the size of the indel event when the sequence "duplicated" or "collapsed" induced by the MH.

The first part is to survey a given sequence for all MH pair occurence within limititions. Here we want to find MH pairs with homologous k-mers' size not smaller than 4 bps, interspace greater than 3 pbs, and indel size not over 100 bps.

#### Source

	find_mh.c

#### Usage

	find_mh [output prefix] < [input fasta file]
	e.g:
	find_mh test_data/test < test_data/test.fa

changing the macros for parameter definition will change the limitition settings.

#### Output Format

The program outputs file(s) with name `[output prefix].number.sequence_name.out` contains three fields seperated by tabs. These fields are:

1. starting position of the first MH in a pair
2. starting position of the second MH in a pair
3. the size of the MH (k-mer's k)

### 2. Generate signatures for MTD reads detection

#### Description

The second part is to generate a list of signature sequence features for sequencing reads that carry a MH induced tandem duplication/collapse. Here we used the GNU's `gawk` for scripting.

#### Script

	generate_signatures.awk

#### Usage

Change the `#!` interpreter indicator on the first line in the script to your `gwak` path, and run `chmod u+x` on the script, then you can

	generate_signatures.awk <SZ=[feature size]> <TH=[feature diff threshold]> reference.fa mh.1.out <mh.2.out ...>

Or directly call `gawk` for

	gawk -f generate_signatures.awk <SZ=[feature size]> <TH=[feature diff threshold]> reference.fa mh.1.out <mh.2.out ...>


.....

still editing
