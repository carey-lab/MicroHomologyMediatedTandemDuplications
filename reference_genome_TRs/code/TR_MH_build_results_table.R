#!/usr/bin/env Rscript
library('optparse')
#source("~/Develop/MicroHomologyMediatedTandemDuplications/reference_genome_TRs/code/genome_util_lite.R")
#source("~/Develop/MicroHomologyMediatedTandemDuplications/reference_genome_TRs/code/tr_util.R")
source("~/Develop/MicroHomologyMediatedIndels/reference_genome_TRs/code/genome_util_lite.R")
source("~/Develop/MicroHomologyMediatedIndels/reference_genome_TRs/code/tr_util.R")

option_list = list( 
    make_option( c( "-f" , "--fasta_file") , type="character" , default=NULL , help='FASTA file name') , 
    make_option( c( "-t" , "--trf_output_file") , type="character" , default=NULL , help='TRF output txt table file') , 
    make_option( c( "-o" , "--output_file") , type="character" , default=NULL , help='save output table as this file (tab delim)') ,
    make_option( c( "-i" , "--included_fasta_seq_in_output_table") , type="logical" , default=TRUE , help='optional : include the entire sequence in the .tab output file') 
    
    )
opt = parse_args( OptionParser( option_list=option_list ) )

#opt$fasta_file = '/Users/lcarey/Downloads/test/micro-homology-pair-finding-catching/reference_genome_TRs/data/s_pombe.fa'
#opt$trf_output_file = '/Users/lcarey/Downloads/test/micro-homology-pair-finding-catching/reference_genome_TRs/results/s_pombe.out'
#opt$output_file = 'test.txt'

## prepare genome utils
fa <- initGenomeFromFasta( opt$fasta_file ) # 
chr_lengths = fa$getLength()
message(opt$fasta_file  , ' contains ' , as.character(length(fa$getLength()))  , '  fasta sequence entries')

## read TR utils
# read TRF output summarized table produced by extract.awk
df <- read.delim( opt$trf_output_file , stringsAsFactors = FALSE ) # "results/s_pombe.out"
message( colnames(df) )
message(head(df))

# revise the result using homemade TR util
for (i in 1:nrow(df)) {
    if( i %% 10 ==0){
        message( as.character(i) , ' '  , df[i, "chr"] ,  ' ' ,  as.character(df[i, "pos"]) , "\t" , as.character(df[i, "k"]) )
    }
    possibleError <- tryCatch({
    summary <- annotateTR(fa, df[i, "chr"], df[i, "pos"], df[i, "k"])
    df$posOriginal[i] <- summary$posOriginal
    df$pos[i] <- summary$pos
    df$k[i] <- summary$k
    df$rep[i] <- summary$rep;
    df$consensus[i] <- summary$cons
    df$entropy[i] <- summary$entropy
    df$mh[i] <- summary$mh != "" 
    df$mh_seq[i] <- summary$mh
    if(opt$included_fasta_seq_in_output_table){ df$sequence[i] <- fa$extract(df[i, "chr"],1,as.numeric(chr_lengths[df[i, "chr"]])) }
    } , error = function(e) { e
        print(paste("Oops! --> Error in Loop ",i,sep = ""))
    }
    )
    if( inherits(possibleError , "error") ) next

}
df <- unique(df)

write.table( df , file=opt$output_file , sep="\t" , row.names=FALSE , col.names=TRUE , quote = FALSE)