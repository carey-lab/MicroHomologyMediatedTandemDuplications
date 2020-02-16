#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
start_time = Sys.time()

# TRF_2_MHtxt(  genome_file , trf_processed_out_file , df_output_txt_file_name ,  [chr_to_proc] )
# previously run
# trf micro-homology-pair-finding-catching/reference_genome_TRs/data/s_pombe.fa 1 5 1000 80 10 20 100 -f -d -l 2 -ngs -h > s_pombe.trf.txt
# gawk -f micro-homology-pair-finding-catching/reference_genome_TRs/results/extract.awk s_pombe.trf.txt > s_pombe.trf.out


# load data & functions  ---------------------------------------------------------------

source("~/Develop/MicroHomologyMediatedIndels/reference_genome_TRs/code/genome_util_lite.R")
source("~/Develop/MicroHomologyMediatedIndels/reference_genome_TRs/code/tr_util.R")


# get cmd line args
genome_file <- args[1]
trf_processed_out_file <- args[2]
df_output_txt_file_name <- args[3]
if (length(args)==3) {
  chr_to_proc = FALSE
} else{
  chr_to_proc = args[3]
}

#  genome_file = "SS/Port.scaffolds"
# trf_processed_out_file = "SS/Port.trf.out"
# df_output_txt_file_name <-  "SS/Port.txt"

#genome_file <- "micro-homology-pair-finding-catching/reference_genome_TRs/data/s_pombe.fa"
#trf_processed_out_file <- "s_pombe.trf.out"
message( Sys.time() )
message( 'loading genome_file  : ' , genome_file )
gu_sp <- initGenomeFromFasta( genome_file ) 
message(genome_file , ' contains ' , as.character(length(gu_sp$getLength()))  , '  fasta sequence entries')
message( Sys.time() )

# example result:
# GAATTGCAAGCCTCCTCACTGGTAGCCTCACTGGCAGCCTCAACGATG
# consensus = CCTCACTGG[C/T]AG
# mhseq = CCTCA 


# process trf ouytput -----------------------------------------------------
message( 'loading  trf_processed_out_file :   ' , trf_processed_out_file )
pombe_orig <- read.delim( trf_processed_out_file , stringsAsFactors = FALSE )
message( colnames(pombe_orig) )
message(head(pombe_orig))
message( Sys.time() )


# summary <- annotateTR(gu_sp, pombe[i, "chr"], pombe[i, "pos"], pombe[i, "k"])
#  to
# summary <- annotateTR(gu_sp,  as.character(pombe$chr[i]) , pombe[i, "pos"], pombe[i, "k"])
#  Yutze: Use stringsAsFactors = F in the read.delim() may also helps with this

# optionally, restrict to a single chromosome, because the human genome is big
if(chr_to_proc!=FALSE) {
  pombe <- pombe_orig[ pombe_orig$chr == chr_to_proc , ] # for human
} else{
  pombe <- pombe_orig
}

message( 'Finding MH around tandem duplications: df has ' , nrow(pombe) , ' rows' )
print( pombe[1:10,] )
for (i in 1:nrow(pombe)) {
	if( i %% 10 ==0){
  message( as.character(i) , ' '  , pombe[i, "chr"] ,  ' ' ,  as.character(pombe[i, "pos"]) , "\t" , as.character(pombe[i, "k"]) )
	}
  #if( i %% 10 ==0){message(i , "/" , as.character(nrow(pombe)) , "\t" , as.character(round(i/nrow(pombe)*100)) , "%" )}
  #summary <- annotateTR(gu_sp,  as.character(pombe$chr[i]) , pombe[i, "pos"], pombe[i, "k"])
possibleError <- tryCatch({
	summary <- annotateTR(gu_sp, pombe[i, "chr"], pombe[i, "pos"], pombe[i, "k"] )
	pombe$posOriginal[i] <- summary$posOriginal ; pombe$pos[i] <- summary$pos; pombe$k[i] <- summary$k; pombe$rep[i] <- summary$rep;
	pombe$consensus[i] <- summary$cons ; pombe$entropy[i] <- summary$entropy
	pombe$mh[i] <- summary$mh != "" ; pombe$mh_seq[i] <- summary$mh
} , error = function(e) { e 
	print(paste("Oops! --> Error in Loop ",i,sep = ""))
}
)
if( inherits(possibleError , "error") ) next
}

## exclude len2 & len3
pombe$mh_len = nchar(pombe$mh_seq) # MHlen
#pombe$mh[pombe$mh_len==2] <- FALSE
#pombe$mh[pombe$mh_len==3] <- FALSE

pombe <- unique(pombe)
pombe_filter <- data.frame(k = pombe$k < 10, e = pombe$entropy < 1.5, r = pombe$rep < 2)
pombe_use <- pombe[!(pombe_filter$k | pombe_filter$e | pombe_filter$r), ]
off <- 0.025*nrow(pombe_use)

message('finished processing, created the dataframe: ')
message(head(pombe))
message( Sys.time() )

message('saving dataframe to to:   ' , df_output_txt_file_name )

write.table( pombe , df_output_txt_file_name , sep='\t' , row.names=FALSE)
