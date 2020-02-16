#Make a single script that loads all the data, saves a single data frame with all the data needed for prediction into one .Rdata file
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\Manuscript-todo\\processeddata")
load("MHfeatures.Rdata")
####I. DATA PREPARE#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
#MHfeatures <- read.table("10k.sign.count.tsv.GSM1023946_MN_Mit.features.txt",header = F,as.is = T)
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features.txt",header = F,as.is = T)
#MHfeatures <- read.table("10k.sign.count.tsv.Nucleosome-density-wtNucWave-reps-median.depth_wl_trimmed_PE2.features.txt",header = F,as.is = T)

withingene <- read.table("10k.sign.count.tsv.sorted_full_MHR__interect100pct__intersectANY.txt",header = F,as.is = T,fill = T)
MHfeatures$withingenes <- c(0)
MHfeatures$withingenes[(withingene[,6]==0)&(withingene[,7]!=0)]=1
MHfeatures$withingenes[(withingene[,6]!=0)&(withingene[,7]!=0)]=2

closest <- read.table("distance_to_closest_MHR_with_dup.bed",header = F,as.is = T,fill = T)
MHfeatures$ntclosestMHR <- closest[,5]
MHfeatures$logntclosestMHR <- log10(closest[,5])
MHfeatures$bintclosestMHR <- 0
MHfeatures$bintclosestMHR[MHfeatures$ntclosestMHR<=100] <- 1

colnames(MHfeatures) <- c("duplication","MHlen","interMH","MHseq","entire_nucleSum","entire_nucleMean","entire_nucleMedian","entire_nucleMin","entire_nucleMax","entire_geneSum","entire_geneMean","entire_geneMedian","entire_geneMin","entire_geneMax") 
MHfeatures$dup <- MHfeatures$duplication
MHfeatures$dup[MHfeatures$dup!=0]=1 #binary dupication: turn 1,2,3 to 0   
attach(MHfeatures)

##f.interMH:transformed intreMH
num_inter_Dup <- table(data.frame(interMH,dup))
x.name <- as.numeric(row.names(num_inter_Dup))
y <- num_inter_Dup[,2]/num_inter_Dup[,1]
which(y==max(y))
constant <- x.name[which(y==max(y))]
MHfeatures$f.interMH <- abs(interMH-constant)#transform interMH feature:monotanic model:f(inter_MHlen) = (inter_MHlen) - 170(which interMHlen has most duplications frequency)
##/f.interMH

##GC content
GC <- apply(as.matrix(MHseq),1,function(x){sum(unlist(strsplit(x,''))=="C")+sum(unlist(strsplit(x,''))=="G")})
MHfeatures$GCcon <- as.numeric(unlist(GC/MHlen))# caculate GC content
#/GC content

##log2(gene expression)
MHfeatures$logentire_geneMedian <- log2(MHfeatures$entire_geneMedian)
MHfeatures$logentire_geneMedian[is.infinite(MHfeatures$logentire_geneMedian)]=NA
MHfeatures$logentire_geneSum <- log2(MHfeatures$entire_geneSum)
MHfeatures$logentire_geneSum[is.infinite(MHfeatures$logentire_geneSum)]=NA
MHfeatures$logentire_geneMean <- log2(MHfeatures$entire_geneMean)
MHfeatures$logentire_geneMean[is.infinite(MHfeatures$logentire_geneMean)]=NA
MHfeatures$logentire_geneMax <- log2(MHfeatures$entire_geneMax)
MHfeatures$logentire_geneMax[is.infinite(MHfeatures$logentire_geneMax)]=NA
MHfeatures$logentire_geneMin <- log2(MHfeatures$entire_geneMin)
MHfeatures$logentire_geneMin[is.infinite(MHfeatures$logentire_geneMin)]=NA
#/log2(gene expression)

detach(MHfeatures)
save(MHfeatures,file = "MHfeatures.Rdata")
###/DATA PREPARE
