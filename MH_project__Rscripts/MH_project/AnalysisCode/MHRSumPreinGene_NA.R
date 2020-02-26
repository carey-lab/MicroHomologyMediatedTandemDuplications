setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHRSumPreinGene_NA <- read.table("MHRSumPreinGene_NA.txt",header = F,fill = T)
colnames(MHRSumPreinGene_NA) <- c("chr","start1","end2","gene_length","gene","systematic_name","sumObs","sumPre","sumFloat")
sumFloat<-c()
sumNA <- c()
NAprecent <- c()
for(i in 1:dim(MHRSumPreinGene_NA)[1]){
  d <- as.numeric(unlist(strsplit(as.character(MHRSumPreinGene_NA$sumFloat[i]),",",fixed = TRUE)))
  sumFloat[i] <- round(sum(d,na.rm = T),4)
  sumNA[i] <- sum(is.na(d))
  NAprecent[i]<- round(sum(is.na(d))/length(d),4)
}
MHRSumPreinGene_NA$sumFloat <- sumFloat
MHRSumPreinGene_NA$NAprecent <- NAprecent
MHRSumPreinGene_NA$sumFloat_div_geneLen <- round(MHRSumPreinGene_NA$sumFloat/MHRSumPreinGene_NA$gene_length,4)
library(dplyr)
MHRSumPreinGene <- arrange(MHRSumPreinGene_NA,desc(sumFloat_div_geneLen))
write.table(MHRSumPreinGene,file = "MHRSumPreinGene.txt",col.names = T,row.names = F,quote = F,sep = "\t")
