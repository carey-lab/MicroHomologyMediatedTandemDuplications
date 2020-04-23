setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\new_MH_project\\Processeddata")
load("MHfeatures.SZ=10.Rdata")
####I. DATA PREPARE#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\new_MH_project\\ProcessedData")
#MHfeatures <- read.table("10k.sign.count.tsv.GSM1023946_MN_Mit.features.txt",header = F,as.is = T)
MHfeatures <- read.table("10k.SZ=10.sign.count.tsv.features.txt",header = F,as.is = T)
colnames(MHfeatures) <- c("MHP_start_position","MHlen","interMH","duplication","MHseq","GCcon","InterMHGC")
MHfeatures$dup <- MHfeatures$duplication
MHfeatures$dup[MHfeatures$dup!=0]=1 #binary dupication: turn 1,2,3 to 0   
attach(MHfeatures)
detach(MHfeatures)
save(MHfeatures,file = "MHfeatures.SZ=10.Rdata")
