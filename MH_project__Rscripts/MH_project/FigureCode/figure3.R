##dupication frequency as the function of MH length#####
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHfeature <- read.table("10k.sign.count.tsv.features.txt",header = F,as.is = T)
colnames(MHfeature) <- c("duplication","MHlen","interMH","MHseq","entire_nucleSum","entire_nucleMean","entire_nucleMedian","entire_nucleMin","entire_nucleMax") 
attach(MHfeature)

num_Dup <- table(MHlen,dup)
#fre_Dup <- apply(num_Dup[,-1],1,sum)/apply(num_Dup,1,sum)
fre_Dup <- num_Dup[,-1]/apply(num_Dup,1,sum)
x <- 4:(length(fre_Dup)+3)
ba <- plot(x,fre_Dup,type = "l",pch = 20,col = 1, ylab = "Duplication frequency",xlab = "MH length(bp)",main = "Duplication frequency")


