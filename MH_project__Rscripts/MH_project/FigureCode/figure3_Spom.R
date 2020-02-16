setwd("E:\\work\\SynologyDrive\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHlen_Dup <- read.table("10k.sign.count.tsv.MHlen_Dup.txt",header = F, comment.char = "#",as.is = T)
num_Dup <- table(MHlen_Dup)
fre_Dup <- apply(num_Dup[,-1],1,sum)/apply(num_Dup,1,sum)
ba <- plot(4:(length(fre_Dup)+3),fre_Dup,type = "l",pch = 20,col = 1, ylab = "Duplication frequency",xlab = "MH length(bp)",main = "Duplication frequency")


pombegenome <- read.table("pombegenome.fasta",header = F, comment.char = "#",as.is = T)





seq <- c("GCAG","AGAA","TTGAA","GAAAAG","AAAAG","AGCTAG","GATA","AGTA","AGTA","TAGAT","AATT","GCTA","GCTA","TATG","TGTTAA","TAAA","AATT","AATT","TACA","ACAT","ATAAA","TAAAA","AAATT","AATT","AATTC","TTCA","GTAT","ATAA","ATAA","TAATA","ATATT","ATTA","TTAA","TTAAT")
seq <- as.matrix(rep(seq,1000),ncol=1)
MHlen_Dup[1:34,1]
unlist(strsplit(seq[1],''))

MH_D <- MHlen_Dup[1:34000,]
MH_D[MH_D[,2]!=0,2]=1

GC <- apply(seq,1,function(x){sum(unlist(strsplit(x,''))=="C")+sum(unlist(strsplit(x,''))=="G")})
GC_Dup <- data.frame(GC/MH_D[,1],MH_D[,2])
GC_pre <- table(GC_Dup)
GC_pre[,2]/GC_pre[,1]

ba <- plot(as.numeric(row.names(GC_pre)),GC_pre[,2]/GC_pre[,1],type = "l",pch = 20,col = 1, ylab = "Duplication frequency",xlab = "MH length(bp)",main = "Duplication frequency")


