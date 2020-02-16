##I.DATA PREPARE#######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
#MHfeatures <- read.table("10k.sign.count.tsv.GSM1023946_MN_Mit.features.txt",header = F,as.is = T)
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features.txt",header = F,as.is = T)
#MHfeatures <- read.table("10k.sign.count.tsv.Nucleosome-density-wtNucWave-reps-median.depth_wl_trimmed_PE2.features.txt",header = F,as.is = T)

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
###/DATA PREPARE

#write.table(MHfeatures,file = "10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features_pro.txt",quote = F,sep = "\t")

##II.dupication frequency as the function of GC content ######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features_pro.txt",header = F,as.is = T)
attach(MHfeatures)
load("MHfeatures.Rdata")

##GC content sort bin
##generate average GC content and frequency data
library(dplyr)
aveGCpre <- function(){
aveGC_pre <- list()
#MH length=4
i=4
GC_Dup <- data.frame(GCcon,dup,MHseq,MHlen)
GC_Dup <- GC_Dup[GC_Dup$MHlen == i,]
GC_Dup <-arrange(GC_Dup,GC_Dup[,"GCcon"])
GC_Dup[,"GCcon"] <- round(GC_Dup[,"GCcon"]*10)/10
cons<- 100000
lab <- c(rep(1:floor(dim(GC_Dup)[1]/cons),each = cons),rep(floor(dim(GC_Dup)[1]/cons)+1,each = dim(GC_Dup)[1] - floor(dim(GC_Dup)[1]/cons)*cons))
GC_Dup_lab <- cbind(GC_Dup,lab)
GC_pre_tab <- table(GC_Dup_lab[,c("GCcon","dup","lab")])

GC_pre <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
aveGC<- c()
for(i in 1:dim(GC_pre_tab)[3]){
  seq <- as.matrix(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHseq"])
  GCnum <- sum(apply(as.matrix(seq),1,function(x){sum(unlist(strsplit(x,''))=="C")+sum(unlist(strsplit(x,''))=="G")}))
  lensum <- sum(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHlen"])
  aveGC <- c(aveGC,GCnum/lensum)
  GC_pre <- rbind(GC_pre,colSums(GC_pre_tab[,,i]))
}

fre <- GC_pre[,2]/(GC_pre[,1]+GC_pre[,2])
GC_pre <- cbind(GC_pre,aveGC,fre)
aveGC_pre$MHlen4 <- GC_pre
##/MHlength = 4

#MH length=5
i=5
GC_Dup <- data.frame(GCcon,dup,MHseq,MHlen)
GC_Dup <- GC_Dup[GC_Dup$MHlen == i,]
GC_Dup <-arrange(GC_Dup,GC_Dup[,"GCcon"])
GC_Dup[,"GCcon"] <- round(GC_Dup[,"GCcon"]*10)/10
cons<- 100000
lab <- c(rep(1:floor(dim(GC_Dup)[1]/cons),each=cons),rep(floor(dim(GC_Dup)[1]/cons)+1,each=(dim(GC_Dup)[1]-floor(dim(GC_Dup)[1]/cons)*cons)))
GC_Dup_lab <- cbind(GC_Dup,lab)
GC_pre_tab <- table(GC_Dup_lab[,c("GCcon","dup","lab")])

GC_pre <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
aveGC<- c()
for(i in 1:dim(GC_pre_tab)[3]){
  seq <- as.matrix(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHseq"])
  GCnum <- sum(apply(as.matrix(seq),1,function(x){sum(unlist(strsplit(x,''))=="C")+sum(unlist(strsplit(x,''))=="G")}))
  lensum <- sum(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHlen"])
  aveGC <- c(aveGC,GCnum/lensum)
  GC_pre <- rbind(GC_pre,colSums(GC_pre_tab[,,i]))
}

fre <- GC_pre[,2]/(GC_pre[,1]+GC_pre[,2])
GC_pre <- cbind(GC_pre,aveGC,fre)
aveGC_pre$MHlen5 <- GC_pre

##MH length = 6
i=6
GC_Dup <- data.frame(GCcon,dup,MHseq,MHlen)
GC_Dup <- GC_Dup[GC_Dup$MHlen == i,]
GC_Dup <-arrange(GC_Dup,GC_Dup[,"GCcon"])
GC_Dup[,"GCcon"] <- round(GC_Dup[,"GCcon"]*10)/10
cons<- 10000
lab <- c(rep(1:floor(dim(GC_Dup)[1]/cons),each=cons),rep(floor(dim(GC_Dup)[1]/cons)+1,each=(dim(GC_Dup)[1]-floor(dim(GC_Dup)[1]/cons)*cons)))
GC_Dup_lab <- cbind(GC_Dup,lab)
GC_pre_tab <- table(GC_Dup_lab[,c("GCcon","dup","lab")])

GC_pre <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
aveGC<- c()
for(i in 1:dim(GC_pre_tab)[3]){
  seq <- as.matrix(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHseq"])
  GCnum <- sum(apply(as.matrix(seq),1,function(x){sum(unlist(strsplit(x,''))=="C")+sum(unlist(strsplit(x,''))=="G")}))
  lensum <- sum(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHlen"])
  aveGC <- c(aveGC,GCnum/lensum)
  GC_pre <- rbind(GC_pre,colSums(GC_pre_tab[,,i]))
}

fre <- GC_pre[,2]/(GC_pre[,1]+GC_pre[,2])
GC_pre <- cbind(GC_pre,aveGC,fre)
aveGC_pre$MHlen6 <- GC_pre
##/MH length =6

##MH length >= 7
i=7
GC_Dup <- data.frame(GCcon,dup,MHseq,MHlen)
GC_Dup <- GC_Dup[GC_Dup$MHlen >= i,]
GC_Dup <-arrange(GC_Dup,GC_Dup[,"GCcon"])
GC_Dup[,"GCcon"] <- round(GC_Dup[,"GCcon"]*10)/10
cons<- 10000
lab <- c(rep(1:floor(dim(GC_Dup)[1]/cons),each=cons),rep(floor(dim(GC_Dup)[1]/cons)+1,each=(dim(GC_Dup)[1]-floor(dim(GC_Dup)[1]/cons)*cons)))
GC_Dup_lab <- cbind(GC_Dup,lab)
GC_pre_tab <- table(GC_Dup_lab[,c("GCcon","dup","lab")])

GC_pre <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
aveGC<- c()
for(i in 1:dim(GC_pre_tab)[3]){
  seq <- as.matrix(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHseq"])
  GCnum <- sum(apply(as.matrix(seq),1,function(x){sum(unlist(strsplit(x,''))=="C")+sum(unlist(strsplit(x,''))=="G")}))
  lensum <- sum(GC_Dup_lab[GC_Dup_lab[,"lab"]==i,"MHlen"])
  aveGC <- c(aveGC,GCnum/lensum)
  GC_pre <- rbind(GC_pre,colSums(GC_pre_tab[,,i]))
}

fre <- GC_pre[,2]/(GC_pre[,1]+GC_pre[,2])
GC_pre <- cbind(GC_pre,aveGC,fre)
aveGC_pre$MHlen7 <- GC_pre
##/MH length >=7
aveGC_pre
}
aveGC_pre <- aveGCpre()
##/generate average GC content and frequency data

##generate aveGC content and frequency data frame
aveGC_frame <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for(i in 1:4){
  label<- rep(i, each=dim(aveGC_pre[[i]])[1])
  #lab_ave <- cbind(label,aveGC_pre[[i]][,3:4])
  aveGC_frame <- rbind(aveGC_frame,cbind(label,aveGC_pre[[i]][,3:4]))
}
##/generate aveGC content and frequency data frame

##plot
plotGC <- function(){
lm.reg<- list()
x <- as.numeric(aveGC_pre[[1]][,3]*100)
y <- log10(aveGC_pre[[1]][,4]*100)
ba <- plot(x,y,pch = 20,col = 1,ylim = c(),ylab = "Duplication %",xlab = "average GC content (%)",main = "GC content-log10(Duplication %)")
#lines(x,y,type = "l",pch = 20,lwd=1.5,col = 1, ylab = "Duplication frequency",xlab = "GC content",main = "Duplication frequency")

for(i in 2:length(aveGC_pre))
{
  x <- as.numeric(aveGC_pre[[i]][,3]*100)
  y <- 10^(log10(aveGC_pre[[i]][,4]*100))
  points(x,y,pch = 20,col = i, lwd=1,ylab = "Duplication frequency",xlab = "GC content",main = "Duplication frequency")
  lm.reg[i]<-lm(y~x)#linear regression to calculate slope
  abline(lm(y~x),col = i)
  }
legend("bottomright", legend=c("MHlen=4","MHlen=5","MHlen=6","MHlen>=7"), col=c(1:length(aveGC_pre)), lty=1)
lm.reg
}
plotGC()
##/plot

##plot
library(ggplot2)
p1 <- function(){
  ggplot(data = aveGC_frame, aes(x = aveGC*100, y = fre*100))+ggtitle("GC content-Duplication percent")+
  geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
  geom_smooth(data=aveGC_frame,aes( colour = paste0("MHlen=",label+3)),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid") +
  labs(x="average GC content (%)",y="Duplication percent (%)",fill = "MH length (bp)")+
  #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
   scale_y_log10(breaks=c(.0001,.001,.01,.1,1),labels=c(.0001,.001,.01,.1,1)
                 ,expand = c(0.2,0.00001)
     #breaks = trans_breaks("log10", function(x) 10^x),
     #labels = trans_format("log10", math_format(10^.x))
                ) +
  theme_bw() +
  
  theme(
    plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
    axis.title.x = element_text(color="black", size=18),
    axis.title.y = element_text(color="black", size=18),
    #delete background
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
    # panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    #???????????????
    axis.line = element_line(colour = "black",size=0.5),
    #?????????
    axis.ticks = element_line(size=0.5),
    axis.ticks.length=unit(-0.2,"cm"),
    
    #x?????????
    axis.text.x = element_text(size=18,color='black',margin = margin(t=0.3, unit = "cm")),
    #y?????????
    axis.text.y = element_text(size=18,color='black',margin = margin(r=0.3, unit = "cm")),
    #??????
    #legend.title = element_text(colour="black", size=14),
    legend.title = element_blank(),
    #legend.text = element_blank()
    legend.text = element_text(colour=c("black"), size=16)
    # remove legend
    ,legend.justification=c(1,0), legend.position=c(0.95,0.05),
    legend.background = element_rect(fill="white",
                                     size=0.5, linetype="solid", 
                                     colour ="black")
    #,legend.position="none"
  ) 
}
p1()
##/plot

##doubles duplication probility
dupfold <- list()
for(i in 1:4){
  a <- aveGC_frame[aveGC_frame$label==i,]
 dupfold_ <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for (j in 1:dim(a)[1]) {
    for (k in 1:dim(a)[1]) {
      if(a[k,3]/a[j,3]>=2){
      dupfold_ <- rbind(dupfold_,c(a[k,2],a[j,2],a[k,3],a[j,3],a[k,3]/a[j,3]))
      colnames(dupfold_) <- c("GCcon1","GCcon2","dupfre1","dupfre2","fold") 
      }
      k=k+1
    }
    j=j+1
  }
  
  dupfold[[i]]<- dupfold_  
}
dupfold
##/doubles duplication probility

##/GC content bin

##/dupication frequency as the function of GC content



##IV.dupication frequency as the function of nicleosome occupancy######
nucle_Dup <- data.frame(entire_nucleMedian,dup)
nucle_Dup[,1] <- round(nucle_Dup[,1]*2)/2
nucle_pre <- table(nucle_Dup)
nucle_pre[,2]/nucle_pre[,1]
x <- as.numeric(row.names(nucle_pre))
y <- nucle_pre[,2]/(nucle_pre[,1]+nucle_pre[,2])
ba <- plot(x[1:30],y[1:30],type = "l",pch = 20,col = 1, ylab = "Duplication frequency",xlab = "nucleosome occupancy",main = "Duplication frequency")
##/dupication frequency as the function of gene expression


