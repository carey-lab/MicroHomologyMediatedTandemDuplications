setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData\\SRR")
load("MN_Spombe_972h-_Rep1.MHfeatures.Rdata") 
#load("GSM1023946_MN_Mit.MHfeatures.Rdata")
####I. DATA PREPARE#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData\\SRR")
#MHfeatures <- read.table("SRR7817502_rm.sign.count.tsv.GSM1023946_MN_Mit.features.txt",header = F,as.is = T,fill=T)
MHfeatures <- read.table("SRR7817502_rm.sign.count.tsv.MN_Spombe_972h-_Rep1.features.txt",header = F,as.is = T,fill = T)
#MHfeatures <- read.table("SRR7817502_rm.sign.count.tsv.Nucleosome-density-wtNucWave-reps-median.depth_wl_trimmed_PE2.features.txt",header = F,as.is = T,fill=T)
colnames(MHfeatures) <- c("duplication","MHlen","interMH","MHseq","entire_nucleSum","entire_nucleMean","entire_nucleMedian","entire_nucleMin","entire_nucleMax","entire_geneSum","entire_geneMean","entire_geneMedian","entire_geneMin","entire_geneMax") 

closest <- read.table("distance_to_closest_MHR_with_dup_SRR.bed",header = F,as.is = T,fill = T)
MHfeatures$ntclosestMHR <- closest[,5]
MHfeatures$logntclosestMHR <- log10(closest[,5])
MHfeatures$bintclosestMHR <- 0
MHfeatures$bintclosestMHR[MHfeatures$ntclosestMHR<=100] <- 1

#colnames(MHfeatures) <- c("duplication","MHlen","interMH","MHseq","entire_nucleSum","entire_nucleMean","entire_nucleMedian","entire_nucleMin","entire_nucleMax","entire_geneSum","entire_geneMean","entire_geneMedian","entire_geneMin","entire_geneMax") 
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
save(MHfeatures,file = "MN_Spombe_972h-_Rep1.MHfeatures.Rdata")
###/DATA PREPARE



####II.LOGISTICS REFRESSION#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\Manuscript-todo\\processeddata")
#load("MHfeatures.Rdata")

####(1)balanced data
library(ROCR)
library(caret)
library(pROC)
library(plotROC)
predata <- function(x){
  #prepare data
  event <- MHfeatures[MHfeatures$dup==1,]
  set.seed(x)
  a <- table(MHfeatures$dup)
  num <- as.numeric(a[names(a)==1])
  nonevent <- MHfeatures[sample(which(MHfeatures$dup==0),num),]
  MHfeatures_te <- data.frame(rbind(event,nonevent))
  attach(MHfeatures_te)
  dup<- as.factor(MHfeatures_te$dup)
  MHfeatures_te <- data.frame(scale(subset(MHfeatures_te,select=-c(dup,MHseq))))#MHseq:MHR sequence scale:normalise
  MHfeatures_te <- data.frame(dup,MHfeatures_te)
  colnames(MHfeatures_te) 
  detach(MHfeatures_te)
  MHfeatures_te
  #/prepare data
}
MHfeatures_te <- predata(1)#randomely choose some non-events,seed =1
##/balanced dat

##(2)train the model and draw ROC
fold_pre1 <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+logentire_geneMedian+logntclosestMHR
                 +interGCcon,
                 data = MHfeatures_te,family = "binomial")
fold_predict1 <- predict(fold_pre1,type='response',newdata=MHfeatures_te)
fold_pre2<- glm(formula = dup ~ MHlen + f.interMH + logentire_geneMedian,
                data = MHfeatures_te,family = "binomial")
fold_predict2 <- predict(fold_pre2,type='response',newdata=MHfeatures_te)
auc1 <- auc(as.numeric(MHfeatures_te$dup),fold_predict1)
auc2 <- auc(as.numeric(MHfeatures_te$dup),fold_predict2)
name <- rep(c("full","top3"),each=length(MHfeatures_te$dup))
test <- data.frame(rbind(cbind(as.numeric(MHfeatures_te$dup),fold_predict1),cbind(as.numeric(MHfeatures_te$dup),fold_predict2)),name)
colnames(test)<- c("dup","pre","name")                   

library(ggplot2)
p8 <- function(){
  ggplot(test, aes(d = dup, m = pre, color = name)) + ggtitle("haploid data")+
    #geom_text(x= 0.75, y =0.5, label=paste0("AUC of full model is ",round(auc1,4),"\nAUC of top3 model is ",round(auc2,4)))+
    annotate("text", x = 0.65, y=0.5, size=6,label = paste0("AUC of full model is ",round(auc1,4),"\nAUC of top3 model is ",round(auc2,4)))+
    geom_roc(n.cuts = 0,size=1.3)+
    labs(y="Sensitivity", x="1-Specificity")+
    theme_bw() +
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=22),
      axis.title.y = element_text(color="black", size=22),
      #delete background
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      
      #???????????????
      axis.line = element_line(colour = "black",size=0.5),
      #?????????
      axis.ticks = element_line(size=0.5),
      axis.ticks.length=unit(-0.2,"cm"),
      
      #x?????????
      axis.text.x = element_text(size=22,color='black',margin = margin(t=0.3, unit = "cm")),
      #y?????????
      axis.text.y = element_text(size=22,color='black',margin = margin(r=0.3, unit = "cm")),
      #??????
      #legend.title = element_text(colour="black", size=14),
      legend.title = element_blank(),
      #legend.text = element_blank()
      legend.text = element_text(colour=c("black"), size=16)
      # remove legend
      ,legend.justification=c(1,0), legend.position=c(0.95,0.75),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p8()

##/train the model and draw ROC
###copmare2models####
d_p <- as.numeric(coef(summary(diploid))[,4])
h_p <- as.numeric(coef(summary(haploid))[,4])
Padj=c(d_p,h_p)
label=rep("NS",14)
label[Padj<0.05]<-"*"
label[Padj<0.01]<-"**"
label[Padj<0.001]<-"***"


df<- data.frame(supp=rep(c("Dip", "Hap"), each=7),
                features=rep(c("intercept","MHlen","GCcon","interMH","nucleOcc","geneExp","closestMHR"),time=2),
                coef=c(as.numeric(diploid$coefficients),as.numeric(haploid$coefficients)),
                Padj,label)
df

library(ggplot2)
p8 <- function(){
  ggplot(data = df, aes(x = features, y = coef, fill=supp,label=label)) + 
    #geom_text(x= 0.75, y =0.5, label=paste0("AUC of full model is ",round(auc1,4),"\nAUC of top3 model is ",round(auc2,4)))+
    #annotate("text", x = 0.65, y=0.5, size=6,label = paste0("AUC of full model is ",round(auc1,4),"\nAUC of top3 model is ",round(auc2,4)))+
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(data=df,label=label,size = 5)+
    scale_fill_brewer(palette="Paired")+
    theme_minimal()+
    theme_bw() +
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=22),
      #axis.text.x=element_text(angle =-90, vjust = 0.5),
      axis.title.y = element_text(color="black", size=22),
      #delete background
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      
      #???????????????
      axis.line = element_line(colour = "black",size=0.5),
      #?????????
      axis.ticks = element_line(size=0.5),
      axis.ticks.length=unit(-0.2,"cm"),
      
      #x?????????
      axis.text.x = element_text(size=22,color='black',margin = margin(t=0.3, unit = "cm")),
      #y?????????
      axis.text.y = element_text(size=22,color='black',margin = margin(r=0.3, unit = "cm")),
      #??????
      #legend.title = element_text(colour="black", size=14),
      legend.title = element_blank(),
      #legend.text = element_blank()
      legend.text = element_text(colour=c("black"), size=16)
      # remove legend
      ,legend.justification=c(1,0), legend.position=c(0.95,0.75),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p8()

##II.dupication frequency as the function of GC content ######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData\\SRR")
load("MN_Spombe_972h-_Rep1.MHfeatures.Rdata")
attach(MHfeatures)

##GC content sort bin
##(1)generate average GC content and frequency data
library(dplyr)
aveGCpre <- function(){
  aveGC_pre <- list()
  #MH length=4
  i=4
  GC_Dup <- data.frame(GCcon,dup,MHseq,MHlen)
  GC_Dup <- GC_Dup[GC_Dup$MHlen == i,]
  GC_Dup <-arrange(GC_Dup,GC_Dup[,"GCcon"])
  GC_Dup[,"GCcon"] <- round(GC_Dup[,"GCcon"]*10)/10
  cons<- 1000000
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
  cons<- 1000000
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
  aveGC_pre$MHlen6 <- GC_pre
  ##/MH length =6
  
  ##MH length >= 7
  i=7
  GC_Dup <- data.frame(GCcon,dup,MHseq,MHlen)
  GC_Dup <- GC_Dup[GC_Dup$MHlen >= i,]
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
  aveGC_pre$MHlen7 <- GC_pre
  ##/MH length >=7
  aveGC_pre
}
aveGC_pre <- aveGCpre()
##/generate average GC content and frequency data

##(2)generate aveGC content and frequency data frame
aveGC_frame <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for(i in 1:4){
  label<- rep(i, each=dim(aveGC_pre[[i]])[1])
  #lab_ave <- cbind(label,aveGC_pre[[i]][,3:4])
  aveGC_frame <- rbind(aveGC_frame,cbind(label,aveGC_pre[[i]][,3:4]))
}
##/generate aveGC content and frequency data frame

##(3)plot
library(ggplot2)
p1 <- function(){
  ggplot(data = aveGC_frame, aes(x = aveGC*100, y = fre*100))+ggtitle("GC content-Dup(haploid)")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=aveGC_frame,aes( colour = paste0("MHlen=",label+3)),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="average GC content (%)",y="% of MHPs with an MTD",fill = "MH length (bp)")+
    guides(color = guide_legend(reverse=T))+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(.001,.01,.1,1)
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
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
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
      ,legend.justification=c(1,0), legend.position=c(0.95,0.7),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p1()
##/plot



##/GC content bin
##/dupication frequency as the function of GC content



###II.dupication frequency as the function of interMH ######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData\\SRR")
load("MN_Spombe_972h-_Rep1.MHfeatures.Rdata")
attach(MHfeatures)
len_inter_Dup <- data.frame(MHlen,interMH,dup)
len_inter_Dup$group <- cut(len_inter_Dup[,"interMH"], breaks = seq(from=0,to=500,by=5))
len_tab <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for (i in 1:3) {
  if(i+3!=6){
    ind <- MHlen==i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_inter_Dup[ind, ]))
  }
  if(i+3==6){
    ind <- MHlen>=i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_inter_Dup[ind, ]))
  }
  
}
num_inter_Dup <- table(len_tab[,c("label","group","dup")])

##generate interMH and frequency data frame
inter_frame <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for(i in 1:3){
  label<- rep(i, each=dim(num_inter_Dup[i,,])[1])
  inter <- seq(from=5,to=500,by=5)
  fre <- num_inter_Dup[i,,2]/(num_inter_Dup[i,,1]+num_inter_Dup[i,,2])*100
  inter_frame <- rbind(inter_frame,cbind(label,inter,fre,num_inter_Dup[i,,2],num_inter_Dup[i,,1]))
}
##/generate interMH and frequency data frame

##plot
library(ggplot2)
p2 <- function(){
  ggplot(data = inter_frame, aes(x = inter, y = fre))+ggtitle("interMH-Dup(haploid)")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=inter_frame,aes( colour = paste0("MHlen=",label+3)),method='loess',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="nucleotides between MHR(bp)",y="% of MHPs with an MTD",fill = "MH length (bp)")+
    guides(color = guide_legend(reverse=T))+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(.0001,.001,.01,.1,1),labels=c(.0001,.001,.01,.1,1)
                  ,expand = c(0.2,0.0001)
                  #breaks = trans_breaks("log10", function(x) 10^x),
                  #labels = trans_format("log10", math_format(10^.x))
    ) +
    xlim(0, 300)+
    theme_bw() +
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=18),
      axis.title.y = element_text(color="black", size=18),
      #delete background
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      #panel.background = element_blank(),
      
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
      ,legend.justification=c(1,0), legend.position=c(0.97,0.70),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p2()
##/plot

#/dupication frequency as the function of interMH

##III.dupication frequency as the function of gene expression######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData\\SRR")
load("MN_Spombe_972h-_Rep1.MHfeatures.Rdata")
attach(MHfeatures)
len_gene_Dup <- data.frame(MHlen,entire_geneMedian,dup)
len_gene_Dup$group <- cut(len_gene_Dup[,"entire_geneMedian"], breaks = seq(from=0,to=15,by=0.1))
len_tab <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for (i in 1:3) {
  if(i+3!=6){
    ind <- MHlen==i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_gene_Dup[ind, ]))
  }
  if(i+3==6){
    ind <- MHlen>=i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_gene_Dup[ind, ]))
  }
  
}
num_gene_Dup <- table(len_tab[,c("label","group","dup")])

##generate entire_geneMedian and frequency data frame
gene_frame <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for(i in 1:3){
  label<- rep(i, each=dim(num_gene_Dup[i,,])[1])
  gene_median <- seq(from=0.5,to=15,by=0.1)
  fre <- num_gene_Dup[i,,2]/(num_gene_Dup[i,,1]+num_gene_Dup[i,,2])*100
  gene_frame <- rbind(gene_frame,cbind(label,gene_median,fre,num_gene_Dup[i,,2],num_gene_Dup[i,,1]))
}
##/generate entire_geneMedian and frequency data frame


##plot
library(ggplot2)
p3 <- function(){
  ggplot(data = gene_frame, aes(x = gene_median, y = fre))+ggtitle("Gene expression-Dup(haploid)")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=gene_frame,aes( colour = paste0("MHlen=",label+3)),method='loess',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="gene expression",y="% of MHPs with an MTD",fill = "MH length (bp)")+
    guides(color = guide_legend(reverse=T))+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(.0001,.001,.01,.1,1),labels=c(.0001,.001,.01,.1,1)
                  ,expand = c(0.2,0.00001)
                  #breaks = trans_breaks("log10", function(x) 10^x),
                  #labels = trans_format("log10", math_format(10^.x))
    ) +
    #ylim(0, 0.5)+
    theme_bw() +
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=18),
      axis.title.y = element_text(color="black", size=18),
      #delete background
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
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
      ,legend.justification=c(1,0), legend.position=c(0.33,0.75),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p3()
##/plot
##/dupication frequency as the function of gene expression

