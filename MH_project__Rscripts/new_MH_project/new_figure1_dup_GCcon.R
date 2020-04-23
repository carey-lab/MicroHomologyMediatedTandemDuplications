##II.dupication frequency as the function of GC content ######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\new_MH_project\\Processeddata")
load("MHfeatures.SZ=10.Rdata")
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
  ggplot(data = aveGC_frame, aes(x = aveGC*100, y = fre*100))+ggtitle("GC content-Duplication percent\nSZ=10")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=aveGC_frame,aes( colour = paste0("MHlen=",label+3)),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="MH sequence %GC",y="% of MHPs with an MTD",fill = "MH length (bp)")+
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
      ,legend.justification=c(1,0), legend.position=c(0.95,0.05),
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


