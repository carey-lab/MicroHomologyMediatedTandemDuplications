###I.DATA PREPARE#######
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



###II.inter_MH ######
inter_Dup <- data.frame(interMH,dup)
ind <- MHlen>=6 # select the MH length >= 6
inter_Dup <- inter_Dup[ind, ]
num_inter_Dup <- table(inter_Dup)

x.name <- as.numeric(row.names(num_inter_Dup))
y <- num_inter_Dup[,2]/num_inter_Dup[,1]
ba <- plot(x.name,y,type = "l",pch = 20,col = 1, ylab = "Duplication frequency",xlab = "nt between MH repeat",main = "Duplication frequency(MH_len>=6)")

###III.monotanic model#####
x.name <- as.numeric(row.names(num_inter_Dup))
y <- num_inter_Dup[,2]/num_inter_Dup[,1]
which(y==max(y))
constant <- x.name[which(y==max(y))]
x.name <- abs(x.name-constant)
ba <- plot(x.name,y,type = "l",pch = 20,col = 1, ylab = "Duplication frequency",xlab = "(nt between MH repeat) - interMH(which has most duplications)",main = "Duplication frequency(MH_len>=6)")
##/monotanic model



###IV.dupication frequency as the function of GC content ######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features_pro.txt",header = T,as.is = T,fill = T,row.names = 1)
detach(MHfeatures)
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
  ggplot(data = inter_frame, aes(x = inter, y = fre))+ggtitle("interMH-Duplication percent")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=inter_frame,aes( colour = paste0("MHlen=",label+3)),method='loess',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="nucleotides between MHR(bp)",y="Duplication percent(%)",fill = "MH length (bp)")+
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
      ,legend.justification=c(1,0), legend.position=c(0.97,0.80),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p2()
##/plot

##short interMH(0-50nt)
len_inter_Dup <-  data.frame(MHlen,interMH,dup)
len_shortInter_Dup <- len_inter_Dup[interMH<50,]
len_shortInter_Dup$group <- cut(len_shortInter_Dup[,"interMH"], breaks = seq(from=0,to=50,by=1))
len_tab <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for (i in 1:3) {
  if(i+3!=6){
    ind <- MHlen==i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_shortInter_Dup[ind, ]))
  }
  if(i+3==6){
    ind <- MHlen>=i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_shortInter_Dup[ind, ]))
  }
  
}
num_inter_Dup <- table(len_tab[,c("label","group","dup")])

##generate interMH and frequency data frame
inter_frame <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for(i in 1:3){
  label<- rep(i, each=dim(num_inter_Dup[i,,])[1])
  inter <- seq(from=1,to=50,by=1)
  fre <- num_inter_Dup[i,,2]/(num_inter_Dup[i,,1]+num_inter_Dup[i,,2])*100
  inter_frame <- rbind(inter_frame,cbind(label,inter,fre,num_inter_Dup[i,,2],num_inter_Dup[i,,1]))
}
##/generate interMH and frequency data frame

##plot
library(ggplot2)
p2 <- function(){
  ggplot(data = inter_frame, aes(x = inter, y = fre))+ggtitle("shortInterMH-Duplication percent")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=inter_frame,aes( colour = paste0("MHlen=",label+3)),method='loess',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    guides(color = guide_legend(reverse=T))+
    labs(x="nucleotides between MHR(bp)",y="% of MHPs with an MTD",fill = "MH length (bp)")+
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
      ,legend.justification=c(1,0), legend.position=c(0.97,0.85),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p2()
##/plot

##/short interMH(0-50nt)

#/dupication frequency as the function of GC content 
