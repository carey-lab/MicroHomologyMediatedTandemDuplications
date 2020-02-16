##II.dupication frequency as the function of nucleosome occupancy######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features_pro.txt",header = T,as.is = T,fill = T,row.names = 1)
attach(MHfeatures)
nucle_Dup <- data.frame(entire_nucleMedian,dup)
nucle_Dup[,1] <- round(nucle_Dup[,1]*10)/10
nucle_pre <- table(nucle_Dup)
nucle_pre[,2]/nucle_pre[,1]
x <- as.numeric(row.names(nucle_pre))
y <- nucle_pre[,2]/(nucle_pre[,1]+nucle_pre[,2])
ba <- plot(x,y,type = "l",pch = 20,col = 1, ylab = "Duplication frequency",xlab = "nucleosome occupancy",main = "Duplication frequency")

##bin nucle occu value#####
len_nucle_Dup <- data.frame(MHlen,entire_nucleMedian,dup)
len_nucle_Dup$group <- cut(len_nucle_Dup[,"entire_nucleMedian"], breaks = seq(from=1,to=50,by=0.1))
len_tab <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for (i in 1:3) {
  if(i+3!=6){
    ind <- MHlen==i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_nucle_Dup[ind, ]))
  }
  if(i+3==6){
    ind <- MHlen>=i+3 # select the MH length >= 6
    label <- rep(i,each=sum(ind))
    len_tab <- rbind(len_tab,cbind(label,len_nucle_Dup[ind, ]))
  }
  
}
num_nucle_Dup <- table(len_tab[,c("label","group","dup")])

##generate entire_nucleMedian and frequency data frame
nucle_frame <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for(i in 1:3){
  label<- rep(i, each=dim(num_nucle_Dup[i,,])[1])
  nucle_median <- seq(from=1,to=50,by=0.1)
  fre <- num_nucle_Dup[i,,2]/(num_nucle_Dup[i,,1]+num_nucle_Dup[i,,2])*100
  nucle_frame <- rbind(nucle_frame,cbind(label,nucle_median,fre,num_nucle_Dup[i,,2],num_nucle_Dup[i,,1]))
}
##/nuclerate entire_nucleMedian and frequency data frame


##plot
library(ggplot2)
p4 <- function(){
  ggplot(data = nucle_frame, aes(x = nucle_median, y = fre))+ggtitle("Nucleosome Occupancy-Duplication percent")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=nucle_frame,aes( colour = paste0("MHlen=",label+3)),method='loess',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="nucleosome occupancy",y="Duplication percent(%)",fill = "MH length (bp)")+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(.0001,.001,.01,.1,1,10),labels=c(.0001,.001,.01,.1,1,10)
                  ,expand = c(0.2,0.00001)
                  #breaks = trans_breaks("log10", function(x) 10^x),
                  #labels = trans_format("log10", math_format(10^.x))
    ) +
    #ylim(0, 0.5)+
    theme_bw() +
    theme(
      plot.title = element_text(lineheight=.8, size=18,hjust = 0.5),
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
      ,legend.justification=c(1,0), legend.position=c(0.91,0.02),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p4()
##/plot
##/bin nucle occu value



##sort bin nucle occu ####
##generate average GC content and frequency data
library(dplyr)
aveNuclepre <- function(){
  aveNucle_pre <- list()
  for (i in 4:5) {
  #MH length=4
  #i=4
  nucle_Dup <- data.frame(entire_nucleMedian,dup,MHlen)
  nucle_Dup <- nucle_Dup[nucle_Dup$MHlen == i,]
  nucle_Dup <- arrange(nucle_Dup,nucle_Dup[,"entire_nucleMedian"])
  nucle_Dup[,"entire_nucleMedian"] <- round(nucle_Dup[,"entire_nucleMedian"]*100)/100
  cons<- 100000
  lab <- c(rep(1:floor(dim(nucle_Dup)[1]/cons),each = cons),rep(floor(dim(nucle_Dup)[1]/cons)+1,each = dim(nucle_Dup)[1] - floor(dim(nucle_Dup)[1]/cons)*cons))
  nucle_Dup_lab <- cbind(nucle_Dup,lab)
  nucle_pre_tab <- table(nucle_Dup_lab[,c("entire_nucleMedian","dup","lab")])
  
  nucle_pre <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  aveNucle<- c()
  for(j in 1:dim(nucle_pre_tab)[3]){
    aveNucle <- c(aveNucle,median(nucle_Dup_lab[nucle_Dup_lab[,"lab"]==j,"entire_nucleMedian"]))
    nucle_pre <- rbind(nucle_pre,colSums(nucle_pre_tab[,,j]))
  }
  
  fre <- nucle_pre[,2]/(nucle_pre[,1]+nucle_pre[,2])
  nucle_pre <- cbind(nucle_pre,aveNucle,fre)
  #nucle_pre <- data.frame(matrix(unlist(nucle_pre), nrow=dim(nucle_pre)[1], byrow=F))
  aveNucle_pre[[i-3]] <- nucle_pre
  ##/MHlength = 4
  }

  #MH length>=7
  i=6
  nucle_Dup <- data.frame(entire_nucleMedian,dup,MHlen)
  nucle_Dup <- nucle_Dup[nucle_Dup$MHlen >= i,]
  nucle_Dup <- arrange(nucle_Dup,nucle_Dup[,"entire_nucleMedian"])
  nucle_Dup[,"entire_nucleMedian"] <- round(nucle_Dup[,"entire_nucleMedian"]*100)/100
  cons<- 100000
  lab <- c(rep(1:floor(dim(nucle_Dup)[1]/cons),each = cons),rep(floor(dim(nucle_Dup)[1]/cons)+1,each = dim(nucle_Dup)[1] - floor(dim(nucle_Dup)[1]/cons)*cons))
  nucle_Dup_lab <- cbind(nucle_Dup,lab)
  nucle_pre_tab <- table(nucle_Dup_lab[,c("entire_nucleMedian","dup","lab")])
  
  nucle_pre <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  aveNucle<- c()
  for(j in 1:dim(nucle_pre_tab)[3]){
    aveNucle <- c(aveNucle,median(nucle_Dup_lab[nucle_Dup_lab[,"lab"]==j,"entire_nucleMedian"]))##average of the nucleosome occupancy in every bin
    nucle_pre <- rbind(nucle_pre,colSums(nucle_pre_tab[,,j]))
  }
  
  fre <- nucle_pre[,2]/(nucle_pre[,1]+nucle_pre[,2])
  nucle_pre <- cbind(nucle_pre,aveNucle,fre)
  #nucle_pre <- data.frame(matrix(unlist(nucle_pre), nrow=dim(nucle_pre)[1], byrow=F))
  aveNucle_pre[[i-3]] <- nucle_pre
  ##/MH length >=7
  aveNucle_pre
}
aveNucle_pre <- aveNuclepre()
##/generate average GC content and frequency data

##generate aveNucle content and frequency data frame
aveNucle_frame <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
for(i in 1:3){
  label<- rep(i, each=dim(aveNucle_pre[[i]])[1])
  #lab_ave <- cbind(label,aveNucle_pre[[i]][,3:4])
  aveNucle_frame <- rbind(aveNucle_frame,cbind(label,aveNucle_pre[[i]][,3:4]))
}
##/generate aveNucle content and frequency data frame


##plot
library(ggplot2)
p4 <- function(){
  ggplot(data = aveNucle_frame, aes(x = aveNucle, y = fre*100))+ggtitle("Nucleosome Occupancy-Duplication percent")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=aveNucle_frame,aes( colour = paste0("MHlen=",label+3)),method='lm',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="median nucleosome occupancy value",y="Duplication percent (%)",fill = "MH length (bp)")+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(.0001,.001,.01,.1,1),labels=c(.0001,.001,.01,.1,1)
                  ,expand = c(0.14,0.0001)
                  #breaks = trans_breaks("log10", function(x) 10^x),
                  #labels = trans_format("log10", math_format(10^.x))
    ) +
    theme_bw() +
    
    theme(
      plot.title = element_text(lineheight=.8, size=18,hjust = 0.5),
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
p4()
##/plot

##/sort bin nucle occu 

##/dupication frequency as the function of nicleosome occupancy

