##I.DATA PREPARE####
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHlocation <- read.table("10k.sign.count.tsv.sorted.bed",header = F,as.is = T,fill = T)
colnames(MHlocation) <- c("chr","str1","end1","str2","end2","dupReads","collapseReads")
save(MHlocation,file = "MHlocation.Rdata")

pombe_chr <- read.table("pombe_n.fasta.fai",header = F,as.is = T,fill = T)
colnames(pombe_chr) <- c("chr","chrlen","col3","col4","col5")
save(pombe_chr,file = "pombe_chr.Rdata")

load("loc_obs_pre2k.RData")
logprediction <- read.table("logprediction.txt",header = T,as.is = T,fill = T,row.names = 1)

##rescale the prediction
x <- logprediction$predict
x[is.na(x)] <- 0
num <- sum(logprediction$dup)## 6234:the number of 1s in observation 
ind <- tail(order(x),num)
logprediction$predictbi <- 0 
logprediction$predictbi[ind] <- 1
##/rescale the prediction
##/DATA PREPARE
##II.NumMHRs MODEL####
##calculate duplication in different chr(window=2k) and #of MHRs in every window
locObsPre <- function(chr){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- logprediction[ind,]
  sumObs <- c()
  sumPre <- c()
  sumMHR <- c()
  sumFloat <- c()
  step <- 500 #step is 500bp
  window <- 2000 #window is 2kb
  chrLength <- pombe_chr[pombe_chr[,"chr"]==eval(chr),"chrlen"] #length of chromosome
  group1 <- seq(from=0,to=chrLength-500,by=500)
  group2 <- seq(from=2000,to=chrLength,by=500)
  group2 <- c(group2,chrLength)
  group1 <- group1[1:length(group2)]
  for (i in 1:length(group1)) {
    ind <- (chrlocation[,"str1"] >group1[i])&(chrlocation[,"end1"]<group2[i])
    sumObs <- c(sumObs,sum(na.omit(chrlogprediction[ind,"dup"])))
    sumPre <- c(sumPre,sum(na.omit(chrlogprediction[ind,"predictbi"])))
    sumMHR <- c(sumMHR,sum(ind))##calculate the number in every window
    sumFloat <- c(sumFloat,sum(na.omit(chrlogprediction[ind,"predict"])))
  }
  loc <- seq(from=(window/2),to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    #loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
    loc<-loc[1:length(group1)]
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumPre,sumMHR,sumFloat)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_preI2k <-locObsPre("I")
loc_obs_preII2k <-locObsPre("II")
loc_obs_preIII2k <-locObsPre("III")
save(loc_obs_preI2k, loc_obs_preII2k,loc_obs_preIII2k, file = "loc_obs_pre2k.RData")
load("loc_obs_pre2k.RData")

##window=1kb
locObsPre <- function(chr){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- features_logprediction[ind,]
  sumObs <- c()
  sumPre <- c()
  sumMHR <- c()
  sumFloat <- c()
  step <- 500 #step is 500bp
  window <- 1000 #window is 1kb
  chrLength <- pombe_chr[pombe_chr[,"chr"]==eval(chr),"chrlen"] #length of chromosome
  group1 <- seq(from=0,to=chrLength-500,by=500)
  group2 <- seq(from=1000,to=chrLength,by=500)
  group2 <- c(group2,chrLength)
  for (i in 1:length(group1)) {
    ind <- (chrlocation[,"str1"] >group1[i])&(chrlocation[,"end1"]<group2[i])
    sumObs <- c(sumObs,sum(na.omit(chrlogprediction[ind,"dup"])))
    sumPre <- c(sumPre,sum(na.omit(chrlogprediction[ind,"predictbi"])))
    sumMHR <- c(sumMHR,sum(ind))##calculate the number in every window
    sumFloat <- c(sumFloat,sum(na.omit(chrlogprediction[ind,"predict"])))
  }
  loc <- seq(from=500,to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumPre,sumMHR,sumFloat)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_preI1k <-locObsPre("I")
loc_obs_preII1k <-locObsPre("II")
loc_obs_preIII1k <-locObsPre("III")
loc_obs_pre1k <- data.frame(rbind(loc_obs_preI1k,loc_obs_preII1k,loc_obs_preIII1k))
loc_obs_pre1k$chr <- c(rep("chrI",dim(loc_obs_preI2k)[1]),rep("chrII",dim(loc_obs_preII2k)[1]),rep("chrIII",dim(loc_obs_preIII2k)[1]))
save(loc_obs_pre1k, file = "loc_obs_pre1k.RData")
write.table(loc_obs_pre1k,file = "loc_obs_pre1k.txt" ,quote = F,sep = "\t",row.names = F)
load("loc_obs_pre1k.RData")
##/window=1kb
##calculate duplication in different chr(window=2k) and #of MHRs in every window
##/NULL MODEL

##III.PERMUTATION TEST(SHUFFLING CALCULATION)####
load("loc_obs_pre1k.RData")
##(1)calculate average rank in models
library(dplyr)
x <- arrange(loc_obs_pre2k,desc(sumPre))
x$rank_Predicted_Full_Model = 1:nrow(x)
obs_ave_rank <-mean( x$rank_Predicted_Full_Model[x$sumObs>=10])#average rank of(logistic model):711.625
y <- arrange(loc_obs_pre2k,desc(sumMHR))
y$rank_NumMHRs_per_Window = 1:nrow(y)
NumMHRs_ave_rank <-mean( y$rank_NumMHRs_per_Window[y$sumObs>=10])#average rank of(NumMHRs model:MHR counts):5997.179
##/calculate average rank in models

##(2)shuffle 10000 times
shuffled_aver_rank <- c()
for (i in 1:10000) {
  shuff <- mean( x$rank_Predicted_Full_Model[ sample(1:nrow(x),sum(x$sumObs>=10)) ])
  shuffled_aver_rank <- c(shuffled_aver_rank,shuff)
  i=i+1
}
shuffled_aver_rank <-data.frame(c(shuffled_aver_rank,obs_ave_rank,NumMHRs_ave_rank))
colnames(shuffled_aver_rank) <- "aver_rank"
p <- mean(shuffled_aver_rank <=obs_ave_rank)#p value=9.998e-05
p
##/shuffle 10000 times

##(3)plot shuffle
library(ggplot2)
p5 <- function(dat){
  ggplot(data = dat, aes(x = aver_rank))+
    #ggtitle(paste0("Duplication pediction in ",chr))+
    geom_histogram(color="black", fill="white")+
    #geom_point(aes(x =sumPre , y =sumObs), size=2, alpha = 0.8) +
    #geom_point(aes(x = loc, y = sumObs,colour="sumObs"), size=2, alpha = 0.8) +
    #geom_smooth(data = dat, aes(x =sumPre , y =sumObs ),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="",y="count",fill = "MH length (bp)")+
    geom_text(x=2000, y=1900, label="p=9.998e-05")+
    geom_text(x=5800, y=300, label="simplest model")+
    geom_text(x=550, y=300, label="5-feature\n model")+
    theme_bw() +
    geom_segment(aes(x=550, xend=550, y=150, yend=20), 
                 arrow = arrow(length = unit(0.3, "cm")))+
    geom_segment(aes(x=6000, xend=6000, y=150, yend=20), 
                 arrow = arrow(length = unit(0.3, "cm")))+
    # scale_x_continuous(limits=c(0, loc[length(loc)]))+    
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
      ,legend.justification=c(1,0), legend.position=c(0.95,0.8),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p5(shuffled_aver_rank)
##/plot shuffle

##(4)plot enrichment ratio
getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
}
mode <- getmode(shuffled_aver_rank[,1])
ratio <- c(mode/obs_ave_rank,
           mode/NumMHRs_ave_rank)
df <- data.frame(model=c("5-feature model", "simplest model"),
                 ratio=ratio)
p6 <- function(){
  ggplot(df, aes(x=model, y=ratio, fill=model)) +
    geom_bar(stat="identity",width=0.5,show.legend = FALSE)+theme_minimal()+
    labs(x="",y="enrichment ratio:\n(predicted/random)-1")+
    geom_text(x=1, y=10, label="18",size=6)+
    geom_text(x=2, y=1, label="2",size=6)+
    
    theme_bw() +
    # scale_x_continuous(limits=c(0, loc[length(loc)]))+    
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
      ,legend.justification=c(1,0), legend.position=c(0.95,0.8),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p6()
##/plot enrichment ratio

##(5)plot duplication prediction in all chromosomes
loc_obs_pre2k$chr <- c(rep("chrI",dim(loc_obs_preI2k)[1]),rep("chrII",dim(loc_obs_preII2k)[1]),rep("chrIII",dim(loc_obs_preIII2k)[1]))
loc_obs_pre2k_p7 <- loc_obs_pre2k
loc_obs_pre2k_p7$sumPre[loc_obs_pre2k_p7$sumPre>=50]=50
ind <- tail(order(loc_obs_pre2k_p7$sumPre),56)#biggest num values
loc_obs_pre2k_p7_cir <- loc_obs_pre2k_p7[ind,]
p7 <- function(dat,chr){
    ggplot(data = dat, aes(x =sumPre , y =sumObs ))+ggtitle(paste0("Duplication prediction in ",chr,"\naccurancy=25/56 = 44.64%(sumObs>=10)"))+
    geom_abline(slope = 1,size=1, alpha = 0.8)+
    geom_segment(aes(x = 0, y = 9.5, xend = 51, yend = 9.5),color="orange", size=0.5)+
    geom_segment(aes(x = 0, y = 38, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = 0, y = 9.5, xend = 0, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = 51, y = 9.5, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_point(data= loc_obs_pre2k_p7_cir,aes(x =sumPre , y =sumObs),color = "orange" ,size=2.5, alpha = 0.8) +
    geom_point(aes(color = chr), size=1.5, alpha = 0.8) +
    
    #geom_smooth(data = dat, aes(color=chr),method='lm',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="prediction",y="observation",fill = "MH length (bp)")+
    #geom_text(x=x, y=y, label=text)+
        theme_bw() +
    # scale_x_continuous(limits=c(0, loc[length(loc)]))+    
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=18),
      axis.title.y = element_text(color="black", size=18),
      #delete background
      #panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
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
p7(loc_obs_pre2k_p7,"all chromosomes")
sum(loc_obs_pre2k_p7$sumObs>=10)

ind <- tail(order(loc_obs_pre2k_p7$sumPre),442)#biggest num values
loc_obs_pre2k_p7_cir5 <- loc_obs_pre2k_p7[ind,]
p7 <- function(dat,chr){
  ggplot(data = dat, aes(x =sumPre , y =sumObs ))+ggtitle(paste0("Duplication prediction in ",chr,"\naccurancy=73/442 = 16.5%"))+
    geom_abline(slope = 1,size=1, alpha = 0.8)+
    geom_segment(aes(x = -1, y = 4.5, xend = 51, yend = 4.5),color="orange", size=0.5)+
    geom_segment(aes(x = -1, y = 38, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = -1, y = 4.5, xend = -1, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = 51, y = 4.5, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_point(data= loc_obs_pre2k_p7_cir5,aes(x =sumPre , y =sumObs),color = "orange" ,size=2.5, alpha = 0.8) +
    geom_point(aes(color = chr), size=1.5, alpha = 0.8) +
    
    #geom_smooth(data = dat, aes(color=chr),method='lm',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="prediction",y="observation",fill = "MH length (bp)")+
    #geom_text(x=x, y=y, label=text)+
    theme_bw() +
    # scale_x_continuous(limits=c(0, loc[length(loc)]))+    
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=18),
      axis.title.y = element_text(color="black", size=18),
      #delete background
      #panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
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
p7(loc_obs_pre2k_p7,"all chromosomes")
sum(loc_obs_pre2k_p7_cir5$sumObs>=5)
sum(loc_obs_pre2k_p7$sumObs>=5)#442
##/plot duplication prediction in all chromosomes

##/PERMUTATION TEST(SHUFFLING CALCULATION)