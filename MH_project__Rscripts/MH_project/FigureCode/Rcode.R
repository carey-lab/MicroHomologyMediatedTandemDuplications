setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHlocation <- read.table("10k.sign.count.tsv.sorted.bed",header = F,as.is = T,fill = T)
colnames(MHlocation) <- c("chr","str1","end1","str2","end2","dupReads","collapseReads")

pombe_chr <- read.table("pombe_n.fasta.fai",header = F,as.is = T,fill = T)
colnames(pombe_chr) <- c("chr","chrlen","col3","col4","col5")

load("loc_obs_pre2k.RData")
logprediction <- read.table("logprediction.txt",header = T,as.is = T,fill = T,row.names = 1)

##rescale the prediction
x <- logprediction$predict
x[is.na(x)] <- 0
num <- sum(logprediction$dup)## 6234:the number of 1s in observation 
ind <- tail(order(x),num)#biggest num values
logprediction$predictbi <- 0 
logprediction$predictbi[ind] <- 1
##/rescale the prediction
##null model
set.seed(1)
ind <- sample(1:length(logprediction$predict),num)
logprediction$null <- 0 
logprediction$null[ind] <-1
##/null model


##calculate duplication in different chr(window=2k) in null model
load("loc_obs_null.RData")
locObsPre <- function(chr){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- logprediction[ind,]
  sumObs<- c()
  sumNull<- c()
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
    sumNull <- c(sumNull,sum(na.omit(chrlogprediction[ind,"null"])))
  }
  loc <- seq(from=(window/2),to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    #loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
    loc<-loc[1:length(group1)]
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumNull)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_nullI2k <-locObsPre("I")
loc_obs_nullII2k <-locObsPre("II")
loc_obs_nullIII2k <-locObsPre("III")
save(loc_obs_nullI2k, loc_obs_nullII2k,loc_obs_nullIII2k, file = "loc_obs_null2k.RData")
load("loc_obs_null2k.RData")
##/calculate duplication in different chr(window=2k) in null model

library(ggplot2)
loc_obs_pre2k <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
loc_obs_pre2k <- data.frame(rbind(loc_obs_preI2k,loc_obs_preII2k,loc_obs_preIII2k))
loc_obs_pre2k$chr <- c(rep("chrI",dim(loc_obs_preI2k)[1]),rep("chrII",dim(loc_obs_preII2k)[1]),rep("chrIII",dim(loc_obs_preIII2k)[1]))
p7 <- function(dat,chr){
  ggplot(data = dat, aes(x =sumPre , y =sumObs ))+ggtitle(paste0("Duplication pediction in ",chr))+
    geom_point(aes(color = chr), size=1.5, alpha = 0.8) +
    #geom_smooth(data = dat, aes(color=chr),method='lm',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="predction",y="observation",fill = "MH length (bp)")+
    #geom_text(x=x, y=y, label=text)+
    geom_abline(slope = 1,size=1.3, alpha = 0.8)+
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
      ,legend.justification=c(1,0), legend.position=c(0.95,0.05),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p7(loc_obs_pre2k,"all chromosomes")


##null
loc_obs_null2k <- data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
loc_obs_null2k <- data.frame(rbind(loc_obs_nullI2k,loc_obs_nullII2k,loc_obs_nullIII2k))
y <- arrange(loc_obs_null2k,desc(sumNull))
y$rank = 1:nrow(y)
null_ave_rank <-mean( y$rank[y$sumObs>=10])#13074.84
##/null
x <- arrange(loc_obs_pre2k,desc(sumPre))
x$rank = 1:nrow(x)
obs_ave_rank <-mean( x$rank[x$sumObs>=10])#711.625
shuffled_aver_rank <- c()
for (i in 1:10000) {
  shuff <- mean( x$rank[ sample(1:nrow(x),sum(x$sumObs>=10)) ])
  shuffled_aver_rank <- c(shuffled_aver_rank,shuff)
  i=i+1
}#shuffle 10000 times

p <- mean(shuffled_aver_rank <=obs_ave_rank)#9.998e-05
shuffled_aver_rank <-data.frame(c(shuffled_aver_rank,obs_ave_rank,null_ave_rank))
colnames(shuffled_aver_rank) <- "aver_rank"
library(ggplot2)
p5 <- function(dat){
  ggplot(data = dat, aes(x = aver_rank))+
    #ggtitle(paste0("Duplication pediction in ",chr))+
    geom_histogram(color="black", fill="white")+
    #geom_point(aes(x =sumPre , y =sumObs), size=2, alpha = 0.8) +
    #geom_point(aes(x = loc, y = sumObs,colour="sumObs"), size=2, alpha = 0.8) +
    #geom_smooth(data = dat, aes(x =sumPre , y =sumObs ),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="",y="count",fill = "MH length (bp)")+
    geom_text(x=2000, y=1900, label="p=9.999e-05")+
    geom_text(x=800, y=300, label="measured")+
    geom_text(x=13000, y=300, label="null model")+
    theme_bw() +
    geom_segment(aes(x=550, xend=550, y=150, yend=20), 
                 arrow = arrow(length = unit(0.3, "cm")))+
    geom_segment(aes(x=13000, xend=13000, y=150, yend=20), 
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

