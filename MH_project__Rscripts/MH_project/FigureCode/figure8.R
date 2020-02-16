#####I.predict probability of predition model####
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features_pro.txt",header = T,as.is = T,fill = T,row.names = 1)
attach(MHfeatures)
logprediction <- read.table("logprediction.txt",header = T,as.is = T,fill = T,row.names = 1)
dim(logprediction)
##rescale the prediction result
x <- logprediction$predict
x[is.na(x)] <- 0
num <- sum(logprediction$dup)## 6234:the number of 1s in observation 
ind <- tail(order(x),num)
logprediction$predictbi <- 0 
logprediction$predictbi[ind] <- 1
##/rescale the prediction result

MHlocation <- read.table("10k.sign.count.tsv.sorted.bed",header = F,as.is = T,fill = T)
colnames(MHlocation) <- c("chr","str1","end1","str2","end2","dupReads","collapseReads")

pombe_chr <- read.table("pombe_n.fasta.fai",header = F,as.is = T,fill = T)
colnames(pombe_chr) <- c("chr","chrlen","col3","col4","col5")

##calculate duplication in different chr(window=1k))
load("loc_obs_pre1k.RData")
locObsPre <- function(chr){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- logprediction[ind,]
  sumObs<- c()
  sumPre<- c()
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
  }
  loc <- seq(from=500,to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumPre)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_preI1k <-locObsPre("I")
loc_obs_preII1k <-locObsPre("II")
loc_obs_preIII1k <-locObsPre("III")
save(loc_obs_preI1k, loc_obs_preII1k,loc_obs_preIII1k, file = "loc_obs_pre1k.RData")
load("loc_obs_pre1k.RData")
cor<- matrix(0,ncol = 4,nrow = 3)
colnames(cor) <- c("1kPe","2kPe","1kSp","2kSp")
rownames(cor) <- c("chrI","chrII","chrIII")
cor[1,1] <- cor(loc_obs_preI1k[,"sumObs"],loc_obs_preI1k[,"sumPre"],method = c("pearson"))

cor[3,3] <- cor(loc_obs_preIII1k[,"sumObs"],loc_obs_preIII1k[,"sumPre"],method = c("spearman"))
##/calculate duplication in different chr

##calculate duplication in different chr(window=2k)
load("loc_obs_pre.RData")
locObsPre <- function(chr){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- logprediction[ind,]
  sumObs<- c()
  sumPre<- c()
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
  }
  loc <- seq(from=(window/2),to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    #loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
    loc<-loc[1:length(group1)]
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumPre)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_preI2k <-locObsPre("I")
loc_obs_preII2k <-locObsPre("II")
loc_obs_preIII2k <-locObsPre("III")
save(loc_obs_preI2k, loc_obs_preII2k,loc_obs_preIII2k, file = "loc_obs_pre2k.RData")
load("loc_obs_pre2k.RData")
cor[1,2] <- cor(loc_obs_preI2k[,"sumObs"],loc_obs_preI2k[,"sumPre"],method = c("pearson"))
cor[2,2] <- cor(loc_obs_preII2k[,"sumObs"],loc_obs_preII2k[,"sumPre"],method = c("pearson"))
cor[3,2] <- cor(loc_obs_preIII2k[,"sumObs"],loc_obs_preIII2k[,"sumPre"],method = c("pearson"))
cor[1,4] <- cor(loc_obs_preI2k[,"sumObs"],loc_obs_preI2k[,"sumPre"],method = c("spearman"))
cor[2,4] <- cor(loc_obs_preII2k[,"sumObs"],loc_obs_preII2k[,"sumPre"],method = c("spearman"))
cor[3,4] <- cor(loc_obs_preIII2k[,"sumObs"],loc_obs_preIII2k[,"sumPre"],method = c("spearman"))
cor
##/calculate duplication in different chr

##calculate duplication in different chr(window=3k)
locObsPre <- function(chr){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- logprediction[ind,]
  sumObs<- c()
  sumPre<- c()
  step <- 500 #step is 500bp
  window <- 3000 #window is 2kb
  chrLength <- pombe_chr[pombe_chr[,"chr"]==eval(chr),"chrlen"] #length of chromosome
  group1 <- seq(from=0,to=chrLength-500,by=500)
  group2 <- seq(from=window,to=chrLength,by=500)
  group2 <- c(group2,chrLength)
  group1 <- group1[1:length(group2)]
  for (i in 1:length(group1)) {
    #for (i in 1:10) { 
    ind <- (chrlocation[,"str1"] >group1[i])&(chrlocation[,"end1"]<group2[i])
    sumObs <- c(sumObs,sum(na.omit(chrlogprediction[ind,"dup"])))
    sumPre <- c(sumPre,sum(na.omit(chrlogprediction[ind,"predictbi"])))
  }
  loc <- seq(from=(window/2),to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    #loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
    loc<-loc[1:length(group1)]
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumPre)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_preI3k <-locObsPre("I")
loc_obs_preII3k <-locObsPre("II")
loc_obs_preIII3k <-locObsPre("III")
save(loc_obs_preI3k, loc_obs_preII3k,loc_obs_preIII3k, file = "loc_obs_pre3k.RData")
load("loc_obs_pre3k.RData")
##/calculate duplication in different chr(window=3k)

##calculate duplication in different chr(window=4k)
locObsPre <- function(chr,window){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- logprediction[ind,]
  sumObs<- c()
  sumPre<- c()
  step <- 500 #step is 500bp
  window <- window #window is 2kb
  chrLength <- pombe_chr[pombe_chr[,"chr"]==eval(chr),"chrlen"] #length of chromosome
  group1 <- seq(from=0,to=chrLength-500,by=500)
  group2 <- seq(from=window,to=chrLength,by=500)
  group2 <- c(group2,chrLength)
  group1 <- group1[1:length(group2)]
  for (i in 1:length(group1)) {
    #for (i in 1:10) { 
    ind <- (chrlocation[,"str1"] >group1[i])&(chrlocation[,"end1"]<group2[i])
    sumObs <- c(sumObs,sum(na.omit(chrlogprediction[ind,"dup"])))
    sumPre <- c(sumPre,sum(na.omit(chrlogprediction[ind,"predictbi"])))
  }
  loc <- seq(from=(window/2),to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    #loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
    loc<-loc[1:length(group1)]
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumPre)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_preI4k <-locObsPre("I",window = 4000)
loc_obs_preII4k <-locObsPre("II",window = 4000)
loc_obs_preIII4k <-locObsPre("III",window = 4000)
save(loc_obs_preI4k, loc_obs_preII4k,loc_obs_preIII4k, file = "loc_obs_pre4k.RData")
load("loc_obs_pre4k.RData")
##/calculate duplication in different chr(window=4k)


##plot
library(ggplot2)
p5 <- function(dat,chr,x,y,text){
  ggplot(data = dat, aes(x =sumPre , y =sumObs ))+ggtitle(paste0("Duplication pediction in ",chr))+
    geom_point(aes(x =sumPre , y =sumObs), size=2, alpha = 0.8) +
    #geom_point(aes(x = loc, y = sumObs,colour="sumObs"), size=2, alpha = 0.8) +
    #geom_smooth(data = dat, aes(x =sumPre , y =sumObs ),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="predction",y="observation",fill = "MH length (bp)")+
    geom_text(x=x, y=y, label=text)+
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
p5(loc_obs_preI1k,"chromosome I",x=10,y=6,paste0("cor(pearson)=",round(cor[1,1],4),"\ncor(spearman)=",round(cor[1,3],4)))
p5(loc_obs_preII1k,"chromosome II",x=5,y=20,paste0("cor(pearson)=",round(cor[2,1],4),"\ncor(spearman)=",round(cor[2,3],4)))
p5(loc_obs_preIII1k,"chromosome III",x=65,y=28,paste0("cor(pearson)=",round(cor[3,1],4),"\ncor(spearman)=",round(cor[3,3],4)))
p5(loc_obs_preI2k,"chromosome I",x=10,y=6,paste0("cor(pearson)=",round(cor[1,2],4),"\ncor(spearman)=",round(cor[1,4],4)))
p5(loc_obs_preII2k,"chromosome II",x=5,y=20,paste0("cor(pearson)=",round(cor[2,2],4),"\ncor(spearman)=",round(cor[2,4],4)))
p5(loc_obs_preIII2k,"chromosome III",x=65,y=28,paste0("cor(pearson)=",round(cor[3,2],4),"\ncor(spearman)=",round(cor[3,4],4)))

##/plot
##/calculate duplication in different chr

##hotspots
function(){
  step <- 500 #step is 500bp
  window <- 2000 #window is 2kb
  chrLength <- pombe_chr[pombe_chr[,"chr"]=="I","chrlen"] #length of chromosome
  group1 <- seq(from=0,to=chrLength-500,by=500)
  group2 <- seq(from=2000,to=chrLength,by=500)
  group2 <- c(group2,chrLength)
  group1 <- group1[1:length(group2)]
  ind1 <- loc_obs_preI2k[,"sumObs"]>=10
  ind2 <- loc_obs_preI2k[,"sumPre"]>=125 
  hotObs <- data.frame(cbind("I",group1[ind1],group2[ind1],loc_obs_preI2k[ind1,c("sumObs","sumPre")]))
  hotPre <- data.frame(cbind("I",group1[ind2],group2[ind2],loc_obs_preI2k[ind2,c("sumObs","sumPre")]))
  colnames(hotObs) <- c("chr","str","end","sumObs","sumPre")
  colnames(hotPre) <- c("chr","str","end","sumObs","sumPre")
  hotspotsObs <- hotObs
  hotspotsPre <- hotPre
  
  chrLength <- pombe_chr[pombe_chr[,"chr"]=="II","chrlen"] #length of chromosome
  group1 <- seq(from=0,to=chrLength-500,by=500)
  group2 <- seq(from=2000,to=chrLength,by=500)
  group2 <- c(group2,chrLength)
  group1 <- group1[1:length(group2)]
  ind1 <- loc_obs_preII2k[,"sumObs"]>=10
  ind2 <- loc_obs_preII2k[,"sumPre"]>=125 
  hotObs <- data.frame(cbind("II",group1[ind1],group2[ind1],loc_obs_preII2k[ind1,c("sumObs","sumPre")]))
  hotPre <- data.frame(cbind("II",group1[ind2],group2[ind2],loc_obs_preII2k[ind2,c("sumObs","sumPre")]))
  colnames(hotObs) <- c("chr","str","end","sumObs","sumPre")
  colnames(hotPre) <- c("chr","str","end","sumObs","sumPre")
  hotspotsObs <- data.frame(rbind(hotspotsObs,hotObs))
  hotspotsPre<- data.frame(rbind(hotspotsPre,hotPre))

  chrLength <- pombe_chr[pombe_chr[,"chr"]=="III","chrlen"] #length of chromosome
  group1 <- seq(from=0,to=chrLength-500,by=500)
  group2 <- seq(from=2000,to=chrLength,by=500)
  group2 <- c(group2,chrLength)
  group1 <- group1[1:length(group2)]
  group1 <- group1[1:(length(group1)-30)]
  group2 <- group2[1:(length(group2)-30)]
  ind1 <- loc_obs_preIII2k[,"sumObs"]>=10
  ind2 <- loc_obs_preIII2k[,"sumPre"]>=125 
  hotObs <- data.frame(cbind("III",group1[ind1],group2[ind1],loc_obs_preIII2k[ind1,c("sumObs","sumPre")]))
  hotPre <- data.frame(cbind("III",group1[ind2],group2[ind2],loc_obs_preIII2k[ind2,c("sumObs","sumPre")]))
  colnames(hotObs) <- c("chr","str","end","sumObs","sumPre")
  colnames(hotPre) <- c("chr","str","end","sumObs","sumPre")
  hotspotsObs <- data.frame(rbind(hotspotsObs,hotObs))
  hotspotsPre<- data.frame(rbind(hotspotsPre,hotPre))
}  
hotspotsPre
hotspotsObs
write.table(hotspotsPre,file="hotspotsPre.bed",quote = F,sep = "\t",col.names = F,row.names = F)
write.table(hotspotsObs,file="hotspotsObs.bed",quote = F,sep = "\t",col.names = F,row.names = F)
##/hotspots

##top 0.5% 
function(){
top_frame <- list()
x <- loc_obs_preI2k[,"sumObs"]
ind1 <- x > quantile(x, prob=0.995)
y <- loc_obs_preI2k[,"sumPre"]
ind2 <- y > quantile(y, prob=0.995)
top_frame[[1]] <- cbind(rep("other",each=length(x)),rep("other",each=length(y))) 
top_frame[[1]][ind1,1]="top"
top_frame[[1]][ind2,2]="top"
x <- loc_obs_preII2k[,"sumObs"]
ind1 <- x > quantile(x, prob=0.995)
y <- loc_obs_preII2k[,"sumPre"]
ind2 <- y > quantile(y, prob=0.995)
top_frame[[2]] <- cbind(rep("other",each=length(x)),rep("other",each=length(y))) 
top_frame[[2]][ind1,1]="top"
top_frame[[2]][ind2,2]="top"
x <- loc_obs_preIII2k[,"sumObs"]
ind1 <- x > quantile(x, prob=0.995)
y <- loc_obs_preIII2k[,"sumPre"]
ind2 <- y > quantile(y, prob=0.995)
top_frame[[3]] <- cbind(rep("other",each=length(x)),rep("other",each=length(y))) 
top_frame[[3]][ind1,1]="top"
top_frame[[3]][ind2,2]="top"

table(as.data.frame(top_frame[[1]][,c(1,2)]))#V1 is observed,V2 is predicted
fisher.test(table(as.data.frame(top_frame[[1]][,c(1,2)])))##fisher.test of top 0.5% perdiction and observation
table(as.data.frame(top_frame[[2]][,c(1,2)]))#V1 is observed,V2 is predicted
fisher.test(table(as.data.frame(top_frame[[2]][,c(1,2)])))##fisher.test of top 0.5% perdiction and observation
table(as.data.frame(top_frame[[3]][,c(1,2)]))#V1 is observed,V2 is predicted
fisher.test(table(as.data.frame(top_frame[[3]][,c(1,2)])))##fisher.test of top 0.5% perdiction and observation

##whole genome
whole <- rbind(loc_obs_preI2k[,c("sumObs","sumPre")],loc_obs_preII2k[,c("sumObs","sumPre")],loc_obs_preIII2k[,c("sumObs","sumPre")])
x <- whole[,"sumObs"]
ind1 <- x > quantile(x, prob=0.995)
y <- whole[,"sumPre"]
ind2 <- y > quantile(y, prob=0.995)
top_frame[[4]] <- cbind(rep("other",each=length(x)),rep("other",each=length(y))) 
top_frame[[4]][ind1,1]="top"
top_frame[[4]][ind2,2]="top"
table(as.data.frame(top_frame[[4]][,c(1,2)]))#V1 is observed,V2 is predicted
fisher.test(table(as.data.frame(top_frame[[4]][,c(1,2)])))##fisher.test of top 0.5% perdiction and observation
##/whole genome
}
table(as.data.frame(top_frame[[4]][,c(1,2)]))#V1 is observed,V2 is predicted
fisher.test(table(as.data.frame(top_frame[[4]][,c(1,2)])))##fisher.test of top 0.5% perdiction and observation
##/top 0.5% 

##/predict probability of predition model