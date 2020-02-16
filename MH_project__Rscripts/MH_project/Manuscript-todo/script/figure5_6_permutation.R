##I.DATA PREPARE####
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\Manuscript-todo\\processeddata")
load("MHlocation.Rdata")
load("pombe_chr.Rdata")
load("features_logprediction.Rdata")
load("loc_obs_pre1k.RData")
##/DATA PREPARE

##II.NumMHRs MODEL####
##calculate duplication in different chr(window=2k) and #of MHRs in every window
load("loc_obs_pre1k.RData")
##calculate duplication in different chr(window=2k) and #of MHRs in every window
##/NULL MODEL

##III.PERMUTATION TEST(SHUFFLING CALCULATION)####
##(1)calculate average rank in models
library(dplyr)
x <- arrange(loc_obs_pre1k,desc(sumPre))
x$rank_Predicted_Full_Model = 1:nrow(x)
obs_ave_rank <- mean( x$rank_Predicted_Full_Model[x$sumObs>=10])#average rank of(logistic model):711.625
y <- arrange(loc_obs_pre1k,desc(sumMHR))
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
NumMHRs_ave_rank
obs_ave_rank
p5 <- function(dat){
  ggplot(data = dat, aes(x = aver_rank))+
    #ggtitle(paste0("Duplication pediction in ",chr))+
    geom_histogram(color="black", fill="white")+
    #geom_point(aes(x =sumPre , y =sumObs), size=2, alpha = 0.8) +
    #geom_point(aes(x = loc, y = sumObs,colour="sumObs"), size=2, alpha = 0.8) +
    #geom_smooth(data = dat, aes(x =sumPre , y =sumObs ),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="",y="count",fill = "MH length (bp)")+
    geom_text(x=2000, y=1900, label="p=9.998e-05",)+
    geom_text(x=NumMHRs_ave_rank, y=250, label="NumMHRs\n model",size=3.5)+
    geom_text(x=obs_ave_rank, y=250, label="6-feature\n model",size=3.5)+
    theme_bw() +
    geom_segment(aes(x=obs_ave_rank, xend=obs_ave_rank, y=150, yend=20), 
                 arrow = arrow(length = unit(0.3, "cm")))+
    geom_segment(aes(x=NumMHRs_ave_rank, xend=NumMHRs_ave_rank, y=150, yend=20), 
                 arrow = arrow(length = unit(0.3, "cm")))+
    # scale_x_continuous(limits=c(0, loc[length(loc)]))+    
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=18),
      axis.title.y = element_text(color="black", size=18),
      #delete background
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
df <- data.frame(model=c("6-feature model", "numMHRs model"),
                 ratio=ratio)
p6 <- function(){
  ggplot(df, aes(x=model, y=ratio, fill=model)) +
    geom_bar(stat="identity",width=0.5,show.legend = FALSE)+theme_minimal()+
    labs(x="",y="enrichment ratio:\n(predicted/random)-1")+
    geom_text(x=1, y=round(ratio[1])/2, label=round(ratio[1]),size=6)+
    geom_text(x=2, y=round(ratio[2])/2, label=round(ratio[2]),size=6)+
    
    theme_bw() +
    # scale_x_continuous(limits=c(0, loc[length(loc)]))+    
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=18),
      axis.title.y = element_text(color="black", size=18),
      #delete background
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
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

##/PERMUTATION TEST(SHUFFLING CALCULATION)