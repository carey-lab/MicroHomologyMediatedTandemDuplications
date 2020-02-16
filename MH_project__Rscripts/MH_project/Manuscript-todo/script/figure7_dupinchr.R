##I.DATA PREPARE####
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\Manuscript-todo\\processeddata")
load("MHlocation.Rdata")
load("pombe_chr.Rdata")
load("features_logprediction.Rdata")
load("loc_obs_pre1k.RData")
##/DATA PREPARE

##plot duplication prediction in all chromosomes
loc_obs_pre1k$sumPre[loc_obs_pre1k$sumPre>=50]=50
sum(loc_obs_pre1k$sumObs>=10)#25
sum((loc_obs_pre1k$sumObs>=10)&(loc_obs_pre1k$sumPre>=10))#9
ind <- tail(order(loc_obs_pre1k$sumPre),sum(loc_obs_pre1k$sumObs>=10))#biggest num values
loc_obs_pre1k_cir <- loc_obs_pre1k[ind,]
library(ggplot2)
p7 <- function(dat,chr){
  ggplot(data = dat, aes(x =sumPre , y =sumObs ))+ggtitle(paste0("Duplication prediction in ",chr,"\naccurancy=9/25 = 36%(sumObs>=10)"))+
    geom_abline(slope = 1,size=1, alpha = 0.8)+
    geom_segment(aes(x = 0, y = 9.5, xend = 51, yend = 9.5),color="orange", size=0.5)+
    geom_segment(aes(x = 0, y = 38, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = 0, y = 9.5, xend = 0, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = 51, y = 9.5, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_point(data= loc_obs_pre1k_cir,aes(x =sumPre , y =sumObs),color = "orange" ,size=2.5, alpha = 0.8) +
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
p7(loc_obs_pre1k,"all chromosomes")



sum(loc_obs_pre1k$sumObs>=5)#134
ind <- tail(order(loc_obs_pre1k$sumPre),sum(loc_obs_pre1k$sumObs>=5))#biggest num values
loc_obs_pre1k_cir5 <- loc_obs_pre1k[ind,]
sum(loc_obs_pre1k_cir5$sumObs>=5)#44
p7 <- function(dat,chr){
  ggplot(data = dat, aes(x =sumPre , y =sumObs ))+ggtitle(paste0("Duplication prediction in ",chr,"\naccurancy= 44/134 = 32.84%(sumObs>=5)"))+
    geom_abline(slope = 1,size=1, alpha = 0.8)+
    geom_segment(aes(x = -1, y = 4.5, xend = 51, yend = 4.5),color="orange", size=0.5)+
    geom_segment(aes(x = -1, y = 38, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = -1, y = 4.5, xend = -1, yend = 38),color="orange", size=0.5)+
    geom_segment(aes(x = 51, y = 4.5, xend = 51, yend = 38),color="orange", size=0.5)+
    geom_point(data= loc_obs_pre1k_cir5,aes(x =sumPre , y =sumObs),color = "orange" ,size=2.5, alpha = 0.8) +
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
p7(loc_obs_pre1k,"all chromosomes")

##/plot duplication prediction in all chromosomes
