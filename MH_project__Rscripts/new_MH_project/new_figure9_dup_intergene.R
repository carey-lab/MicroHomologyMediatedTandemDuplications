setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\new_MH_project\\ProcessedData")
containedInFeature <- read.table("10k.SZ=10.sign.count.tsv.containedInFeature.txt",header = F,as.is = T)
colnames(containedInFeature) <- c("MHlen","interMH","duplication","interGene","geneName")
containedInFeature$dup <- containedInFeature$duplication
containedInFeature$dup[containedInFeature$dup!=0]=1 #binary dupication: turn 1,2,3 to 0   
containedInFeature$duplen <- containedInFeature$MHlen+containedInFeature$interMH
save(containedInFeature,file = "containedInFeatureSZ=10.Rdata")

  
#####intergenic and within gene MHPs#######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\new_MH_project\\ProcessedData")
load("containedInFeatureSZ=10.Rdata")
attach(containedInFeature)
a <-  containedInFeature[interGene=="intergenic",c("MHlen","duplen")]
a[as.numeric(a$MHlen)>=7,"MHlen"]=7
b <- as.data.frame(table(a))
len_inter_Freq_intergene <- b[as.numeric(b$duplen)<=150&as.numeric(b$duplen)>=137,]

a <-  containedInFeature[interGene=="gene",c("MHlen","duplen")]
a[as.numeric(a$MHlen)>=7,"MHlen"]=7
b <- as.data.frame(table(a))
len_inter_Freq_gene <- b[as.numeric(b$duplen)<=150&as.numeric(b$duplen)>=137,]

##plot
library(ggplot2)
p9 <- function(){
  ggplot(data = len_inter_Freq_intergene, aes(x = duplen, y = Freq,group=MHlen))+ggtitle("intergenic MHPs")+
    geom_point(aes(color = paste0("MHlen=",MHlen)), size=2, alpha = 0.8) +
    geom_line(aes(color = paste0("MHlen=",MHlen)))+
    labs(x="BP that would be duplicated by an MTD",y="# of MH pairs",fill = "MH length (bp)")+
    #guides(color = guide_legend(reverse=T))+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(1,100,1000,10000,100000),labels=c(1,100,1000,10000,100000)
                  ,expand = c(0.1,0.4)
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
      ,legend.justification=c(1,0), legend.position=c(0.97,0.80),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p9()
p10 <- function(){
  ggplot(data = len_inter_Freq_gene, aes(x = duplen, y = Freq, group = MHlen))+ggtitle("within gene MHPs")+
    geom_point(aes(color = paste0("MHlen=",MHlen)), size=2, alpha = 0.8) +
    geom_line(aes(color = paste0("MHlen=",MHlen)))+
    #geom_smooth(data=len_inter_Freq_gene,aes(colour = paste0("MHlen=",MHlen)),formula=y~x,se=F, fullrange=T, level=0.95, linetype = "solid") +
    labs(x="BP that would be duplicated by an MTD",y="# of MH pairs",fill = "MH length (bp)")+
    #guides(color = guide_legend(reverse=T))+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(1,100,1000,10000,100000),labels=c(1,100,1000,10000,100000)
                  ,expand = c(0.1,0.4)
                  #breaks = trans_breaks("log10", function(x) 10^x),
                  #labels = trans_format("log10", math_format(10^.x))
    ) +
    #scale_x_continuous(breaks = seq(145,156,3),labels = seq(145,156,3))+
    #ylim(0, 0.5)+
    theme_bw() +
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
      ,legend.justification=c(1,0), legend.position=c(0.97,0.80),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p10()
##/plot
