###II.dupication frequency as the function of interMH ######
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\new_MH_project\\Processeddata")
load("MHfeatures.SZ=10.Rdata")
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
  ggplot(data = inter_frame, aes(x = inter, y = fre))+ggtitle("interMH-Duplication percent\nSZ=10")+
    geom_point(aes(color = paste0("MHlen=",label+3)), size=2, alpha = 0.8) +
    geom_smooth(data=inter_frame,aes( colour = paste0("MHlen=",label+3)),method='loess',formula=y~x,se=T, fullrange=FALSE, level=0.95, linetype = "solid") +
    labs(x="NT between MH seqs",y="% of MHPs with an MTD",fill = "MH length (bp)")+
    guides(color = guide_legend(reverse=T))+
    #scale_y_log10(breaks=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1),labels=c(.001,.002,.003,.004,.005,.01,.02,.03,.04,.05,.1,1)
    scale_y_log10(breaks=c(.0001,.001,.01,.1,1),labels=c(.0001,.001,.01,.1,1)
                  ,expand = c(0.2,0.0001)
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
p2()
##/plot

#/dupication frequency as the function of interMH
