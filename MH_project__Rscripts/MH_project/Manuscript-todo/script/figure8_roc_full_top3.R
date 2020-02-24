####II.LOGISTICS REFRESSION#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\Manuscript-todo\\processeddata")
load("MHfeatures.Rdata")

####(1)balanced data
library(ROCR)
library(caret)
library(pROC)
library(plotROC)
predata <- function(x){
  #prepare data
  event <- MHfeatures[MHfeatures$dup==1,]
  set.seed(x)
  nonevent <- MHfeatures[sample(which(MHfeatures$dup==0),6000),]
  MHfeatures_te <- data.frame(rbind(event,nonevent))
  attach(MHfeatures_te)
  dup<- as.factor(MHfeatures_te$dup)
  MHfeatures_te <- data.frame(scale(subset(MHfeatures_te,select=-c(dup,MHseq))))#MHseq:MHR sequence scale:normalise
  MHfeatures_te <- data.frame(dup,MHfeatures_te)
  colnames(MHfeatures_te) 
  detach(MHfeatures_te)
  MHfeatures_te
  #/prepare data
}
MHfeatures_te <- predata(1)#randomely choose some non-events,seed =1
##/balanced dat

##(2)train the model and draw ROC
fold_pre1 <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+logentire_geneMedian+logntclosestMHR
                 +interGCcon,
                 data = MHfeatures_te,family = "binomial")
fold_predict1 <- predict(fold_pre1,type='response',newdata=MHfeatures_te)
fold_pre2<- glm(formula = dup ~ MHlen + f.interMH + entire_nucleMedian,
                data = MHfeatures_te,family = "binomial")
fold_predict2 <- predict(fold_pre2,type='response',newdata=MHfeatures_te)
auc1 <- auc(as.numeric(MHfeatures_te$dup),fold_predict1)
auc2 <- auc(as.numeric(MHfeatures_te$dup),fold_predict2)
name <- rep(c("full","top3"),each=length(MHfeatures_te$dup))
test <- data.frame(rbind(cbind(as.numeric(MHfeatures_te$dup),fold_predict1),cbind(as.numeric(MHfeatures_te$dup),fold_predict2)),name)
colnames(test)<- c("dup","pre","name")                   

library(ggplot2)
p8 <- function(){
  ggplot(test, aes(d = dup, m = pre, color = name)) + 
    #geom_text(x= 0.75, y =0.5, label=paste0("AUC of full model is ",round(auc1,4),"\nAUC of top3 model is ",round(auc2,4)))+
    annotate("text", x = 0.65, y=0.5, size=6,label = paste0("AUC of full model is ",round(auc1,4),"\nAUC of top3 model is ",round(auc2,4)))+
    geom_roc(n.cuts = 0,size=1.3)+
    labs(y="Sensitivity", x="1-Specificity")+
    theme_bw() +
    theme(
      plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
      axis.title.x = element_text(color="black", size=22),
      axis.title.y = element_text(color="black", size=22),
      #delete background
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      
      #???????????????
      axis.line = element_line(colour = "black",size=0.5),
      #?????????
      axis.ticks = element_line(size=0.5),
      axis.ticks.length=unit(-0.2,"cm"),
      
      #x?????????
      axis.text.x = element_text(size=22,color='black',margin = margin(t=0.3, unit = "cm")),
      #y?????????
      axis.text.y = element_text(size=22,color='black',margin = margin(r=0.3, unit = "cm")),
      #??????
      #legend.title = element_text(colour="black", size=14),
      legend.title = element_blank(),
      #legend.text = element_blank()
      legend.text = element_text(colour=c("black"), size=16)
      # remove legend
      ,legend.justification=c(1,0), legend.position=c(0.95,0.75),
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")
      #,legend.position="none"
    ) 
}
p8()

##/train the model and draw ROC