####I. DATA PREPARE#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
#MHfeatures <- read.table("10k.sign.count.tsv.GSM1023946_MN_Mit.features.txt",header = F,as.is = T)
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features.txt",header = F,as.is = T)
#MHfeatures <- read.table("10k.sign.count.tsv.Nucleosome-density-wtNucWave-reps-median.depth_wl_trimmed_PE2.features.txt",header = F,as.is = T)

withingene <- read.table("10k.sign.count.tsv.sorted_full_MHR__interect100pct__intersectANY.txt",header = F,as.is = T,fill = T)
MHfeatures$withingenes <- c(0)
MHfeatures$withingenes[(withingene[,6]==0)&(withingene[,7]!=0)]=1
MHfeatures$withingenes[(withingene[,6]!=0)&(withingene[,7]!=0)]=2

closest <- read.table("distance_to_closest_MHR_with_dup.bed",header = F,as.is = T,fill = T)
MHfeatures$ntclosestMHR <- closest[,5]
MHfeatures$logntclosestMHR <- log10(closest[,5])
MHfeatures$bintclosestMHR <- 0
MHfeatures$bintclosestMHR[MHfeatures$ntclosestMHR<=100] <- 1

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
#save(MHfeatures,file = "MHfeatures.Rdata")

####II.LOGISTICS REFRESSION#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
load("MHfeatures.Rdata")
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features_pro.txt",header = T,as.is = T,fill = T,row.names = 1)
attach(MHfeatures)
#write.table(MHfeatures,file = "10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features_pro.txt",quote = F,sep = "\t")
#save(MHfeatures,file = "MHfeatures.Rdata")

####balanced data
library(ROCR)
library(caret)
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
MHfeatures_te <- predata(1)#randomely choose some non-events ,seed =1

#cross validation
set.seed(1)#get the randome number of the predata
folds <- createFolds(y=MHfeatures_te[,MHfeatures_te$dup],k=10)#k-fold cross 
library(pROC)
crossval <- function(folds){  
  max=0
  num=0
  auc_value<-as.numeric()
  pred <-data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for(i in 1:10){
    fold_test <- MHfeatures_te[folds[[i]],]#folds[[i]] for test group
    fold_train <- MHfeatures_te[-folds[[i]],]# 
    #fold_pre <- glm(formula = dup ~ withingenes,
                    
    fold_pre <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+logentire_geneMedian+bintclosestMHR,
                    data = fold_train,family = "binomial")
    print(summary(fold_pre))
    fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
    pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),fold_predict),byrow = F,ncol = 2))
    auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),fold_predict)))
  }
  num<-which.max(auc_value)
  print(auc_value)
  max(auc_value)
  AUC <- auc(pred[,1],pred[,2])
  print(AUC)
}

for(i in 1:5){
  MHfeatures_te <- predata(i)
  crossval(folds)
}


##train the model(whose AUC is highest) and draw ROC
fold_pre1 <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+logentire_geneMedian+logntclosestMHR,
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
library(plotROC)

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
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
      # panel.grid.minor = element_blank(),
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

##/train the model(whose AUC is highest) and draw ROC
##/crosss validation
##/balanced data


##predict MHR duplication probability
log_prediction <- function(){
  MHfeatures_te <- predata(1)#randomely choose some non-events ,seed =1
  pre <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+entire_geneMedian,
             data = MHfeatures_te,family = "binomial")
  print(summary(pre))
  predict <- predict(pre,type='response',newdata=MHfeatures)
  predictbi <- ifelse(predict > 0.5,1,0)
  pred <- matrix(c(as.numeric(MHfeatures[,"dup"]),predictbi,predict),byrow = F,ncol = 3)
  colnames(pred) <- c("dup","predictbi","predict")
  AUC <- as.numeric(auc(as.numeric(MHfeatures[,"dup"]),predict))
  conf <- table(data.frame(pred[,c("dup","predictbi")]))
  pred
}
logprediction <- log_prediction()
##/predict MHR duplication probability

pre <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+logentire_geneMedian,
           data = MHfeatures_te,family = "binomial")
print(summary(pre))
predict <- predict(pre,type='response',newdata=MHfeatures_te)

pred <- matrix(c(as.numeric(MHfeatures_te[,"dup"])-1,predictbi,predict),byrow = F,ncol = 3)
colnames(pred) <- c("dup","predictbi","predict")
AUC <- as.numeric(auc(as.numeric(MHfeatures_te[,"dup"]),predict))
conf <- table(data.frame(pred[,c("dup","predictbi")]))
AUC
conf

#write.table(logprediction,file = "logprediction.txt",quote = F,sep = "\t")

##/LOGISTICS REFRESSION


####III.SVM(SUPPORT VECTOR MEACHINE)####
#data prepare
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
MHfeatures_te <- predata(1)
#/data prepare

library(e1071)
#svm cross validation
set.seed(1)
folds <- createFolds(y=MHfeatures_te[,MHfeatures_te$dup],k=10)#k-fold cross 
library(pROC)
crossval_svm <- function(folds){  
  max=0
  num=0
  auc_value<-as.numeric()
  pred <-data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for(i in 1:10){
    fold_test <- MHfeatures_te[folds[[i]],]#folds[[i]] for test group
    fold_train <- MHfeatures_te[-folds[[i]],]#
    fold_test <- na.omit(fold_test)
    fold_train <- na.omit(fold_train)
    fold_pre <- svm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian,
                    data = fold_train)
    #print(summary(fold_pre))
    fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
    pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),as.numeric(fold_predict)),byrow = F,ncol = 2))
    auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
  }
  num<-which.max(auc_value)
  print(auc_value)
  max(auc_value)
  AUC <- auc(pred[,1],pred[,2])
  print(AUC)
}

for(i in 1:5){
  MHfeatures_te <- predata(i)
  crossval_svm(folds)
}

###/SVM(SUPPORT VECTOR MEACHINE)

####IV.KNN(k-NearestNeighbor)####

library(caret)
fold_pre <- knn3(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian + logentire_geneMedian,
                 data = fold_train, k=5)
fold_predict <- predict(fold_pre,type='prob',newdata=fold_test)[,2]
auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))
pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),as.numeric(fold_predict)),byrow = F,ncol = 2))
auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))

#data prepare
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
MHfeatures_te <- predata(1)
#/data prepare

library(caret)
#KNN cross validation
set.seed(1)
folds <- createFolds(y=MHfeatures_te[,MHfeatures_te$dup],k=10)#k-fold cross 
library(pROC)
crossval_KNN <- function(folds){  
  max=0
  num=0
  auc_value<-as.numeric()
  pred <-data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for(i in 1:10){
    fold_test <- MHfeatures_te[folds[[i]],]#folds[[i]] for test group
    fold_train <- MHfeatures_te[-folds[[i]],]#
    fold_test <- na.omit(fold_test)
    fold_train <- na.omit(fold_train)
    fold_pre <- knn3(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian + logentire_geneMedian,
                    data = fold_train)
    #print(summary(fold_pre))
    fold_predict <- predict(fold_pre,type='prob',newdata=fold_test)[,2]
    pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),as.numeric(fold_predict)),byrow = F,ncol = 2))
    auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
  }
  num<-which.max(auc_value)
  print(auc_value)
  max(auc_value)
  AUC <- auc(pred[,1],pred[,2])
  print(AUC)
}

for(i in 1:5){
  MHfeatures_te <- predata(i)
  crossval_KNN(folds)
}

###/KNN(K-nearest neighbourhood)


####V.Neural Network####
# load the package
library(nnet)
data(iris)
# fit model
fit <- nnet(Species~., data=iris, size=4, decay=0.0001, maxit=500)
# summarize the fit
summary(fit)
# make predictions
predictions <- predict(fit, iris[,1:4], type="class")
# summarize accuracy
table(predictions, iris$Species)

#data prepare
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
MHfeatures_te <- predata(1)
#/data prepare

library(caret)
#KNN cross validation
set.seed(1)
folds <- createFolds(y=MHfeatures_te[,MHfeatures_te$dup],k=10)#k-fold cross 
library(pROC)
crossval_NN <- function(folds){  
  max=0
  num=0
  auc_value<-as.numeric()
  pred <-data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for(i in 1:10){
    fold_test <- MHfeatures_te[folds[[i]],]#folds[[i]] for test group
    fold_train <- MHfeatures_te[-folds[[i]],]#
    fold_test <- na.omit(fold_test)
    fold_train <- na.omit(fold_train)
    fold_pre <- nnet(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian + logentire_geneMedian,
                    data = fold_train,size=4,decay=0.0001, maxit=500)
    #print(summary(fold_pre))
    fold_predict <- predict(fold_pre,type='class',newdata=fold_test)
    pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),as.numeric(fold_predict)),byrow = F,ncol = 2))
    auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
  }
  num<-which.max(auc_value)
  print(auc_value)
  max(auc_value)
  AUC <- auc(pred[,1],pred[,2])
  print(AUC)
}

for(i in 1:5){
  MHfeatures_te <- predata(i)
  crossval_NN(folds)
}

###/Neural Network


####VI.Naive Bayes#####
#data prepare
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
MHfeatures_te <- predata(1)
#/data prepare

library(e1071)
#KNN cross validation
set.seed(1)
folds <- createFolds(y=MHfeatures_te[,MHfeatures_te$dup],k=10)#k-fold cross 
library(pROC)
crossval_naivebayes <- function(folds){  
  max=0
  num=0
  auc_value<-as.numeric()
  pred <-data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for(i in 1:10){
    fold_test <- MHfeatures_te[folds[[i]],]#folds[[i]] for test group
    fold_train <- MHfeatures_te[-folds[[i]],]#
    fold_test <- na.omit(fold_test)
    fold_train <- na.omit(fold_train)
    fold_pre <- naiveBayes(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian + logentire_geneMedian,
                           data = fold_train)
    #print(summary(fold_pre))
    fold_predict <- predict(fold_pre,type='class',newdata=fold_test)
    pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),as.numeric(fold_predict)),byrow = F,ncol = 2))
    auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
  }
  num<-which.max(auc_value)
  print(auc_value)
  max(auc_value)
  AUC <- auc(pred[,1],pred[,2])
  print(AUC)
}

for(i in 1:5){
  MHfeatures_te <- predata(i)
  crossval_naivebayes(folds)
}
###/naive bayes

####VII. MDA###
library(mda)
#KNN cross validation
set.seed(1)
folds <- createFolds(y=MHfeatures_te[,MHfeatures_te$dup],k=10)#k-fold cross 
library(pROC)
crossval_mda <- function(folds){  
  max=0
  num=0
  auc_value<-as.numeric()
  pred <-data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for(i in 1:10){
    fold_test <- MHfeatures_te[folds[[i]],]#folds[[i]] for test group
    fold_train <- MHfeatures_te[-folds[[i]],]#
    fold_test <- na.omit(fold_test)
    fold_train <- na.omit(fold_train)
    fold_pre <- mda(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian + logentire_geneMedian,
                           data = fold_train)
    #print(summary(fold_pre))
    fold_predict <- predict(fold_pre,type='class',newdata=fold_test)
    pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),as.numeric(fold_predict)),byrow = F,ncol = 2))
    auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),as.numeric(fold_predict))))
  }
  num<-which.max(auc_value)
  print(auc_value)
  max(auc_value)
  AUC <- auc(pred[,1],pred[,2])
  print(AUC)
}

for(i in 1:5){
  MHfeatures_te <- predata(i)
  crossval_mda(folds)
}

##/MDA(Mixture Discriminant Analysis)
