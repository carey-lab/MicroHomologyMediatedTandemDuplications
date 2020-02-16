####II.LOGISTICS REFRESSION#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\Manuscript-todo\\processeddata")
load("MHfeatures.Rdata")

####(1)balanced data
library(ROCR)
library(caret)
library(pROC)
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

##(2)cross validation
set.seed(1)#get the randome number of the predata
folds <- createFolds(y=MHfeatures_te[,MHfeatures_te$dup],k=10)#k-fold cross 

crossval <- function(folds){  
  max=0
  num=0
  auc_value<-as.numeric()
  pred <-data.frame(obs= numeric(), pre= numeric(), stringsAsFactors=FALSE)
  for(i in 1:10){
    fold_test <- MHfeatures_te[folds[[i]],]#folds[[i]] for test group
    fold_train <- MHfeatures_te[-folds[[i]],]# 
    fold_pre <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+logentire_geneMedian+logntclosestMHR,
                    data = fold_train,family = "binomial")##logistics regression model
    print(summary(fold_pre))
    fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
    pred <- rbind(pred,matrix(c(as.numeric(fold_test[,1]),fold_predict),byrow = F,ncol = 2))
    auc_value<- append(auc_value,as.numeric(auc(as.numeric(fold_test[,1]),fold_predict)))#10 AUC
  }
  num<-which.max(auc_value)
  print(auc_value)
  max(auc_value)
  AUC <- auc(pred[,1],pred[,2])#1 AUC for all the 1/10 testing MHR samples
  print(AUC)
}

for(i in 1:5){
  MHfeatures_te <- predata(i)
  crossval(folds)
}#run 10-folds cross-validation for 5 times

##/crosss validation
##/balanced data


##(3)predict MHR duplication probability
log_prediction <- function(){
  MHfeatures_te <- predata(1)#randomely choose some non-events ,seed =1
  pre <- glm(formula = dup ~ MHlen + GCcon + f.interMH + entire_nucleMedian+logentire_geneMedian+logntclosestMHR,
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

##(4)rescale the prediction
logprediction <- as.data.frame(logprediction)
x <- logprediction$predict
x[is.na(x)] <- 0
num <- sum(logprediction$dup)## 6234:the number of 1s in observation 
ind <- tail(order(x),num)#set the samples with biggest values as duplication events
logprediction$predictbi <- 0 
logprediction$predictbi[ind] <- 1
##/rescale the prediction

features_logprediction <- cbind(MHfeatures,logprediction$predictbi,logprediction$predict)#combine MHfeatures  and binary and continous 
colnames(features_logprediction) <- c(colnames(MHfeatures),colnames(logprediction[,c(2,3)]))
#save(features_logprediction,file = "features_logprediction.Rdata")

##/LOGISTICS REFRESSION

