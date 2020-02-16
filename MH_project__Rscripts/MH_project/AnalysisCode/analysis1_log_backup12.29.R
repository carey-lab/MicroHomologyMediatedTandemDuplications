####logistics regression#### 
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
#MHfeatures <- read.table("10k.sign.count.tsv.GSM1023946_MN_Mit.features.txt",header = F,as.is = T)
MHfeatures <- read.table("10k.sign.count.tsv.MN_Spombe_972h-_Rep1.features.txt",header = F,as.is = T)
#MHfeatures <- read.table("10k.sign.count.tsv.Nucleosome-density-wtNucWave-reps-median.depth_wl_trimmed_PE2.features.txt",header = F,as.is = T)

colnames(MHfeatures) <- c("duplication","MHlen","interMH","MHseq","entire_nucleSum","entire_nucleMean","entire_nucleMedian","entire_nucleMin","entire_nucleMax","entire_geneSum","entire_geneMean","entire_geneMedian","entire_geneMin","entire_geneMax") 
attach(MHfeatures)
MHfeatures$dup <- duplication
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
#/log2(gene expression)

attach(MHfeatures)
#detach(MHfeatures)

####balanced data
library(ROCR)
event <- MHfeatures[MHfeatures$dup==1,]
log <- function(x,feature){
  set.seed(1)
  nonevent <- MHfeatures[sample(which(MHfeatures$dup==0),6000),]
  MHfeatures_te <- data.frame(rbind(event,nonevent))
  attach(MHfeatures_te)
  dup<- MHfeatures_te$dup
  MHfeatures_te <- data.frame(scale(subset(MHfeatures_te,select=-c(dup,MHseq))))#MHseq:MH sequence
  MHfeatures_te <- data.frame(dup,MHfeatures_te)
  colnames(MHfeatures_te) 
  attach(MHfeatures_te)
  
  ##logstic model
  g.out <- glm(formula = dup ~ MHlen + f.interMH + GCcon + entire_nucleMedian + logentire_geneMedian, data = MHfeatures_te, family = "binomial")
  g.out <- glm(formula = dup ~ MHlen + f.interMH + GCcon + entire_nucleMedian + eval(feature) , data = MHfeatures_te, family = "binomial")
  summary(g.out)
  
  glm_predict <- predict(g.out,newdata=MHfeatures_te,type='response')
  p <- glm_predict
  pr <- prediction(p, dup)
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  list(summary(g.out),auc)
  ##/logstic model
}

##cross validation
event <- MHfeatures[MHfeatures$dup==1,]
set.seed(1)
nonevent <- MHfeatures[sample(which(MHfeatures$dup==0),6000),]
MHfeatures_te <- data.frame(rbind(event,nonevent))
attach(MHfeatures_te)
dup<- as.factor(MHfeatures_te$dup)
MHfeatures_te <- data.frame(scale(subset(MHfeatures_te,select=-c(dup,MHseq))))#MHseq:MH sequence
MHfeatures_te <- data.frame(dup,MHfeatures_te)
colnames(MHfeatures_te) 
attach(MHfeatures_te)

library(boot)
cost<- function(r, pi = 0) mean(abs(r-pi) > 0.5)
cv.err <- cv.glm(MHfeatures_te, g.out,K=10,cost=cost)

library(caret)
library(mlbench)
# define training control
train_control <- trainControl(method = "cv", number = 10,summaryFunction = twoClassSummary,
                              classProbs = T)

# train the model on training set, Fit Naive Bayes Model
a <- data.frame(MHfeatures_te[,c(2,1)]) 
a[,2]<-as.factor(a[,2])
colnames(a) <- c("X0","Class")
a <- as.list(a)
model <- train(Class ~ .,
               data = a,
               trControl = train_control,
               method = "glm",
               family=binomial(),
               na.action = na.pass)
# print cv scores
summary(model)

##/crosss validation

c(log(1)[[2]],log(2)[[2]],log(3)[[2]],log(4)[[2]],log(5)[[2]])
set.seed(5)
nonevent <- MHfeatures[sample(which(MHfeatures$dup==0),6000),]
MHfeatures_te <- data.frame(rbind(event,nonevent))
colnames(MHfeatures_te) 
attach(MHfeatures_te)

##logstic model
g.out <- glm(formula = dup ~ MHlen + GCcon + interMH , data = MHfeatures_te, family = "binomial")
summary(g.out)

glm_predict <- predict(g.out,newdata=MHfeatures_te[,c("MHlen", "GCcon", "interMH" )],type='response')
p <- glm_predict
pr <- prediction(p, dup)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc





conf <- t(table(ifelse(glm_predict > 0.5,1,0),"dup"=dup)) 
conf
TN <- conf[1,1]
FP <- conf[1,2]
FN <- conf[2,1]
TP <- conf[2,2]
MCC = (TP * TN - FP * FN)/sqrt((TP+FP) * (TP+FN) * (FP+TN) * (TN+FN))

library(ROCR)
p <- glm_predict
pr <- prediction(p, dup)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf,main="balanced data with f(inter_MHlen)")

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
####/balanced data



####relogit######
library(Zelig)
z.out <- zelig(formula = dup ~ MHlen + GCcon, data = MHlen_dup_GC, model = "relogit", cite = FALSE)
x.out <- setx(z.out)
s.out <- sim(z.out, x = x.out)

zel_predict <- predict(z.out,newdata=MHlen_dup_GC[,2:3],type='response')
zel.results <- ifelse(zel_predict[[1]] > 0.5,1,0)

conf <- t(table(ifelse(zel_predict[[1]] > 0.0002,1,0),dup)) 
conf

MCC = (TP * TN - FP * FN)/sqrt((TP+FP) * (TP+FN) * (FP+TN) * (TN+FN))

##roc
library(ROCR)
p <- zel_predict

pr <- prediction(p, dup)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

conf <- t(table(ifelse(zel_predict > 0.0002,1,0),dup)) 
conf

MCC = (TP * TN - FP * FN)/sqrt((TP+FP) * (TP+FN) * (FP+TN) * (TN+FN))


