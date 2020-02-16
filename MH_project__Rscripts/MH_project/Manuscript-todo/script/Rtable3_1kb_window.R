#Make, a table, in text file format, with each 2kb window, and the # of MHRs, the # of observed duplications, the predicted # of duplications, and the sum of the floating-point predictions

##I.DATA PREPARE####
setwd("E:\\work\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\Manuscript-todo\\processeddata")
load("MHlocation.Rdata")
load("pombe_chr.Rdata")
load("features_logprediction.Rdata")
##/DATA PREPARE
##II.NumMHRs MODEL####
##calculate duplication in different chr(window=1k) and #of MHRs in every window
##window=1kb
locObsPre <- function(chr){
  ind <- MHlocation[,"chr"]==eval(chr)
  chrlocation <- MHlocation[ind,c("str1","end1")]
  chrlogprediction <- features_logprediction[ind,]
  sumObs <- c()
  sumPre <- c()
  sumMHR <- c()
  sumFloat <- c()
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
    sumMHR <- c(sumMHR,sum(ind))##calculate the number in every window
    sumFloat <- c(sumFloat,sum(na.omit(chrlogprediction[ind,"predict"])))
  }
  loc <- seq(from=500,to=chrLength,by=500)
  if(length(loc)!=length(group2)){
    loc <- c(loc,(group1[length(group1)]+group2[length(group2)])/2)
  }
  loc_obs_pre <- data.frame(loc,sumObs,sumPre,sumMHR,sumFloat)#loc_obs_pre dataframe
  loc_obs_pre
}
loc_obs_preI1k <-locObsPre("I")
loc_obs_preII1k <-locObsPre("II")
loc_obs_preIII1k <-locObsPre("III")
loc_obs_pre1k <- data.frame(rbind(loc_obs_preI1k,loc_obs_preII1k,loc_obs_preIII1k))
loc_obs_pre1k$chr <- c(rep("chrI",dim(loc_obs_preI1k)[1]),rep("chrII",dim(loc_obs_preII1k)[1]),rep("chrIII",dim(loc_obs_preIII1k)[1]))
#save(loc_obs_pre1k, file = "loc_obs_pre1k.RData")
#write.table(loc_obs_pre1k,file = "loc_obs_pre1k.txt" ,quote = F,sep = "\t",row.names = F)
#load("loc_obs_pre1k.RData")
##/window=1kb
##calculate duplication in different chr(window=1k) and #of MHRs in every window
##/NULL MODEL
