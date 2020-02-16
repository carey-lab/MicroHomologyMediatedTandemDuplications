#2019.12.03 figure2
#!/apps/bioinfo/R-3.6.0share/bin/Rscipt
#source  /appsnew/source/R-3.6.0share.sh
###calculate log2(MH/kb-feature/wholegonome) & barplot

#table of different length of MH in 6 features(intron, CDS, 5-UTR...) and whole_genome, length from 4bp-12bp
args=commandArgs(T)
inter_feature <- read.table(args[1],header = F,as.is = T, fill =T)#bedtools intersect features(CDS,intron...) of each MH, file1 is pombe_n.intersect.txt
num <- table(inter_feature)

location <- read.table(args[2],header = F, comment.char = "#")#gff file2 is pombe_n.intersect
location$nt <- location$V3-location$V2+1 #V3 is start, V4 is stop 
total_nt <- as.data.frame(aggregate(location$nt,by=list(location$V4), sum)) #total nt of each feature
#total_nt#total nt of each feature

overlap <- intersect(row.names(num),as.character(total_nt$Group.1))
feature_length<- num[match(overlap,row.names(num)),]
feature_length#different length of MH in 6 features(intron, CDS, 5-UTR...) and whole_genome, length from 4bp-12bp
nt<- total_nt[match(overlap,as.character(total_nt$Group.1)),]
#nt


WGnt <- total_nt[[2]][match("chromosome",as.character(total_nt$Group.1))]#nt of whole_genome
featurent <- nt[-which(nt[,1]=="chromosome"),2]#nt of each feature


pdf(paste0(args[3],'.figure2.pdf')) #profiles_Spom/pombe.figure2

#figure of log2(MH/kb_whole_genome)
freWG <- feature_length[match("chromosome",as.character(total_nt$Group.1)),]*1000/WGnt
ba <- plot(4:(length(freWG)+3),log2(freWG),pch=20,type = "l",col = 1,ylab = "log2(#MH/kb_whole_genome)",xlab = "MH length(bp)",main = paste0(args[3],"_log2(#MH/kb_whole_genome)"))
ba <- plot(4:(length(freWG)+3),freWG,pch=20,type = "l",col = 1,ylab = "#MH/kb_whole_genome",xlab = "MH length(bp)",main = paste0(args[3],"_#MH/kb_whole_genome"))

#figure of log2(MH/kb_of_feature)
frefeature <- feature_length[-which(rownames(feature_length)=="chromosome"),]*1000/featurent
ba <- plot(4:(length(freWG)+3),frefeature[1,],type = "l",pch = 20,col = 1,ylim = c(0, 1.2*max(frefeature)), ylab = "#MH/kb_of_feature",xlab = "MH length(bp)",main = paste0(args[3],"_#MH/kb_of_feature"))
for(i in 2:dim(frefeature)[1])
{
	  lines(4:(length(freWG)+3),frefeature[i,],type = "o",pch = 20,cex=0.8,col = i,ylab = "#MH/kb_of_feature",xlab = "MH length(bp)")
}
legend("topright", legend=rownames(frefeature), col=c(1:dim(frefeature)[1]), lty=1)


ba <- plot(4:(length(freWG)+3),log2(frefeature[1,]),type = "l",pch = 20,col = 1,ylim = c(-2*max(log2(frefeature)),2*max(log2(frefeature))), ylab = "log2(#MH/kb_of_feature)",xlab = "MH length(bp)",main = paste0(args[3],"_log2(#MH/kb_of_feature)"))
for(i in 2:dim(frefeature)[1])
{
	  lines(4:(length(freWG)+3),log2(frefeature[i,]),type = "o",pch = 20,cex=0.8,col = i,ylab = "#MH/kb_of_feature",xlab = "MH length(bp)")
}
legend("topright", legend=rownames(frefeature), col=c(1:dim(frefeature)[1]), lty=1)


#figure of log2((#MH/kb_of_feature)/(#MH/kb_whole_genome))
frefea_WG <- log2(t(t(frefeature)/freWG))
ba <- plot(4:(length(freWG)+3),frefea_WG[1,],type = "l",col = 1,ylim= c(-2*max(frefea_WG),2*max(frefea_WG)),ylab = "log2(Fre-feature/Fre-WG)",xlab = "MH length(bp)",main = paste0(args[3],"feature-MHlength in features"))
for(i in 2:dim(frefeature)[1])
{
	  lines(4:(length(freWG)+3),frefea_WG[i,],type = "o",pch = 20,cex=0.8,col = i,ylab = "log2(Fre-feature/Fre-WG)",xlab = "MH length(bp)")
}
legend("bottomleft", legend=rownames(frefeature), col=c(1:dim(frefeature)[1]), lty=1)



WGnt <- total_nt[[2]][match("chromosome",as.character(total_nt$Group.1))]#nt of whole_genome
types <- c("intergenic","intron","CDS","five_prime_UTR","three_prime_UTR","ncRNA","nuclear_mt_pseudogene")

featurent <- na.omit(nt[match(types,as.character(nt$Group.1)),2])#nt of each feature

if(!is.null(featurent)){
	#figure of log2(#MH/kb_whole_genome)
	freWG <- feature_length[match("chromosome",as.character(total_nt$Group.1)),]*1000/WGnt
	ba <- plot(4:(length(freWG)+3),log2(freWG),pch=20,type = "l",col = 1,ylab = "log2(#MH/kb_whole_genome)",xlab = "MH length(bp)",main = paste0(args[3],"_log2(#MH/kb_whole_genome)"))
	ba <- plot(4:(length(freWG)+3),freWG,pch=20,type = "l",col = 1,ylab = "#MH/kb_whole_genome",xlab = "MH length(bp)",main = paste0(args[3],"_#MH/kb_whole_genome"))


	#figure of log2(#MH/kb_of_feature)
	frefeature <- feature_length[na.omit(match(types,rownames(feature_length))),]*1000/featurent
	ba <- plot(4:(length(freWG)+3),frefeature[1,],ylim = c(0,1.2*max(frefeature)),type = "l",pch = 20,col = 1, ylab = "#MH/kb_of_feature",xlab = "MH length(bp)",main = paste0(args[3],"_#MH/kb_of_feature"))
	for(i in 2:dim(frefeature)[1])
	{
		  lines(4:(length(freWG)+3),frefeature[i,],type = "o",pch = 20,cex=0.8,col = i,ylab = "#MH/kb_of_feature",xlab = "MH length(bp)")
	}
	legend("topright", legend=rownames(frefeature), col=c(1:dim(frefeature)[1]), lty=1)


	ba <- plot(4:(length(freWG)+3),log2(frefeature[1,]),type = "l",pch = 20,col = 1,ylim = c(-2*max(log2(frefeature)),2*max(log2(frefeature))), ylab = "log2(#MH/kb_of_feature)",xlab = "MH length(bp)",main = paste0(args[3],"_log2(#MH/kb_of_feature)"))
	for(i in 2:dim(frefeature)[1])
	{
		  lines(4:(length(freWG)+3),log2(frefeature[i,]),type = "o",pch = 20,cex=0.8,col = i,ylab = "#MH/kb_of_feature",xlab = "MH length(bp)")
	}
	legend("topright", legend=rownames(frefeature), col=c(1:dim(frefeature)[1]), lty=1)


	#figure of log2((#MH/kb_of_feature)/(#MH/kb_whole_genome))
	frefea_WG <- log2(t(t(frefeature)/freWG))
	ba <- plot(4:(length(freWG)+3),frefea_WG[1,],type = "l",col = 1,ylim= c(-2*max(frefea_WG),2*max(frefea_WG)),ylab = "log2(Fre-feature/Fre-WG)",xlab = "MH length(bp)",main = paste0(args[3],"feature-MHlength in features"))
	for(i in 2:dim(frefeature)[1])
	{
		  lines(4:(length(freWG)+3),frefea_WG[i,],type = "o",pch = 20,cex=0.8,col = i,ylab = "log2(Fre-feature/Fre-WG)",xlab = "MH length(bp)")
	}
	legend("bottomleft", legend=rownames(frefeature), col=c(1:dim(frefeature)[1]), lty=1)
	 
}


dev.off()

