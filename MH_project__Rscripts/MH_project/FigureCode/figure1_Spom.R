#!/apps/bioinfo/R-3.6.0share/bin/Rscipt
#source  /appsnew/source/R-3.6.0share.sh
args=commandArgs(T)
#inter_feature <- read.table(args[1],header = F) #file1 is pombe_n.intersect.txt
inter_feature <- read.table("pombe_n.intersect.txt",header = F,fill = T) #file1 is pombe_n.intersect.txt
num <- table(inter_feature[,1])
types <- names(num)
num <- as.data.frame(num)
num#MH of each feature

#location <- read.table(args[2],header = F, comment.char = "#") #pombe_n.gff3.bed
location <- read.table("pombe_n.gff3.bed",header = F, comment.char = "#") #pombe_n.gff3.bed
location$nt <- location$V3-location$V2+1 #V3 is start, V4 is stop 
total_nt <- as.data.frame(aggregate(location$nt,by=list(location$V4), sum)) #total nt of each feature
total_nt#total nt of each feature

overlap <- intersect(as.character(num$Var1),as.character(total_nt$Group.1))
ind1 <- match(overlap,as.character(total_nt$Group.1))
nt<- total_nt[[2]][ind1]
ind2 <- match(overlap,as.character(num$Var1))
MH <- num[[2]][ind2]
MH_kb <- MH*1000/nt
type_MHkb <- data.frame(overlap,MH,nt,MH_kb)
type_MHkb#type_MHkb: table of each feature and nt#


######plot "intron","CDS","five_prime_UTR","three_prime_UTR","ncRNA","nuclear_mt_pseudogene"
types <- c("intron","CDS","five_prime_UTR","three_prime_UTR","ncRNA","nuclear_mt_pseudogene")
types_tab <- na.omit(type_MHkb[match(types,type_MHkb[,1]),])
par(las=2,mar=c(12, 4.1, 4.1, 2.1))
if(!is.null(types_tab)){
  ##barplot: overlap of MH+1bp in each feature
  count<-types_tab[,4]
  names(count) <- types_tab[,1]
  ba <- barplot(count,col = "white",ylab = "MH pairs / Kilobases",space = 0.5,main =paste0("Spom_MH/kb"))
  
  #barplot: log2(Fre-feature/Fre-WG) overlap of MH+1bp in each feature
  nt_wg <- total_nt[[2]][match("chromosome",as.character(total_nt$Group.1))]    #the length of whole genome
  fre_fre <- log2(count/(1000*sum(types_tab[,2])/nt_wg))
  names(count) <- types_tab[,1]
  ba <- barplot(fre_fre,col = "white",ylab = "log2(Fre-feature/Fre-WG)",space = 0.5,main = paste0("Spom_log2(Fre-feature/Fre-WG)")) #,ylim = c(-0.8,0.8)
}




dev.off()

#pdf(paste0(args[3],'.figure1.pdf')) #profiles_Spom/figure1.pdf

##barplot: overlap of MH+1bp in each feature
name_dele <- c("chromosome","snoRNA", "snRNA", "tRNA", "rRNA","five_prime_UTR","pseudogenic_transcript","gene")
count <- type_MHkb[-match(name_dele,type_MHkb$types),4]
names(count) <- type_MHkb[-match(name_dele,type_MHkb$types),1]
ba <- barplot(count,col = "white",ylab = "MH pairs / Kilobases",space = 0.5,main = paste0("Scer_MH/kb"))

#barplot: log2(Fre-feature/Fre-WG) overlap of MH+1bp in each feature
nt_wg <- total_nt[[2]][match("chromosome",as.character(total_nt$Group.1))]    #the length of whole genome
fre_fre <- log2(count/(1000*sum(type_MHkb[-match(name_dele,type_MHkb$types),2])/nt_wg))
names(count) <- type_MHkb[-match(name_dele,type_MHkb$types),1]
ba <- barplot(fre_fre,col = "white",ylab = "log2(Fre-feature/Fre-WG)",ylim = c(-0.5,0.5),space = 0.5,main = paste0("Scer_log2(Fre-feature/Fre-WG)")) #,ylim = c(-0.8,0.8)

