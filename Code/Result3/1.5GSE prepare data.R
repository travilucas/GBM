library(xlsx)
#GSE108474
exp <- read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/GSE108474_series_matrix.txt",comment.char = "!",sep = "\t",header=T,row.names = 1,fill=TRUE)
info<-read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/GSE108474_REMBRANDT_clinical.data.txt",sep = "\t",header=T,fill=TRUE)
gsm<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/GSE108474_info.xlsx",sheetIndex = 1)
gsm<-gsm[,-1]
identical(colnames(gsm),colnames(exp))
colnames(exp)<-gsm[1,]
info_new<-info[match(gsm[1,],info$SUBJECT_ID),]
write.xlsx(info_new,"D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/info_sort.xlsx")
###
gpl<-read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/GPL570-55999.txt",header = T,sep = "\t",fill = T)
gpl<-gpl[,c(1,11)]
write.xlsx(gpl,"D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/gpl_sort.xlsx")

exp<-exp[which(rownames(exp) %in% gpl$ID),]
gpl<-gpl[which(gpl$ID %in% rownames(exp) ),]
gpl$Gene.Symbol<-sub("\\///.*", "", gpl$Gene.Symbol)

#
info_new<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/info_sort.xlsx",sheetIndex = 2)
exp_new<-exp[,which(colnames(exp) %in% info_new$SUBJECT_ID)]

######################################
#GSE42669
library(xlsx)
exp <- read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE42669/GSE42669_series_matrix.txt",comment.char = "!",sep = "\t",header=T,row.names = 1,fill=TRUE)
gpl<-read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE42669/GPL6244-17930.txt",sep = "\t",header=T,fill=TRUE)
info<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE42669/GSE42669_info.xlsx",sheetIndex = 1)

gpl<-gpl[,c(1,10)]
gpl$symbol <- sapply(gpl$gene_assignment, function(x) strsplit(x, "//")[[1]][2])
gpl<-na.omit(gpl)
exp<-exp[which(rownames(exp) %in% gpl$ID),]
#
tmp = by(exp,
         gpl$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])

probes = as.character(tmp)
dim(exp)
exp = exp[rownames(exp) %in% probes,] # ?????ж???̽???Ļ???
dim(exp)

rownames(exp)=(gpl[match(rownames(exp),gpl$ID),3])

info<-as.data.frame(t(info))
info_new<-info[-1,12:13]

info_new$OS.time <- sapply(info_new$V12, function(x) strsplit(x, ":")[[1]][2])
info_new$OS <- sapply(info_new$V13, function(x) strsplit(x, ":")[[1]][2])

info_new$OS.time<-as.numeric(trimws(info_new$OS.time))
info_new$OS<-as.numeric(trimws(info_new$OS))
info_new<-na.omit(info_new)
info_new$SUBJECT_ID<-rownames(info_new)
exp_new<-exp[,which(colnames(exp) %in% info_new$SUBJECT_ID)]
info_new<-info_new[match(colnames(exp_new),info_new$SUBJECT_ID),]
identical(colnames(exp_new),info_new$SUBJECT_ID)
exp<-cbind(info_new$OS.time,info_new$OS,t(exp_new))
colnames(exp)[1:2]<-c("OS.time","OS")
exp<-as.data.frame(exp)
exp$OS.time<-exp$OS.time *7
save(exp,file = "D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE42669/GSE42669exp_over.Rdata")
###########

#######################
#GSE16011
library(xlsx)
exp <- read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE16011/GSE16011_series_matrix.txt",comment.char = "!",sep = "\t",header=T,row.names = 1,fill=TRUE)
gpl<-read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE16011/GPL8542.txt",sep = "\t",header=T,fill=TRUE)
info<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE16011/GSE16011_info.xlsx",sheetIndex = 1)
gsm<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE16011/GSE16011_gsm.xlsx",sheetIndex = 1)
gsm<-gsm[,-1]
identical(colnames(gsm),colnames(exp))
gsm[1,]<-sapply(gsm[1,], function(x) strsplit(x, " ")[[1]][2])
colnames(exp)<-gsm[1,]
info_new<-info[match(gsm[1,],info$Database.number),c(1,9,10)]
info_new<-na.omit(info_new)
gpl<-na.omit(gpl[,1:2])

library(clusterProfiler)
library("org.Hs.eg.db")
gene = bitr(gpl$ORF, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
id = gene$ENTREZID
gpl<-gpl[match(id,gpl$ORF),]
gpl<-cbind(gpl,gene)

exp<-exp[which(rownames(exp) %in% gpl$ID),]

tmp = by(exp,
         gpl$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))])

rownames(exp)=(gpl[match(rownames(exp),gpl$ID),4])
#
info_new$Survival..years.<-gsub(",", ".",info_new$Survival..years.)
info_new$Alive<-gsub("Dead", "1",info_new$Alive)
info_new$Alive<-gsub("Alive", "0",info_new$Alive)
#
info_new$Alive<-as.numeric(info_new$Alive)
info_new$Survival..years.<-as.numeric(info_new$Survival..years.)
info_new<-na.omit(info_new)
exp_new<-exp[,which(colnames(exp) %in% info_new$Database.number)]
exp_new<-as.data.frame(t(exp_new))
exp<-cbind(info_new[,c(1,3,2)],exp_new)
rownames(exp)<-rownames(exp_new)
colnames(exp)[1:3]<-c("ID","OS.time","OS")
exp$OS.time<-exp$OS.time *365
save(exp,file ="D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE16011/GSE16011exp_over.Rdata" )
###########
#GSE7696
exp <- read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE7696/GSE7696_series_matrix.txt",comment.char = "!",sep = "\t",header=T,row.names = 1,fill=TRUE)
gpl<-read.table("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE7696/GPL570-55999.txt",sep = "\t",header=T,fill=TRUE)
library(xlsx)
info<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE7696/GSE7696_info.xlsx",sheetIndex = 1)

gpl<-gpl[,c(1,11)]
gpl$symbol <- sapply(gpl$Gene.Symbol, function(x) strsplit(x, "///")[[1]][1])

gpl<-na.omit(gpl)
exp<-exp[which(rownames(exp) %in% gpl$ID),]
gpl<-gpl[which(gpl$ID %in% rownames(exp) ),]

#
tmp = by(exp,
         gpl$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])

probes = as.character(tmp)
dim(exp)
exp = exp[rownames(exp) %in% probes,] # ?????ж???̽???Ļ???
dim(exp)

rownames(exp)=(gpl[match(rownames(exp),gpl$ID),3])

info<-as.data.frame(t(info[,-1]))
info_new<-info[,1:2]

info_new$OS.time <- sapply(info_new$V2, function(x) strsplit(x, ":")[[1]][2])
info_new$OS <- sapply(info_new$V1, function(x) strsplit(x, ":")[[1]][2])
info_new$OS.time<-as.numeric(trimws(info_new$OS.time))
info_new$OS<-as.numeric(trimws(info_new$OS))
info_new<-na.omit(info_new)

exp_new = exp[,colnames(exp) %in% rownames(info_new)] # ?????ж???̽???Ļ???
exp_new<-as.data.frame(t(exp_new))
exp_new<-exp_new[order(row.names(exp_new)),] 
info_new<-info_new[order(row.names(info_new)),] 
identical(rownames(exp_new),rownames(info_new))
exp<-cbind(info_new[,3:4],exp_new)
exp$OS.time <- exp$OS.time *31
save(exp,file = "D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE7696/GSE7696exp_over.Rdata")
