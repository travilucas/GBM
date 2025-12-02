#prepare data
load("D:/A-BS/UCSC-GBM/GBM_Tumor_merge_Combatover.Rdata")
tumor_combat<-as.data.frame(tumor_combat)
library(xlsx)
info<-read.xlsx("D:/A-BS/UCSC-GBM/GBM-Tumor_info.xlsx",sheetIndex = 1)

info$Sample<-gsub('[-]', '.', info$Sample)
colnames(tumor_combat)<-gsub('[-]', '.',colnames(tumor_combat))

info_new<-info[info$Sample %in% colnames(tumor_combat),]
tumor_combat<-tumor_combat[,colnames(tumor_combat) %in% info_new$Sample]

exp<-t(tumor_combat)
exp<-as.data.frame(exp)
info_new<-info_new[match(rownames(exp),info_new$Sample),]
identical(rownames(exp),info_new$Sample)
exp_new<-cbind(rownames(exp),info_new$OS.time,info_new$OS,exp)
colnames(exp_new)[1:3]<-c("ID","OS.time","OS")
save(exp_new,file = "D:/A-BS/mime1/TCGA_CGGA_input_exp.Rdata")
##################
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


#
library(xlsx)
info_new<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/info_sort.xlsx",sheetIndex = 2)
exp_new<-exp[,which(colnames(exp) %in% info_new$SUBJECT_ID)]

gpl<-read.xlsx("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/gpl_sort.xlsx",sheetIndex = 1)

gpl$Gene.Symbol<-sub("\\///.*", "", gpl$Gene.Symbol)
gpl<-gpl[match(rownames(exp_new),gpl$ID),]
gpl<-na.omit(gpl)
exp_new<-exp_new[which(rownames(exp_new) %in% gpl$ID),]
identical(rownames(exp_new),gpl$ID)

tmp = by(exp_new,
         gpl$Gene.Symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])

probes = as.character(tmp)
dim(exp)
exp_new = exp_new[rownames(exp_new) %in% probes,] # ?????ж???̽???Ļ???
dim(exp)

rownames(exp_new)=(gpl[match(rownames(exp_new),gpl$ID),3])

exp_new<-t(exp_new)
identical(rownames(exp_new),info_new$SUBJECT_ID)
exp<-cbind(info_new[,c(8,7)],exp_new)
exp$OVERALL_SURVIVAL_MONTHS<-exp$OVERALL_SURVIVAL_MONTHS*31
colnames(exp)[1:2]<-c("OS.time","OS")
save(exp,file = "D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/GSE108474exp_over.Rdata")
