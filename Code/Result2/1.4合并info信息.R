#merge tumor info
#order: LGGGBM GBM CGGA325 693
library(xlsx)
LG_info<-read.xlsx("D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA_LGGGBM_ALLinfo2.xls",sheetIndex = 1)
LG_info<-LG_info[,-4]
colnames(LG_info)<-c("Sample","Age","Gender","OS","OS.time")
LG_info$Platform<-c(rep("TCGA-LGGGBM",399))
##
GBM_info<-read.xlsx("D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA_LGGGBM_ALLinfo_over.xlsx",sheetIndex = 1)
GBM_info<-GBM_info[,-4]
colnames(GBM_info)<-c("Sample","Age","Gender","OS","OS.time")
GBM_info$Platform<-c(rep("TCGA-GBM",143))

##########
setwd("D:/A-BS/UCSC-GBM/TCGA/GBM/")
sur<-read.table("TCGA-GBM.survival.tsv",header = T,sep = "\t")
info<-read.csv("TCGA-GBM.GDC_phenotype.tsv",header = T,sep = "\t",fill = T)
group<-read.csv("TCGA-GBM_sample_id.csv",header = T)
##
group$sample_id<-gsub('[-]', '.', group$sample_id)
info$submitter_id.samples<-gsub('[-]', '.', info$submitter_id.samples)
sur$sample<-gsub('[-]', '.', sur$sample)
sur_new<-sur[match(group$sample_id,sur$sample),]
sur_new<-sur_new[,c(1,2,4)]
sur_new<-na.omit(sur_new)
info_new<-info[which(info$submitter_id.samples %in% group$sample_id),c(1,2,51,53)]
info_new<-na.omit(info_new)
info_new<-info_new[match(sur_new$sample,info_new$submitter_id.samples),]
identical(sur_new$sample,info_new$submitter_id.samples)
info_over<-cbind(info_new,sur_new)
write.xlsx(info_over,file = "D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA_LGGGBM_ALLinfo_over.xlsx")

##########
CGGA325_info<-read.xlsx("D:/A-BS/UCSC-GBM/CGGA/CGGA325_GBM_info.xlsx",sheetIndex = 1)  
CGGA693_info<-read.xlsx("D:/A-BS/UCSC-GBM/CGGA/CGGA693_GBM_info.xlsx",sheetIndex = 1)
CGGA325_info<-CGGA325_info[,c(1,5,6,7,8)]
CGGA693_info<-CGGA693_info[,c(1,5,6,7,8)]
#alive=0 dead=1
colnames(CGGA325_info)[c(1,4:5)]<-c("Sample","OS.time","OS")
colnames(CGGA693_info)[c(1,4:5)]<-c("Sample","OS.time","OS")
CGGA325_info$Platform<-c(rep("CGGA325",71))
CGGA693_info$Platform<-c(rep("CGGA693",105))

#
LG_info <- LG_info[,order(colnames(LG_info)) ]
GBM_info<- GBM_info[,order(colnames(GBM_info))]
CGGA325_info<-CGGA325_info[,order(colnames(CGGA325_info))]
CGGA693_info<-CGGA693_info[,order(colnames(CGGA693_info))]
#
Tumor_info<-rbind(LG_info,GBM_info,CGGA325_info,CGGA693_info)
Tumor_info$age1 <- ifelse((Tumor_info$Age<60)==TRUE,"<60",">=60")

library(stringr)
Tumor_info$Gender<-str_to_title(Tumor_info$Gender) 
write.xlsx(Tumor_info,file = "D:/A-BS/UCSC-GBM/GBM-Tumor_info.xlsx")
#############



