library(devtools)
devtools::install_github('jmzeng1314/AnnoProbe')

devtools::install_github('xjsun1221/tinyarray',upgrade = F)
library(data.table)
library(tidyverse)
library(AnnoProbe)
library(tinyarray)

setwd("D:/A-BS/UCSC-GBM/")
dat<-data.table::fread("gtex_gene_expected_count.gz",data.table = F)

exp<-column_to_rownames(dat,"sample") %>% as.matrix()

rownames(exp)<-rownames(exp) %>% str_split("\\.",simplify = T) %>% .[,1]

an<-annoGene(rownames(exp),ID_type = "ENSEMBL")

exp <- trans_array(exp,ids = an,from="ENSEMBL",to="SYMBOL")
#
clinical<-data.table::fread("GTEX_phenotype.gz")
table(clinical$`_primary_site`)

clinical<-clinical[clinical$`_primary_site` != "<not provided>"]

colnames(clinical)[3]="site"
clinical.subset<-subset(clinical,site=="Brain")
s<-intersect(colnames(exp),clinical.subset$Sample)
clinical.subset<-clinical.subset[match(s,clinical.subset$Sample),]
exp<-exp[,s]

identical(clinical.subset$Sample,colnames(exp))
dim(exp)

#count是通过log2(expected_count+1)得到，转化回来
exp<-round(2^exp-1,4)
#写出的是所有疾病的正常样本
write.table(data.frame(ID=rownames(exp),exp),file = "GTEX_Brain_Normal_sample.txt",sep = "\t",quote = F,row.names = F)
write.table(clinical.subset,file = "GTEX_Brain_Normal_sample_phenotype.txt",sep = "\t",quote = F,row.names = F)
#############
#GTEX normal
gtex_exp<-read.table("D:/A-BS/UCSC-GBM/GTEx/GTEX_Brain_Normal_sample.txt",sep = '\t',row.names = 1,header = T)
gtex_info<-read.table("D:/A-BS/UCSC-GBM/GTEx/GTEX_Brain_Normal_sample_phenotype.txt",sep = '\t',row.names = 1,header = T)

gtex_info1<-gtex_info[which(gtex_info$body_site_detail..SMTSD.=="Brain - Cortex"),]
gtex_sample<-as.matrix(colnames(gtex_exp))

rownames(gtex_info1)<-gsub('[-]', '.', rownames(gtex_info1))
gtex_sample<-gtex_sample[which(gtex_sample %in% rownames(gtex_info1))]

gtex_exp<-gtex_exp[,which(colnames(gtex_exp) %in% gtex_sample)]
save(gtex_exp,file="D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105.Rdata")

load("D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105.Rdata")
