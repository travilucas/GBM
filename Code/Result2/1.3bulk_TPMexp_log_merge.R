#load tumor exp
load("D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA_GBMexp_Tumor_TPM.Rdata")
LG_TPMs<-TPMs
#
#change symbol
ENSG<-as.matrix(rownames(LG_TPMs))
#去除.后面的字符
a <- gsub("\\..*", "",ENSG)
rownames(LG_TPMs)<-a
library(clusterProfiler)
library(org.Hs.eg.db)

gene <- bitr(a, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
LG_TPMs = LG_TPMs[rownames(LG_TPMs) %in% gene$ENSEMBL,]

gene=gene[match(rownames(LG_TPMs),gene$ENSEMBL),]
rownames(LG_TPMs)<-gene$SYMBOL

tmp = by(LG_TPMs,
         gene$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))])

probes = as.character(tmp)
LG_TPMs =LG_TPMs[rownames(LG_TPMs) %in% probes,] # ?????ж???̽???Ļ???
rownames(LG_TPMs)=gene[match(rownames(LG_TPMs),gene$ENSEMBL),2]
save(LG_TPMs,file ="D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA_GBMexp_TumorGene_TPM.Rdata")
load("D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA_GBMexp_TumorGene_TPM.Rdata")
######
LG_TPM <- log2(LG_TPMs+1)
save(LG_TPM,file ="D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA_GBMexp_TumorGene_TPM_log.Rdata")

#####TCGA GBM log_TPM
#ENSG to symbol
#change symbol
tumor_exp<-read.table("D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA-GBM-Tumorexp_TPM.txt",sep="\t",header = T)
ENSG<-as.matrix(tumor_exp$V1)

gene <- bitr(ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
tumor_exp = tumor_exp[which(tumor_exp$V1 %in% gene$ENSEMBL),]

gene=gene[match(tumor_exp$V1,gene$ENSEMBL),]
gene<- gene[!duplicated(gene$SYMBOL), ]
tumor_exp<-tumor_exp[match(gene$ENSEMBL,tumor_exp$V1),]
identical(gene$ENSEMBL,tumor_exp$V1)

rownames(tumor_exp)<-gene$SYMBOL
exp<-as.data.frame(lapply(tumor_exp[,-1],as.numeric))
rownames(exp)<-rownames(tumor_exp)
save(exp,file ="D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA_GBMexp_TumorGene_TPM.Rdata" )
#######
#CGGA325
load("D:/A-BS/UCSC-GBM/CGGA/data/CGGA.mRNAseq_325Tumorexp_TPM.Rdata")
TPMS325<-log2(TPMs+1)
save(TPMS325,file ="D:/A-BS/UCSC-GBM/CGGA/data/CGGA.mRNAseq_325Tumorexp_TPM_log.Rdata" )
##CGGA693
load("D:/A-BS/UCSC-GBM/CGGA/data/CGGA.mRNAseq_693Tumorexp_TPM.Rdata")
TPMS693<-log2(TPMs+1)
save(TPMS693,file ="D:/A-BS/UCSC-GBM/CGGA/data/CGGA.mRNAseq_693Tumorexp_TPM_log.Rdata" )
####normal log
#load normal exp 110
load("D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105_TPM.Rdata")
gtex<-log2(TPMs+1)
save(gtex,file ="D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105_TPM_log.Rdata" )
#############
#log over start merge
#############
#tumor
load("D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA_GBMexp_TumorGene_TPM_log.Rdata")
load("D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA_GBMexp_TumorGene_TPM_log.Rdata")
load("D:/A-BS/UCSC-GBM/CGGA/data/CGGA.mRNAseq_325Tumorexp_TPM_log.Rdata")
load("D:/A-BS/UCSC-GBM/CGGA/data/CGGA.mRNAseq_693Tumorexp_TPM_log.Rdata")
TPMS325<-as.data.frame(TPMS325)
TPMS693<-as.data.frame(TPMS693)
#####################
LG_TPM<-LG_TPM[which(rownames(LG_TPM) %in% rownames(exp)),]
exp_new<-exp[which(rownames(exp) %in% rownames(LG_TPM)),]

LG_TPM <- LG_TPM[order(row.names(LG_TPM)), ]
exp_new <- exp_new[order(row.names(exp_new)), ]

identical(rownames(LG_TPM),rownames(exp_new))
TCGA_merge<-cbind(LG_TPM,exp_new)
###################
TPMS325<-TPMS325[which(rownames(TPMS325) %in% rownames(TPMS693)),]
TPMS693<-TPMS693[which(rownames(TPMS693) %in% rownames(TPMS325)),]

TPMS325 <- TPMS325[order(row.names(TPMS325)), ]
TPMS693 <- TPMS693[order(row.names(TPMS693)), ]
identical(rownames(TPMS325),rownames(TPMS693))
CGGA_merge<-cbind(TPMS325,TPMS693)
#################
TCGA_merge<-TCGA_merge[which(rownames(TCGA_merge) %in% rownames(CGGA_merge)),]
CGGA_merge<-CGGA_merge[which(rownames(CGGA_merge) %in% rownames(TCGA_merge)),]

TCGA_merge<- TCGA_merge[order(row.names(TCGA_merge)), ]
CGGA_merge <- CGGA_merge[order(row.names(CGGA_merge)), ]

identical(rownames(TCGA_merge),rownames(CGGA_merge))
Tumor_merge<-cbind(TCGA_merge,CGGA_merge)

write.table(exp1,file = "TCGA_Tumor_merge_over.txt",sep = "\t",quote = FALSE,row.names = FALSE)

save(Tumor_merge,file ="D:/A-BS/UCSC-GBM/GBM_Tumor_merge_over.Rdata" )
######################
#change gene symbol
#tcga-gbm-normal######
load("D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA-GBM-Normalexp_TPM_log.Rdata")
ENSG<-as.matrix(rownames(Normal_exp))
a <- gsub("\\..*", "",ENSG)
rownames(Normal_exp)<-a
library(clusterProfiler)
library(org.Hs.eg.db)
exp_new<-cbind(a,Normal_exp)
#
gene <- bitr(a, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")

ENSG<-as.matrix(tumor_exp$V1)

gene <- bitr(ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
exp_new = exp_new[which(exp_new$a %in% gene$ENSEMBL),]

gene=gene[match(exp_new$a,gene$ENSEMBL),]
gene<- gene[!duplicated(gene$SYMBOL), ]
exp_new<-exp_new[match(gene$ENSEMBL,exp_new$a),]
identical(gene$ENSEMBL,exp_new$a)

rownames(exp_new)<-gene$SYMBOL
exp<-as.data.frame(lapply(exp_new[,-1],as.numeric))
rownames(exp)<-rownames(exp_new)
save(exp,file ="D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA-GBM-NormalexpGene_TPM_log.Rdata" )
##########
load("D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105_TPM_log.Rdata")
#GTEX##########
#change symbol
ENSG<-as.matrix(rownames(gtex))
#去除.后面的字符
library(clusterProfiler)
library(org.Hs.eg.db)

gene <- bitr(ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
gtex =gtex[rownames(gtex) %in% gene$ENSEMBL,]

gene=gene[match(rownames(gtex),gene$ENSEMBL),]
rownames(gtex)<-gene$SYMBOL

tmp = by(gtex,
         gene$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))])

probes = as.character(tmp)
gtex =gtex[rownames(gtex) %in% probes,] # ?????ж???̽???Ļ???
rownames(gtex)=gene[match(rownames(gtex),gene$ENSEMBL),2]
#
save(gtex,file ="D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105_TPMGene_log.Rdata")
################################
#merge normal
#############
load("D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105_TPMGene_log.Rdata")
load("D:/A-BS/UCSC-GBM/TCGA/GBM/TCGA-GBM-NormalexpGene_TPM_log.Rdata")
exp<-exp[which(rownames(exp) %in% rownames(gtex)),]
gtex<-gtex[which(rownames(gtex) %in% rownames(exp)),]

exp <- exp[order(row.names(exp)), ]
gtex <- gtex[order(row.names(gtex)), ]
identical(rownames(exp),rownames(gtex))
Normal_merge<-cbind(exp,gtex)
save(Normal_merge,file ="D:/A-BS/UCSC-GBM/GBM_Normal_merge_over.Rdata" )
