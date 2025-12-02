#### bulk去除批次效应
# 安装 install.packages("devtools")
#devtools::install_github("zhangyuqing/sva-devel")
#################################
#将所有表达谱的数据转换成tpm
#############################
#Gtex normal count
load("D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105.Rdata")
############get Count_gene efflen
#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
library(txdbmaker)
txdb <- makeTxDbFromGFF("gencode.v36.annotation.gtf",format="gtf")
# 获取每个基因id的外显子数据
exons.list.per.gene <- exonsBy(txdb,by="gene")
# 对于每个基因，将所有外显子减少成一组非重叠外显子，计算它们的长度(宽度)并求和
exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))
# 得到geneid和长度数据
gfe <- data.frame(gene_id=names(exonic.gene.sizes),
                  length=exonic.gene.sizes)
head(gfe)[1:5,1:2]
#                               gene_id length
# ENSG00000000003.15 ENSG00000000003.15   4536
# ENSG00000000005.6   ENSG00000000005.6   1476
# ENSG00000000419.13 ENSG00000000419.13   1207
# ENSG00000000457.14 ENSG00000000457.14   6883
# ENSG00000000460.17 ENSG00000000460.17   5970
save(gfe,file = "gfe_length.Rdata")
#################
#count to tpm
load("D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA-GBM-Origion_count_Tumor2.Rdata")
load("D:/A-BS/UCSC-GBM/TCGA/LGGGBM/gfe_length.Rdata")

gene<-rownames(exp_new1)
library(clusterProfiler)
library(org.Hs.eg.db)
gene_ENSG <- bitr(gene, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")

exp_new1<-exp_new1[which(rownames(exp_new1) %in% gene_ENSG$SYMBOL),]
gene_ENSG<-gene_ENSG[match(rownames(exp_new1),gene_ENSG$SYMBOL),]
identical(gene_ENSG$SYMBOL,rownames(exp_new1))
rownames(exp_new1)<-gene_ENSG$ENSEMBL
# 去除重复的映射，只保留每个ENSG编号的第一个符号
genes_unique <- gene_ENSG[!duplicated(gene_ENSG$ENSEMBL), ]

exp<-exp_new1[which(rownames(exp_new1) %in% genes_unique$SYMBOL),]
genes_unique<-genes_unique[match(rownames(exp),genes_unique$SYMBOL),]
identical(genes_unique$SYMBOL,rownames(exp))
rownames(exp)<-genes_unique$ENSEMBL

gfe$ENSG <- sub("\\..*", "", gfe$gene_id)
gfe<-gfe[which(gfe$ENSG %in% rownames(exp)),]
exp<-exp[which(rownames(exp) %in% gfe$ENSG),]
gfe_new<-gfe[match(rownames(exp),gfe$ENSG),]
identical(rownames(exp),gfe_new$ENSG)
#[1] TRUE 行名是能够对上的
gfe_new$gene_id<-gfe_new$ENSG
rownames(gfe_new)<-gfe_new$ENSG
gfe_new<-gfe_new[,-3]

effLen = gfe_new$length
# 查看转换的结果
Counts2TPM <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
TPMs <- apply(exp,2,Counts2TPM, effLen = effLen)
TPMs<-as.data.frame(TPMs)
save(TPMs,file = "D:/A-BS/UCSC-GBM/TCGA/LGGGBM/TCGA_GBMexp_Tumor_TPM.Rdata")
save(TPMs,file = "D:/A-BS/UCSC-GBM/GTEx/GTEx_Normal105_TPM.Rdata")
###############################
setwd("D:/A-BS/UCSC-GBM/CGGA/data/")
exp<-read.table("CGGA.mRNAseq_693.RSEM-genes.20200506.txt",header = T,sep = "\t",row.names = 1)
library(xlsx)
info<-read.xlsx("D:/A-BS/UCSC-GBM/CGGA/CGGA693_GBM_info.xlsx",sheetIndex = 1)
exp_new<-exp[,which(colnames(exp) %in% info$CGGA_ID)]
#FKPM to tpm
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
TPMs <- apply(exp_new,2,FPKM2TPM)
save(TPMs,file = "CGGA.mRNAseq_693Tumorexp_TPM.Rdata")
