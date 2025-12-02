#TF risk gene 相关性分析

load("D:/A-BS/GSE182109/result/182109_annotation_over.Rdata")
expression_matrix <- as(as.matrix(sce@assays$RNA@counts), 'sparseMatrix')

exprSet<- as.matrix(GetAssayData(object = sce, layer  = "counts"))
exp_new<-exprSet[which(rownames(exprSet) %in% gene),]
y <- exprSet[which(rownames(exprSet) %in% TF1),]

data<-rbind(exp_new,y)
data<-as.data.frame(t(data))
#############################
setwd("D:/A-BS/SCENIC")
load("4Cell_expMatrix.Rdata")
exp<-as.data.frame(t(exprMat))
save(exp,file = "4Cell_expMatrix2cor.Rdata")
#####################TCGA-CGGA exp
load("D:/A-BS/mime1/TCGA_CGGA_input_exp.Rdata")
exp<-as.data.frame(scale(t(exp_new[,4:14457])))
exp<-as.data.frame(t(exp))
#########################################
###################sc-level tumor exp
load("D:/A-BS/CellChat/182109_sce_chat_umap_over.Rdata")
scRNA<-subset(x=sce,idents=c("Glia and neuronal cell(Tumor-associated)","Oligodendrocyte(Tumor-associated)"))
exp  <-  as.matrix(scRNA@assays$RNA@data)
exp<-as.data.frame(t(exp))
#filter cell
info<-ct[which(ct$celltype %in% c("Glia and neuronal cell(Tumor-associated)","Oligodendrocyte(Tumor-associated)")),]
exp_new<-exp[which(rownames(exp) %in% info$cell),]

data<-exp_new[,which(colnames(exp_new) %in% TF1)]
data2 <- data[apply(data!= 0 , 1 , all),]


save(exp,file = "TumorCell_expMatrix2cor.Rdata")
#############sc-level T exp
load("D:/A-BS/CellChat/182109_sce_chat_umap_over.Rdata")
scRNA<-subset(x=sce,idents=c("CD8_T_EX","CD8_T_EM"))
exp  <-  as.matrix(scRNA@assays$RNA@data)
exp<-as.data.frame(t(exp))
##filter cell
library(xlsx)
ct<-read.xlsx("D:/A-BS/SCENIC/4cell_filtered_info.xlsx",sheetIndex = 1)
info<-ct[which(ct$celltype %in% c("CD8_T_EX","CD8_T_EM")),]
exp_new<-exp[which(rownames(exp) %in% info$cell),]

exp_new<-exp[,which(colnames(exp) %in% gene)]

#data<-exp_new[which(rowSums(exp_new) > 0),]

data1 <- data[apply(data!= 0 , 1 , all),]

save(exp_new,file = "T2Cell_expMatrix2cor.Rdata")

##############
TF1<-unique(c("SOX11","CEBPD","SOX4","EGR1"))
ICI<-c("ETS1","ELF1","RUNX3")
#gene<-c("AEBP1","ASF1A","DCC","HDAC5",
#        "IL13RA2","OPHN1","PRPS1")
gene<-ICI


#42669

colnames(exp)<-trimws(colnames(exp))
exp_new<-exp[,which(colnames(exp) %in% gene)]
exprSet<-exp
y <- exprSet[,which(colnames(exprSet) %in% TF1)]
data<-cbind(exp_new,y)

data<-exp_new[,which(colnames(exp_new) %in% TF1)]

#exprSet<-as.data.frame(t(exprMat))#转置

# 2. 计算某个基因和其它基因的相关性（以S100A8为例）-----#####
exprSet1<-exprSet[,which(colnames(exprSet) %in% gene)]

y1<-as.numeric(y$SOX4)
colnames<-colnames(exprSet1)
cor_data_df1<- data.frame(colnames)
for(i in 1:length(colnames)){
  test<-cor.test(as.numeric(exprSet1[,i]),y1,type="spearman") #可更换pearson
  cor_data_df1[i,2]<- test$estimate
  cor_data_df1[i,3]<- test$p.value
  
}
names(cor_data_df1)<-c("symbol","correlation","pvalue")
cor_data_df %>% head()

# 3. 筛选有意义的正相关和负相关的基因-----####
library(dplyr)
library(tidyr)
cor_data_sig_pos <- cor_data_df %>%
  dplyr::filter(pvalue <0.01)%>%dplyr::filter(correlation >0)%>%
  dplyr::arrange(desc(correlation))

cor_data_sig_neg <- cor_data_df %>%
  dplyr::filter(pvalue <0.01)%>%dplyr::filter(correlation <0)%>%
  dplyr::arrange(desc(abs(correlation)))

# 4. 随机选取正相关和负相关基因，分别作图验证----######
#1）S100A9正相关####
# install.packages("ggstatsplot") ##linux安装过程中报错，需先安装mpfr工具
# $conda install anaconda::mpfr
library(ggstatsplot)
pdf(paste0("GSE7696相关性ETS1-CD44.pdf"), width = 7, height = 6)
ggscatterstats(data = data2 ,
               y =SOX11,
               x =SOX4,
               centrality.para = "mean",               
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram")# #类型可以换成density,boxplot,violin,densigram                
               #title = "Relationship between S100A8 and S100A9")
dev.off()
