setwd("D:/A-BS/SCENIC/SCENIC_result")
library(Seurat) 
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

load("D:/A-BS/SCENIC/4Cell_expMatrix.Rdata")
load("D:/A-BS/SCENIC/4Cell_InfoMatrix.Rdata")
load("D:/A-BS/CellChat/182109_sce_chat_umap_over.Rdata")
scRNAsub<-subset(x=sce,idents=c("Glia and neuronal cell(Tumor-associated)","Oligodendrocyte(Tumor-associated)","CD8_T_EX","CD8_T_EM"))



AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
library(Seurat)
scRNAauc <- AddMetaData(scRNAsub, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'D:/A-BS/SCENIC/scRNAauc.rds')

##??????????regulonAUC????
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNAsub, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'D:/A-BS/SCENIC/scRNAbin.rds')

library(pheatmap)

celltype = subset(cellInfo,select = 'CellType')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#??ѡ???ָ???Ȥ??regulons
setwd("D:/A-BS/SCENIC/SCENIC_result")
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

load("D:/A-BS/SCENIC/4Cell_InfoMatrix.Rdata")
sub_regulonAUC <- regulonAUC[,match(rownames(cellInfo),colnames(regulonAUC))]

### 4.1. TF活性均值
# 看看不同单细胞亚群的转录因子活性平均值
# Split the cells by cluster:
selectedResolution <- "CellType" # select resolution
cellsPerGroup <- split(rownames(cellInfo), 
                       cellInfo[,selectedResolution])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 

rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellInfo[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss,zThreshold = 0.1,thr=0.1)
rssPlot

a<-rssPlot[["rowOrder"]]

a <- gsub(' \\(','_',a)
a<- gsub('\\)','',a)
my.regulons=a

myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]

#ʹ??regulonԭʼAUCֵ??????ͼ
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype)
         #filename = 'scenic_seurat/myAUCmatrix_heatmap.png',
         #width = 6, height = 5)
ct<-cbind(rownames(celltype),celltype$CellType)
ct<-ct[match(colnames(BINmatrix),ct[,1]),]
rownames(ct)<-ct[,1]
ct<-ct[,-1]
ct<-as.data.frame(ct)
ann_colors<-list(
  ct=c("CD8_T_EM"="#3778ad","CD8_T_EX"="#4ea74a",
       "Glia and neuronal cell(Tumor-associated)"="#e89118",
       "Oligodendrocyte(Tumor-associated)"="#56a8d7")
)
pheatmap(myBINmatrix,
         show_colnames = F,
         cluster_rows =T,#是否对行聚类
         annotation_col=ct,annotation_colors =ann_colors,
         color = colorRampPalette(colors = c("white","#9F0000"))(100),
         )
                  #filename = 'scenic_seurat/myBINmatrix_heatmap.png',
                  #color = colorRampPalette(colors = c("white","black"))(100),
                  #width = 6, height = 5)