library(Seurat)
library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(glmGamPoi)
library(Rfast2)
##################
# Load the expression data
#该路径下有h5文件以及spatial文件夹(5个文件)
sce = Load10X_Spatial(data.dir ="D:/A-BS/ST_GBM/GSE237183/GSE237183_RAW/zh8812/new"
                      ,slice = "zh8812")
##########
st <- SCTransform(sce, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
##
#st <- FindSpatiallyVariableFeatures(st, assay = "SCT", selection.method = "moransi")
st <- RunPCA(st, assay = "SCT", verbose = FALSE)

pct<-st[["pca"]]@stdev/sum(st[["pca"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu>90 & pct<5)[1]
co1
co2<-sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1),decreasing = T)[1]+1
co2
pcs<-min(co1,co2)
pcs

st <- FindNeighbors(st, reduction = "pca", dims = 1:pcs)
st <- FindClusters(st, verbose = FALSE)
st <- RunUMAP(st, reduction = "pca", dims = 1:pcs)

saveRDS(st,file = "D:/A-BS/ST_GBM_new/data/SpatialData237183_zh1019inf.rds")
