library(dplyr)
library(Seurat)
library(patchwork)
library(irlba)
#setwd("D:/A-BS/GSE182109/result/")
###########
##读取数据
########
setwd("D:/A-BS/GSE182109_RAW_ALL/")

scRNAlist<-list()
for(i in 1:length(list.files("./"))){
  counts<-Read10X(data.dir =list.files("./")[i],gene.column = 2)
  scRNAlist[[i]]<-CreateSeuratObject(counts,project = list.files("./")[i],min.cells=3,min.features=200)
  print(i)
}
############
###分别质控
###########
for(i in 1:length(scRNAlist)){
  scRNA<-scRNAlist[[i]]
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  scRNAlist[[i]]<-scRNA
}
######整合为一个Seurat对象
scRNA<-merge(scRNAlist[[1]],y=c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]],scRNAlist[[7]],scRNAlist[[8]],scRNAlist[[9]],scRNAlist[[10]],scRNAlist[[11]],scRNAlist[[12]],scRNAlist[[13]],scRNAlist[[14]],scRNAlist[[15]]))
ncol(scRNA)#细胞数
nrow(scRNA)#基因数
################
#导入其他样本信息
#############
library(readxl)
info<-read_xlsx("D:/A-BS/GSE182109/data/info.xlsx",sheet=2)
scRNA@meta.data$sample = info[match(scRNA@meta.data$orig.ident,info$GSM),2]
scRNA@meta.data[["sample"]]<-scRNA@meta.data[["sample"]][["sample_id"]]
#
save(scRNA,file="D:/A-BS/GSE182109/result/182109_qc_over.Rdata")
#########################qc over
load("D:/A-BS/GSE182109/result/182109_qc_over.Rdata")
############
#标准化
############
scRNA <- NormalizeData(scRNA)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA<-ScaleData(scRNA,features = VariableFeatures(object = scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
ElbowPlot(scRNA)
#准确确定线性降维的维数
sce=scRNA
#
pct<-sce[["pca"]]@stdev/sum(sce[["pca"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu>90 & pct<5)[1]
co1
co2<-sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1),decreasing = T)[1]+1
co2
pcs<-min(co1,co2)
pcs
plot_df<-data.frame(pct=pct,cumu=cumu,rank=1:length(pct))

library(ggplot2)
ggplot(plot_df,aes(cumu,pct,label=rank,color=rank>pcs))+
  geom_text()+
  geom_vline(xintercept = 90,color='grey')+
  geom_hline(yintercept = min(pct[pct>5]),color="grey")+
  theme_bw()
###########确定维数 修改dims
scRNA <- FindNeighbors(scRNA, dims = 1:16)
scRNA <- FindClusters(scRNA, resolution = c(1:10/10))
#细胞聚类#############
library(clustree)
p1<-clustree(scRNA@meta.data,prefix='RNA_snn_res.')+coord_flip()
p1
###
#UMAP
scRNA<- RunUMAP(scRNA, dims = 1:16)
#####
#去批次处理
######
library(harmony)
Idents(scRNA)<-'sample'
sce=scRNA %>% RunHarmony("orig.ident",plot_convergence=TRUE)
#########
#再次umap
########
sce <- sce %>% 
  RunUMAP(reduction="harmony",dims=1:16) %>%
  FindNeighbors(reduction="harmony",dims=1:16) %>%
  FindClusters(resolution=c(1:10/10)) %>%
  identity()
#
save(sce,file="D:/A-BS/GSE182109/result/182109_harmony_over.Rdata")
########################################
#start annotation
###########
load("D:/A-BS/GSE182109/result/182109_harmony_over.Rdata")

Idents(sce)<-"RNA_snn_res.0.2"
Idents(sce)<-"celltype"

library(readxl)
ct<-read_xlsx("D:/A-BS/GSE182109/data/cluster_annotation.xlsx",sheet=1)
clusters <- sce@meta.data$RNA_snn_res.0.2

sce@meta.data$celltype = ct[match(clusters,ct$ClusterID),'celltype']

sce@meta.data$celltype<-sce@meta.data[["celltype"]][["celltype"]]

scRNA.markers <- FindAllMarkers(scRNA,logfc.threshold = 0.5,only.pos = T)

scRNA.markers_used<- scRNA.markers %>%
  group_by(cluster)%>%
  slice_max(n=30,order_by = avg_log2FC)%>%
  dplyr::filter(pct.1/pct.2>1.5)

marker_list<-data.frame(cluster=unique(scRNA.markers_used$cluster),gene="")
for(i in 1:length(unique(scRNA.markers_used$cluster))){
  marker_list$gene[i]<-paste(scRNA.markers_used$gene[which(scRNA.markers_used$cluster==i-1)],collapse = ',')
}

#############准备inferCNV的输入数据
sce1 = sce[, Idents(sce) %in% c("B cell","Endothelial cell","Glia and neuronal cell","Oligodendrocyte","Pericyte","T cell")]
save(sce1,file="D:/A-BS/GSE182109/result/182109_annotation_no_myeloid.Rdata")
#
load("D:/A-BS/GSE182109/result/182109_annotation_no_myeloid.Rdata")
write.table(meta, "D:/A-BS/GSE182109/result/182109_metadata_no_myeloid.txt",sep = '\t',quote = F,col.names = F)
