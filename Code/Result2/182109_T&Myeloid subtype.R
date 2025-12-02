load("D:/A-BS/GSE182109/result/182109_annotation_over.Rdata")
library(dplyr)
library(Seurat)
library(patchwork)
#subset T cells
scRNA<-subset(x=sce,idents="T cell")

save(scRNA,file = "D:/A-BS/GSE182109/result/182109_Tcell_sce.Rdata")
#subset myeloid cells
scRNA<-subset(x=sce,idents="Myeloid cell")

save(scRNA,file = "D:/A-BS/GSE182109/result/182109_Myeloidcell_sce.Rdata")

##########
#数据降维
#RunPCA
load("D:/A-BS/GSE182109/result/182109_Myeloidcell_sce.Rdata")
scRNA1<-ScaleData(scRNA)
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(object = scRNA1))
#绘制肘部图可视化确定聚类个数
ElbowPlot(scRNA1)
#
#准确确定线性降维的维数
sce=scRNA1
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
#remember change 参数
scRNA1 <- FindNeighbors(scRNA1, dims = 1:14)
scRNA1 <- FindClusters(scRNA1, resolution = c(1:10/10))
#细胞聚类#############
library(clustree)
p1<-clustree(scRNA1@meta.data,prefix='RNA_snn_res.')+coord_flip()
p1
colnames(scRNA1@meta.data)
# Look at cluster IDs of the first 5 cells
head(Idents(scRNA1), 5)
#UMAP
scRNA1<- RunUMAP(scRNA1, dims = 1:14)

save(scRNA1,file = "D:/A-BS/GSE182109/result/182109_Myeloidcell_sce_not_annotation.Rdata")
Idents(scRNA1)<-scRNA1$RNA_snn_res.0.8
#################
library(SingleR)
library(celldex)
library(dplyr)
library(Seurat)
library(patchwork)

load("D:/A-BS/marker/ref_Human_all.RData")
#celldex::HumanPrimaryCellAtlasData()
hpca.se<-ref_Human_all
hpca.se@metadata
hpca.se@colData@metadata

#bpe.se<-celldex::BlueprintEncodeData() 
#bpe.se@metadata
#Idents(scRNA_combined)<-"integrated_snn_res.0.3"
sce1<-scRNA1

sce_for_SingleR <- GetAssayData(sce1, layer="data")
sce_for_SingleR
clusters <- scRNA1@meta.data$RNA_snn_res.0.8
pred.hesc <- SingleR(test = sce_for_SingleR, 
                     ref = hpca.se, 
                     labels = hpca.se$label.main,
                     method = "cluster", 
                     clusters = clusters, 
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")

table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 

#celltype$celltype <-read.delim("clipboard",header=FALSE)

sce1@meta.data$celltype = celltype[match(clusters,celltype$ClusterID),'celltype']

P9 <- DimPlot(sce1, reduction = "umap", group.by = "celltype")
P9

plotScoreHeatmap(pred.hesc)
# 基于 per-cell “deltas”诊断，Delta值低，说明注释结果不是很明确
plotDeltaDistribution(pred.hesc, ncol = 3)

write.csv(celltype,file = "D:/A-BS/GSE182109/result/182109_Myeloidcell_sce_annoTable.csv",row.names = F,quote = F)

######################
load("D:/A-BS/GSE182109/result/182109_Tcell_sce_not_annotation.Rdata")
Idents(scRNA1)<-scRNA1$RNA_snn_res.0.9

cells_to_pick <- WhichCells(scRNA1, idents = c(0:9,13:16))

# 获取这些细胞的数据
scRNA1<- subset(scRNA1, cells = cells_to_pick)

library(ggplot2)

celltype<-read_xlsx("D:/A-BS/marker/sce/sce_T_subtype.xlsx",sheet = 2)
clusters<-levels(Idents(sce))

sce@meta.data$celltype = celltype[match(clusters,celltype$ClusterID),'celltype']

sce@meta.data$celltype<-sce@meta.data[["celltype"]][["celltype"]]
sce$celltype<-Idents(sce)
Idents(sce)<-sce$celltype

list_genes<-list(
  DC=c("HLA-DRA","HLA-DPB1","HLA-DPA1"),
  Mac=c("CSF1R","VSIG4","SIGLEC1"),
  Neut=c("IL1R2","ITGAM","FPR2")

)

list_genes<-list(
  CD8_T_EX=c("CD3D","CD3E","CD8A","CD8B","LAG3","PDCD1","TOX"),
  Naive_T=c("IL7R","CCR7","LEF1"),
  T_EFF=c("CCR5","CXCR3"),
  Treg_T=c("TNFRSF4","BATF","TNFRSF18","FOXP3","IL2RA")
)

dot <- DotPlot(sce, features = list_genes)$data
dot$id <- factor(dot$id,levels = unique(as.character(dot$id))[order(unique(as.character(dot$id)))])
#dot <- dot[order(dot$id),]
####做注释文件
data.anno <- data.frame(
  features.plot = rev(unique(dot$features.plot)),
  label = rev(c(rep("CD8_T_EX",7),rep("Naive_T",3),
                rep("T_EFF",2),
                rep("Treg_T",5)
                
  )
  )
)

data.anno <- data.frame(
  features.plot = rev(unique(dot$features.plot)),
  label = rev(c(rep("DC",3),rep("Mac",3),
                rep("Neut",3)
  )
  )
)

df.plot <- plyr::join(dot,data.anno)
df.plot$features.plot <- factor(df.plot$features.plot,levels = rev(levels(df.plot$features.plot)))

p <- ggplot(df.plot,aes(x=features.plot,y =  id,size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("% detected", range = c(0,8)) + #调整绘图点的相对大小
  scale_color_gradient2(low = "darkgrey",mid = "white",high = "red",midpoint = 0,
                        name ="Average\nexpression" 
  ) +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id))+coord_flip()+#一开始设定的x轴和y轴错误，现在调转
  facet_grid(~id, scales="free_x",space = "free")+#theme_classic() +
  theme(
    axis.text.x = element_text(size=8, angle=0, hjust=0.5, color="black",face="bold"),#x轴标签样式
    #axis.text.y = element_blank(),####去掉y轴刻度
    axis.title.x = element_text(size=20,colour = 'black',vjust = -0.8,hjust = 0.5),#坐标轴标题
    axis.ticks.y = element_blank(),#坐标轴刻度
    axis.text.y.right = element_blank(),#坐标轴标签隐藏
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',linewidth = 0.2), #坐标轴轴线样式
    
    panel.spacing=unit(1, "mm"), #分面图图块间距
    strip.text.x = element_text(size=8, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),#分面标签样式
    strip.background.x = element_rect(linewidth = 1)
  )+scale_y_discrete(breaks=NULL)

p
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))

cols<-c("#DF9E9B","#99BADF","#D8E7CA","#99CDCE","#999ACD")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
pdf(file = "D:/A-BS/TCell_dotplot.pdf",width = 7,height = 8)
plot(g)
dev.off()

save(scRNA1,file = "D:/A-BS/GSE182109/result/182109_Tcell_sce_annotation_over.Rdata")




