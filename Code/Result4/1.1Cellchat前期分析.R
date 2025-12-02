#library(devtools)
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(Seurat)
library(dplyr)
library(future)

load("D:/A-BS/CellChat/182109_sce_chat_over.Rdata")

sce_chat<-ScaleData(sce_chat)
sce_chat <- FindVariableFeatures(sce_chat, selection.method = "vst", nfeatures = 2000)

sce_chat<- RunPCA(sce_chat,features = VariableFeatures(object = sce_chat))
#绘制肘部图可视化确定聚类个数
ElbowPlot(scRNA)

#准确确定线性降维的维数
sce=sce_chat
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
#UMAP
sce<- RunUMAP(sce, dims = 1:14)
save(sce,file="D:/A-BS/CellChat/182109_sce_chat_umap_over.Rdata")
###############################
##############
#cellchat
###########
load("D:/A-BS/CellChat/182109_sce_chat_umap_over.Rdata")

cellchat <- createCellChat(object = sce,
                           meta = sce@meta.data,
                           group.by = "celltype")
#如果想用全部的用于cellchat分析
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB 
##################
##取出相应分类用作分析数据库
# 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.human
#查看数据库的组成比例
showDatabaseCategory(CellChatDB)
# 查看数据库具体信息
#CellChatDB$interaction[1:4,1:4]
#head(CellChatDB$cofactor)
#head(CellChatDB$complex)
#head(CellChatDB$geneInfo)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
CellChatDB.use <- CellChatDB
#也可以不取出分类
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中，相当于设置好接下来要参考的数据库
##################

# This step is necessary even if using the whole database
cellchat <- subsetData(cellchat) 
# do parallel ，根据配置设置
#future::plan("multiprocess", workers = 1) 

#识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)

#project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

cellchat@data.project[1:4,1:4]

cellchat <- computeCommunProb(cellchat)
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)


#cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#cellchat <- filterCommunication(cellchat, min.cells = 10)
#cellchat <- computeCommunProbPathway(cellchat)
#cellchat <- aggregateNet(cellchat)

save(cellchat,file = "D:/A-BS/CellChat/Tumor&subtype_Cellchat.Rdata")
#可视化
###############
load("D:/A-BS/CellChat/Tumor&subtype_Cellchat.Rdata")

groupSize <- table(cellchat@idents) %>% as.numeric()

par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Interaction weights/strength")

################
groupSize <- table(cellchat@idents) %>% as.numeric()

par(mfrow = c(4, 3), xpd = TRUE)

for(i in 1:nrow(cellchat@net$weight)){
  mat <- matrix(0, 
                nrow = nrow(cellchat@net$weight),
                ncol = ncol(cellchat@net$weight),
                dimnames = dimnames(cellchat@net$weight))
  mat[i, ] <- cellchat@net$weight[i, ]
  netVisual_circle(mat, 
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(cellchat@net$weight),
                   title.name = rownames(mat)[i])
}
##信号通路层次的可视化
cellchat@netP[["pathways"]]
levels(cellchat@idents)
vertex.receiver<-c(1:3)
#nlevels(sce)=10 vertex.receiver<-c(1:5)
vertex.receiver<-seq(1,5)
c('MIF', 'PTN', 'CD99', 'MK', 'THBS', 'CDH','PDGF','NCAM','TNF')
pathways.show <- "CD45"
netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    vertex.receiver = c(5,4,6,8,9),
                    layout = 'hierarchy')

#
netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    vertex.receiver = vertex.receiver,
                    layout = 'chord')


list<-cellchat@netP[["pathways"]]
for (i in 1:length(list)) {
  p <-netVisual_aggregate(cellchat, 
                          signaling = list[i], 
                          vertex.receiver = c(5,4,6,8,9),
                          layout = 'hierarchy')
  d <- paste("D:/A-BS/CellChat/pathway_new/",list[i],".pdf",sep="")
  pdf(file=d,family = "Times",width=12,height = 10)
  print(p)
  dev.off()
}


###################
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netVisual_chord_gene(cellchat, sources.use = c(5,9), targets.use = c(1,2,3,7,10), lab.cex = 0.5,legend.pos.y = 12)

netVisual_chord_gene(cellchat, sources.use = c(5,9), targets.use = c(1,2,3,7,10), slot.name = "netP", legend.pos.x = 10)

netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    sources.use = c(5,9), #naive CD4 T, Memory CD4 T
                    targets.use = c(2:4,6:8), #CD8 T, FCGR3A+ Mono
                    layout = 'chord')

plotGeneExpression(cellchat, signaling = "CD99")
##################
netVisual_aggregate(cellchat, 
                    signaling = pathways.show,
                    layout = 'circle')

netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = 'Reds')

##################
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(5,9), #tumor cell
                     targets.use = c(1,2,3,7,10), #T B cell 
                     remove.isolate = FALSE)
ggsave("Tumor_TB_ligand-receptor pair.pdf", p, width = 8, height = 12) 
# save as bubble.pdf
##########################
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8)
ht1 + ht2
# save as TIL/SNA_SignalingPattern.pdf
########################
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
#计算分解成几个因子(pattern)比较合适（这一步运行比较慢 。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
selectK(cellchat, pattern = "outgoing")
# save as TIL/NMF_outgoing_selectK.pdf
nPatterns = 4 # 挑选曲线中第一个出现下降的点
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, 
                                          width = 5, height = 9, font.size = 8)
# save as TIL/NMF_outgoing_comPattern_heatmap.pdf
netAnalysis_dot(cellchat, pattern = "outgoing")
# save as TIL/NMF_outgoing_comPattern_dotplot.pdf
selectK(cellchat, pattern = "incoming") 
# save as TIL/NMF_incoming_selectK.pdf
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, 
                                          width = 5, height = 9, font.size = 8)
# save as TIL/NMF_incoming_comPattern_heatmap.pdf
netAnalysis_river(cellchat, pattern = "incoming")
# save as TIL/NMF_incoming_comPattern_river.pdf
netAnalysis_dot(cellchat, pattern = "incoming")
# save as TIL/NMF_incoming_comPattern_dotplot.pdf
####################################
#sci dotplot
## 1.1挑选合适通路 ----
library(Seurat)
library(CellChat)
library(patchwork)
setwd("D:/A-BS/CellChat/")
outdir <- './results/files/'
outdir2 <- './results/plots/'
# 读取cellchat对象
cellchat <- readRDS('cellchat.rds')
table(cellchat@meta[["celltype"]])#查看有哪些细胞类型
## 挑选细胞类型和关注的通路
pdf(file = paste0(outdir2, 'dotployt_inflammation.pdf'),width=20, height= 5)
sources = c(4,8)#设置sources cells对应的位置顺序
targets = c(1,2,6,9)#设置target cells的
p <- netVisual_bubble(cellchat, sources.use = sources, targets.use = targets,  
                      signaling = c('MIF', 'PTN', 'CD99', 'MK', 'THBS', 'CDH','PDGF','NCAM','TNF'),#关注的通路
                      angle.x = 90)
p
dev.off()
## 1.2提取数据 ----
# 提取画图数据
LR_data <- p[["data"]][,1:4]
LR_data <- LR_data[!is.na(LR_data$ligand), ] 
table(LR_data$ligand)
table(LR_data$receptor)
# 更改信息，复杂的受体只保留后者
LR_data$receptor <- gsub('CD74_CD44', 'CD44', LR_data$receptor)
LR_data$receptor <- gsub('CD74_CXCR4', 'CXCR4', LR_data$receptor)
LR_data$receptor <- gsub('ITGA4_ITGB1', 'ITGB1', LR_data$receptor)
#LR_data$receptor <- gsub('ITGAM_ITGB2', 'ITGB2', LR_data$receptor)
#LR_data$receptor <- gsub('ITGAX_ITGB2', 'ITGB2', LR_data$receptor)
table(LR_data$ligand)
table(LR_data$receptor)
write.table(LR_data, file = file.path(outdir, "LR_data.csv"),
            quote = F, sep = ",", row.names = F)
# 整理数据
data <- data.frame(sou=LR_data$ligand,
                   x1=rep(2,length(LR_data$ligand)),
                   net_y1='',
                   tar=LR_data$receptor,
                   x2=rep(3,length(LR_data$receptor)),
                   net_y2='')
library(readxl)
LR_data_name<-read_xlsx('D:/A-BS/CellChat/results/files/LR_data_num.xlsx',sheet=1)
colnames(LR_data_name) <- c('sou', 'number')
# 找到 LR_data_name 中 sou 对应的行索引
match_indices <- match(data$sou, LR_data_name$sou)
# 将 LR_data_name 中对应行的 num 列的值赋给 LR_data 中的 net_y1 列
data$net_y1 <- LR_data_name$number[match_indices]
# 同理
colnames(LR_data_name) <- c('tar', 'number')
match_indices <- match(data$tar, LR_data_name$tar)
data$net_y2 <- LR_data_name$number[match_indices]
data <- unique(data)
write.table(data, file = file.path(outdir, "data_dotplot.csv"),
            quote = F, sep = ",", row.names = T)

## 连线图
p2 <- ggplot(data)+
  geom_segment(aes(x1,net_y1,xend=x2,yend=net_y2),
               size=0.5,color=c('black','black','black','black','#CC0033',
                                'black','black','black','black','black',
                                'black','black','black'))+#虚线 ,linetype = "dashed"
  #设置了红色为了突出显示某个配受体对
  geom_point(aes(x=x1,y=net_y1),size=1, fill="#3cb346", color="#3cb346",
             stroke=1, shape = 21)+
  geom_point(aes(x=x2,y=net_y2),size=1, fill="#44c1f0", color="#44c1f0",
             stroke=1, shape = 21)+
  scale_y_continuous(limits = c(1, 10),expand = expansion(add=c(0.5,0.7)))+
  scale_x_continuous(expand = expansion(0,0.1))+
  theme_void()
p2

ggsave(file.path(outdir2, 'mid_line.pdf'), 
       p2, width = 2, height = 5)
## 3.1左侧
features_1 <- c(LR_data_name$genes[1:9])
p1 <- DotPlot(object = seurat_obj[, seurat_obj[["cell_type"]] == c('NK CD56highCD16low', 'NK CD56lowCD16high', 'CD14+ Monocytes',  'CD14+CD16+ Monocytes', 'CD16+ Monocytes',  'Dendritic cells' )], #可以用全部的细胞类型，也可以挑选一部分
              features = features_1, assay = "RNA") +
  guides(color = guide_colorbar(title = 'Average  Expression')) +
  scale_color_gradientn(colours = c('white','grey',"black"))+#主题颜色 
  coord_flip() +
  theme(axis.title = element_blank(),#标题为空
        legend.position = "left",
        plot.margin = unit(c(0.5,0,0.5,0.5), 'cm'),
        axis.line = element_blank(),#去除坐标轴线
        axis.ticks = element_line(size = 0.2),#刻度线
        panel.border = element_rect(color = "black",fill=NA, size = 0.1),#图形边框
        panel.grid.major = element_line(size = 0.1,color = 'black',linetype = 2),# 显示主网格线
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) +#坐标轴字体
  scale_x_discrete(position = "top") + scale_y_discrete(position = "right")
p1
ggsave(file.path(outdir2, 'left_target.pdf'), 
       p1, width = 4.5, height = 5.5)

## 3.2右侧
features_2 <- c(LR_data_name$genes[10:19])
p3 <- DotPlot(object = seurat_obj, features = features_2, assay = "RNA") +
  guides(color = guide_colorbar(title = 'Average  Expression')) +
  scale_color_gradientn(colours = c('white','grey',"black"))+#主题颜色 
  coord_flip() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0), 'cm'),
        axis.line = element_blank(),#去除坐标轴线
        axis.ticks = element_line(size = 0.2),#刻度线
        panel.border = element_rect(color = "black",fill=NA, size = 0.1),#图形边框
        panel.grid.major = element_line(size = 0.1,color = 'black',linetype = 2),# 显示主网格线
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) +
  scale_y_discrete(position = "right")  + 
  NoLegend()
p3
ggsave(file.path(outdir2, 'right_source.pdf'), 
       p3, width = 7, height = 6)
