load("D:/A-BS/CellChat/Tumor&subtype_Cellchat.Rdata")
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(Seurat)
library(dplyr)
library(future)
library(ReactomePA)
# 执行富集分析
gene<-cellchat@data.signaling@Dimnames[[1]]

library(clusterProfiler)
library("org.Hs.eg.db")

ligand<-unique(cellchat@LR[["LRsig"]][["ligand"]])
receptor<-unique(cellchat@LR[["LRsig"]][["receptor"]])
gene<-c(ligand,receptor)

gene1 = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

pathway_enrichment <- enrichPathway(gene1$ENTREZID,organism = "human",readable=TRUE)

# 查看富集分析结果
head(pathway_enrichment)
# 使用CellChat的可视化函数绘制富集结果
library(enrichplot)
barplot(pathway_enrichment)
#####################################
p1<-cnetplot(pathway_enrichment,layout="kk",node_label="all",shadowtext="all",
             color.params=list(foldChange=NULL,edge=FALSE,category="#E5C494",gene="#B3B3B3"),
             cex.params=list(category_node=1,gene_node=1,category_label=1,gene_label=1),
             hilight.params=list(category=NULL,alpha_hilight=1,alpha_no_hilight=0.3)
             )
p1
####################################

# 绘制细胞间的通信网络图
netVisual(cellchat, type = "network")
# 绘制通路的可视化图
pathway_plot <- netVisual_aggregate(cellchat, signaling = "all")


###########################
# 提取通信强度数据
net.weight <- cellchat@net[["weight"]]
# 设定一个阈值
threshold <- 0.5
# 筛选强烈的通信
strong_interactions <- net.weight[net.weight > threshold]
library(ComplexHeatmap)
net_weight_matrix_dense <- as.matrix(net.weight)
library(RColorBrewer)
library(circlize)
# 创建热图
Heatmap(matrix=net_weight_matrix_dense, 
        name = "Communication Weight", 
        col = colorRamp2(c(min(net_weight_matrix_dense), 0, max(net_weight_matrix_dense)), 
                         c("#88CEE6", "white", "#C25759")), # 根据数据范围选择颜色
        clustering_distance_rows = "euclidean", # 行聚类距离
        clustering_distance_columns = "euclidean", # 列聚类距离
      
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_side = "left",
        column_names_side = "top",
        cell_fun = function(j, i, x, y, w, h, col) { # 自定义单元格内容
          grid.text(sprintf("%.2f", net_weight_matrix_dense[i, j]), x, y, gp = gpar(fontsize = 10))
        }
)


Heatmap(matrix=net_weight_matrix_dense, 
        name = "Communication Weight", 
        col = colorRamp2(c(min(net_weight_matrix_dense), 0, max(net_weight_matrix_dense)), 
                         c("#88CEE6", "white", "#C25759")), # 根据数据范围选择颜色
        clustering_distance_rows = "euclidean", # 行聚类距离
        clustering_distance_columns = "euclidean", # 列聚类距离
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_side = "right",
        column_names_side = "bottom",
        
)

# 假设你已经完成了CellChat分析，并且cellchat_obj是CellChat对象
# 获取细胞间的通信概率
communication_prob <- cellchat@net[["weight"]]
# 找出互作最强烈的细胞对
strongest_interactions <- sort(communication_prob, decreasing = TRUE)
# 设定一个阈值，例如top 10%的互作
threshold <- quantile(strongest_interactions, probs = 0.9)

