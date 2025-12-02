library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
###################
#load
load("D:/A-BS/GSE182109/result/182109_Tcell_sce_annotation_over.Rdata")
expression_matrix <- as(as.matrix(sce@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

#cds
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


###################
###################
cds <- preprocess_cds(cds, num_dim = 10)
plot_pc_variance_explained(cds)
cds <- align_cds(cds, alignment_group = "orig.ident", preprocess_method = "PCA")
cds <- reduce_dimension(cds,reduction_method = 'UMAP',preprocess_method="Aligned",umap.min_dist = 0.5,umap.n_neighbors = 15)


plot_cells(cds, color_cells_by="celltype",cell_size=0.5,group_label_size=5)
plot_cells(cds, color_cells_by="celltype",cell_size=0.5,group_label_size=5)
plot_cells(cds, genes=c("S100a8", "S100a9", "Il1b", "Ly6g"))
# plot_cells(cds, color_cells_by="orig.ident",cell_size=0.5,group_label_size=5)

#ϸ???ִ?
cds <- cluster_cells(cds)
plot_cells(cds,cell_size=0.5,group_label_size=5)

cds <- learn_graph(cds,
                   verbose=T,
                   learn_graph_control=list(minimal_branch_len=20,#default is ten
                                            euclidean_distance_ratio=10#default is one
                   ))

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, celltype="Naive_T"){
  cell_ids <- which(colData(cds)[, "celltype"] == celltype)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

time<-pseudotime(cds)
colData(cds)$pseudotime<-time

p1 <- plot_cells(cds, label_cell_groups = F, 
                 color_cells_by = "pseudotime", 
                 label_leaves = F, 
                 label_branch_points = F, 
                 graph_label_size = 4, 
                 cell_size=0.5, 
                 trajectory_graph_segment_size = 1)
p2 <- plot_cells(cds, color_cells_by="celltype",cell_size=0.5,group_label_size=5)+
  scale_fill_manual(values = c("Naive_T" = "#B699C6", "Treg_T" = "#202D70", "CD8_T_EM" = "#3778AD", "CD8_T_EX" = "#4EA74A"))
p3 <- p1+p2
ggsave()
saveRDS(cds,"cds.rds")

###################

setwd("D:/A-BS/Monocle/")
cds <- readRDS("./cds.rds")
load("monocle3.RData")
cds_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
genes <- row.names(subset(cds_res, q_value< 0.01 & morans_I > 0.1))

plot_matrix <- exprs(cds)[match(genes,
                                rownames(cds)),
                          order(pseudotime(cds))]
plot_matrix
#
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- genes;
colnames(plot_matrix) <- colnames(cds)[order(pseudotime(cds))]
dim(plot_matrix)


callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


p1 <- pheatmap::pheatmap(plot_matrix, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=5,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)


#
annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 5)))
row.names(annotation_row) <- rownames(plot_matrix)

rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3','purple') 
names(rowcolor) <- c("1","2","3","4","5") #??????ɫ
#
cds_meta <- as.data.frame(colData(cds))
cds_meta[colnames(plot_matrix),]$celltype
annotation_col <- data.frame(celltype=factor(cds_meta[colnames(plot_matrix),]$celltype))
row.names(annotation_col) <- colnames(plot_matrix)
colcolor <- c("#3778AD","#B699C6","#4EA74A",'#202D70') 
names(colcolor) <- c("CD8_T_EM","Naive_T","CD8_T_EX","Treg_T") #??????ɫ
#

ann_colors <- list(Cluster=rowcolor,celltype = colcolor) #??ɫ????

p3 <- pheatmap::pheatmap(plot_matrix, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         cutree_rows=5,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         annotation_col = annotation_col,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         main="Pseudotime")
pdf("monocle3heatmap1129.pdf",width = 6,height = 10)
print(p3)
dev.off()

library(ggplot2)
library(ggridges)


annotation_col1 <- annotation_col
annotation_col1$Count <- pseudotime(cds)[order(pseudotime(cds))]  


p <- ggplot(annotation_col1, aes(x = Count, y = factor(celltype, levels = c("CD8_T_EX", "CD8_T_EM", "Treg_T","Naive_T")),
                                                       fill = celltype
                                                       )) +
  geom_density_ridges(scale = 5, rel_min_height = 0.01) +  # scale ????ɽ???Ŀ???
  theme_ridges() +  # ʹ??ɽ??ͼ??????
  theme(legend.position = "none", 
        axis.title.x = element_text(hjust = 0.5)) +  # ????X??????????
  labs(title = "", x = "Pseudotime", y = "") +
  scale_fill_manual(values = c("Naive_T" = "#B699C6", "Treg_T" = "#202D70", "CD8_T_EM" = "#3778AD", "CD8_T_EX" = "#4EA74A"))  # ?Զ?????ɫ
p
ggsave(filename = "ridge.pdf",plot = p,width = 6,height = 6)

p <- ggplot(annotation_col1, aes(x = Count, y = factor(celltype, levels = c("Naive_T", "Treg_T", "CD8_T_EM", "CD8_T_EX")), 
                                 fill = celltype, 
                                 height =after_stat(density))) +
  geom_density_ridges(scale = 1,  # Set scale=1 for equal height ridges
                      stat = "density",  # Ensure density calculation
                      rel_min_height = 0.01) +
  theme_ridges(font_size = 12) +  # 移除网格线会视觉上增强高度感
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5)) +
  labs(x = "Pseudotime", y = "") +
  scale_fill_manual(values = c("Naive_T" = "#B699C6", 
                               "Treg_T" = "#202D70", 
                               "CD8_T_EM" = "#3778AD", 
                               "CD8_T_EX" = "#4EA74A")) +
  scale_y_discrete(limits = rev)  # Reverse to show Naive at top
p
##################################
library(clusterProfiler,lib.loc = "E:\\R-4.4.0\\library")
library(org.Hs.eg.db)

###
module_gene <- as.data.frame(cutree(p3$tree_row, k=5))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)

# module 5
module_gene_final <- module_gene[which(module_gene$Module=="5"),]
write.table(module_gene_final,"sig_gene_TEX.txt",row.names = F,col.names = T,quote = F)

save.image("monocle3.RData")

MG<-module_gene[order(module_gene$Module),]
write_xlsx(MG,"Tmodule_gene.xlsx")

########################
#对每个模块做GO功能富集分析
mg5<- module_gene[which(module_gene$Module=="5"),]
library(clusterProfiler)
library(org.Hs.eg.db)

degs.list=mg5$gene
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
ego_result_BP <- as.data.frame(erich.go.BP)
write.csv(ego_result_BP,file = "GO_result_BP_moduleG5.csv",row.names = T)
