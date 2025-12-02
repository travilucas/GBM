library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#3d 安包
#devtools::install_github("t-carroll/monocle3")

get_earliest_principal_node <- function(cds, time_bin="Naive_T"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

load("D:/A-BS/GSE182109/result/182109_Tcell_sce_annotation_over.Rdata")
#######
expression_matrix <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_data = expression_matrix,
                         cell_metadata =  cell_metadata,
                         gene_metadata = gene_annotation)
# 确定维度数
# 数据预处理=标准化+归一化+PCA
cds <- preprocess_cds(cds,num_dim = 15)
#plot_pc_variance_explained(cds)
#3.去批次----
# 去批次(逐个样本就不用了) = Runharmony
cds <- align_cds(cds, alignment_group = "orig.ident", preprocess_method = "PCA")

cds_3d <- reduce_dimension(cds, max_components = 3,reduction_method = 'UMAP',
                             preprocess_method="Aligned",umap.min_dist =0.5, 
                             umap.n_neighbors = 15)
plot_cells_3d(cds_3d, color_cells_by="celltype",
              color_palette = c("#3778AD","#B699C6","#4EA74A","#202D70"))


cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d,verbose=T,
                      learn_graph_control=list(minimal_branch_len=20,#default is ten
                                               euclidean_distance_ratio=10#default is one
                      )
                      )

cds_3d <- order_cells(cds_3d,root_pr_nodes=get_earliest_principal_node(cds_3d))
time<-pseudotime(cds_3d)
colData(cds_3d)$pseudotime<-time
#cds_3d <- order_cells(cds_3d)
plot_cells_3d(cds_3d, color_cells_by="celltype",
                show_trajectory_graph = F)
plot_cells_3d(cds_3d, color_cells_by="pseudotime",show_trajectory_graph = TRUE) 
