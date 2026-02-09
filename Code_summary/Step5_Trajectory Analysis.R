# ==============================================================================
# Step 5: Single-cell Trajectory Analysis (Monocle3)
# Project: GBM TME Analysis
# Description: This script performs pseudotime trajectory analysis using Monocle3,
#              visualizes differentiation paths in 2D/3D, and plots gene expression
#              changes along the trajectory.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(patchwork)
library(magrittr)

# Set working directory (adjust as needed)
# setwd("./") 

# Create output directories
if(!dir.exists("./results/step5_trajectory")) dir.create("./results/step5_trajectory", recursive = TRUE)

# ==============================================================================
# 1. Monocle3 Data Preparation & 2D Trajectory
# Source: 1.1monocle3_2D.R
# ==============================================================================
cat(">>> [Step 5.1] Initializing Monocle3 Analysis (2D)...\n")

# Load Seurat object ( subset T cells or Myeloid cells)

sce_sub <- subset(sce, idents = c("CD8+ T", "Naive T", "Treg"))

# Pre-processing (Retrieve UMAP coordinates from Seurat)
# Monocle3 uses partitions to separate disjoint trajectories
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# --- Define Root Node & Order Cells ---
# Option A: Manual Selection (Interactive - paused for automated script)
# cds <- order_cells(cds)

# Option B: Automated Root Selection (Based on specific cell type, e.g., Naive T)
# Function to find the principal node closest to the median of a cell cluster
get_earliest_principal_node <- function(cds, time_bin = "Naive T") {
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

# Example: Set root to 'T cell' or specific cluster (Adjust 'T cell' to your naive cluster name)
# root_node <- get_earliest_principal_node(cds, time_bin = "T cell")
# cds <- order_cells(cds, root_pr_nodes = root_node)

# If no specific root is known, use the interactive line locally:
# cds <- order_cells(cds) 
# For reproduction script, we assume pseudotime is already calculated or use a default:
cds <- order_cells(cds, root_pr_nodes = NULL) # Will prompt if run interactively, or start from arbitrary root

# Save Monocle Object
save(cds, file = "./results/step5_trajectory/Monocle3_CDS_2D.Rdata")

# --- Visualization (2D) ---
# 1. Trajectory Plot (Colored by Cell Type)
p1 <- plot_cells(cds, 
                 color_cells_by = "celltype", 
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 graph_label_size = 1.5) +
  theme(legend.position = "right") +
  ggtitle("Trajectory by Cell Type")

# 2. Pseudotime Plot
p2 <- plot_cells(cds, 
                 color_cells_by = "pseudotime", 
                 label_cell_groups = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE, 
                 graph_label_size = 1.5) +
  ggtitle("Pseudotime")

# Combine and Save
p_combined <- p1 + p2
print(p_combined)
ggsave("./results/step5_trajectory/Monocle3_2D_Trajectory.pdf", p_combined, width = 14, height = 6)

# ==============================================================================
# 2. 3D Trajectory Analysis
# Source: 1.1monocle3_3D.R
# ==============================================================================
cat(">>> [Step 5.2] Performing 3D Trajectory Analysis...\n")

# Note: 3D plotting requires re-processing dimensions with max_components = 3
cds_3d <- cds

# Reduce dimensions to 3 components
cds_3d <- reduce_dimension(cds_3d, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)

# Order cells (using the same logic as 2D or re-calculating)
# cds_3d <- order_cells(cds_3d) 
cds_3d <- order_cells(cds_3d, root_pr_nodes = NULL) 

# Save 3D CDS
save(cds_3d, file = "./results/step5_trajectory/Monocle3_CDS_3D.Rdata")

# --- Visualization (3D) ---
# Note: plot_cells_3d generates an interactive HTML widget.
# It will not render in a static PDF but can be saved as HTML or viewed in RStudio.
p_3d <- plot_cells_3d(cds_3d, color_cells_by = "celltype")

# To save 3D plot (requires htmlwidgets)
# htmlwidgets::saveWidget(p_3d, file = "./results/step5_trajectory/Trajectory_3D.html")
print("3D Plot generated. View 'p_3d' object or check output HTML.")

# ==============================================================================
# 3. Gene Expression Visualization (Boxplots)
# Source: 1.2sc-gene-箱线图.R
# ==============================================================================
cat(">>> [Step 5.3] Visualizing Gene Expression Trends...\n")

# Define genes of interest (Markers or Risk Genes)
# Replace with your specific genes
my_genes <- c("CD3D", "CD8A", "SPP1", "CD44", "PTN", "MIF") 

# Check if genes exist in the dataset
valid_genes <- my_genes[my_genes %in% rownames(cds)]

if(length(valid_genes) > 0) {
  # Extract expression data and metadata
  # Using normalized counts
  expr_data <- exprs(cds)[valid_genes, , drop = FALSE]
  expr_data <- as.matrix(expr_data)
  expr_df <- as.data.frame(t(expr_data))
  
  # Add metadata (Cell Type and Pseudotime)
  meta <- colData(cds)
  expr_df$celltype <- meta$celltype
  expr_df$pseudotime <- pseudotime(cds, reduction_method = "UMAP")
  
  # Reshape for ggplot (Long format)
  expr_long <- reshape2::melt(expr_df, id.vars = c("celltype", "pseudotime"), 
                              variable.name = "Gene", value.name = "Expression")
  
  # Plot 1: Boxplot by Cell Type
  p_box <- ggplot(expr_long, aes(x = celltype, y = Expression, fill = celltype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.1, alpha = 0.3) +
    facet_wrap(~Gene, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cell Type", y = "Expression Level")
  
  print(p_box)
  ggsave("./results/step5_trajectory/Gene_Expression_Boxplots.pdf", p_box, width = 10, height = 8)
  
  # Plot 2: Expression vs Pseudotime (Scatter + Smooth)
  p_pseudo <- ggplot(expr_long, aes(x = pseudotime, y = Expression, color = celltype)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, color = "black") +
    facet_wrap(~Gene, scales = "free_y") +
    theme_bw() +
    labs(x = "Pseudotime", y = "Expression")
  
  print(p_pseudo)
  ggsave("./results/step5_trajectory/Gene_Pseudotime_Trends.pdf", p_pseudo, width = 10, height = 8)
  
} else {
  warning("None of the specified genes were found in the dataset.")
}

cat(">>> Step 5 Analysis Completed.\n")