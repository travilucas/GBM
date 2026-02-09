# ==============================================================================
# Step 7: Spatial Transcriptomics Analysis Pipeline
# Project: GBM TME Analysis
# Description: This script processes 10X Visium spatial transcriptomics data,
#              performs clustering, and applies CARD for spatial deconvolution 
#              using scRNA-seq reference.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(CARD)
library(MuSiC) # CARD dependency
library(ggcorrplot)

# Set working directory (adjust as needed)
# setwd("./") 

# Create output directories
if(!dir.exists("./results/step7_spatial")) dir.create("./results/step7_spatial", recursive = TRUE)

# ==============================================================================
# 1. Spatial Data Loading & Preprocessing
# Source: 1.1 read-STdata.R
# ==============================================================================
cat(">>> [Step 7.1] Loading and Processing Spatial Data...\n")

# Define Data Directory
# Ensure your folder structure matches: ./data/GSE237183/sample_name/
# Containing: filtered_feature_bc_matrix.h5 and spatial/ folder
data_dir <- "./data/GSE237183_Spatial/" 
samples <- list.files(data_dir)

# Example: Load one sample (adjust loop if processing multiple)
# Assuming loading the first sample found or specific GSM ID
sample_name <- samples[1] 
print(paste("Processing sample:", sample_name))

# Load 10X Spatial Data
# slice = "slice1" is standard, image type can be "lowres" or "hires"
se_st <- Load10X_Spatial(
  data.dir = file.path(data_dir, sample_name),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE
)

# Quality Control (Optional visualization)
# plot1 <- VlnPlot(se_st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# plot2 <- SpatialFeaturePlot(se_st, features = "nCount_Spatial") + theme(legend.position = "right")
# print(plot1 | plot2)

# Normalization (SCTransform)
se_st <- SCTransform(se_st, assay = "Spatial", verbose = FALSE)

# Dimensionality Reduction
se_st <- RunPCA(se_st, assay = "SCT", verbose = FALSE)
se_st <- RunUMAP(se_st, reduction = "pca", dims = 1:30)

# Clustering
se_st <- FindNeighbors(se_st, reduction = "pca", dims = 1:30)
se_st <- FindClusters(se_st, verbose = FALSE, resolution = 0.5)

# Visualization
p_dim <- DimPlot(se_st, reduction = "umap", label = TRUE)
p_spatial <- SpatialDimPlot(se_st, label = TRUE, label.size = 3)

pdf("./results/step7_spatial/Spatial_Clustering.pdf", width = 12, height = 6)
print(p_dim + p_spatial)
dev.off()

# Manual Annotation (Placeholder for user input)
# Replace 'new.cluster.ids' with actual names based on marker genes or histology
# current.cluster.ids <- levels(se_st)
# new.cluster.ids <- c("Tumor Core", "Invasive Margin", "Vascular", ...)
# names(new.cluster.ids) <- levels(se_st)
# se_st <- RenameIdents(se_st, new.cluster.ids)

save(se_st, file = "./results/step7_spatial/Seurat_Spatial_Processed.Rdata")

# ==============================================================================
# 2. CARD Deconvolution
# Source: 1.2CARD.R
# ==============================================================================
cat(">>> [Step 7.2] Performing CARD Deconvolution...\n")

# Load scRNA-seq Reference (from Step 1)
# Ensure this file exists
load("./results/182109_annotation_no_myeloid.Rdata") # Variable: sce1 (or sce)

# Prepare scRNA-seq Input for CARD
# 1. Counts Matrix (Sparse)
sc_count <- GetAssayData(sce1, slot = "counts")
# 2. Metadata (Cell Type and Sample)
sc_meta <- sce1@meta.data[, c("celltype", "sample")]
colnames(sc_meta) <- c("cellType", "sampleID")

# Prepare Spatial Input for CARD
# 1. Spatial Counts
spatial_count <- GetAssayData(se_st, slot = "counts", assay = "Spatial")
# 2. Spatial Locations (x, y coordinates)
spatial_location <- GetTissueCoordinates(se_st)
colnames(spatial_location) <- c("x", "y")

# Create CARD Object
# Note: minCountGene and minCountSpot can be adjusted
CARD_obj <- createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleID",
  minCountGene = 100,
  minCountSpot = 5
)

# Run Deconvolution
# This step involves optimizing the spatial correlation, might take time
CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)

# Save CARD Results
save(CARD_obj, file = "./results/step7_spatial/CARD_Deconvolution_Result.Rdata")

# ==============================================================================
# 3. CARD Visualization
# Source: 1.2CARD.R
# ==============================================================================
cat(">>> [Step 7.3] Visualizing Deconvolution Results...\n")

# A. Scatter Pie Plot (Cell Type Proportions at each spot)
# colors: Define a color vector matching your cell types if needed
# colors <- c("#FF0000", "#00FF00", ...) 

p_pie <- CARD.visualize.pie(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location, 
  # colors = colors, 
  radius = 0.5 # Adjust radius based on spot density
)

pdf("./results/step7_spatial/CARD_ScatterPie.pdf", width = 8, height = 8)
print(p_pie)
dev.off()

# B. Spatial Distribution of Specific Cell Types
# Extract cell types
cell_types <- colnames(CARD_obj@Proportion_CARD)

# Plot each cell type
pdf("./results/step7_spatial/CARD_CellType_Distribution.pdf", width = 10, height = 8)
for(ct in cell_types){
  p_map <- CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location, 
    ct.visualize = ct,
    cmap = "YlOrRd",
    pointSize = 2.0
  )
  print(p_map + ggtitle(paste("Spatial Distribution:", ct)))
}
dev.off()

# C. Correlation between Cell Types (Spatial Co-occurrence)
prop_matrix <- CARD_obj@Proportion_CARD
cor_matrix <- cor(prop_matrix)

pdf("./results/step7_spatial/CellType_Correlation.pdf")
print(ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower", lab = TRUE))
dev.off()

cat(">>> Step 7 Spatial Analysis Completed.\n")