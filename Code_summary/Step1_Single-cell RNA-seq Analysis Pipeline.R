# ==============================================================================
# Step 1: Integrated Single-cell RNA-seq Analysis Pipeline
# Project: GBM TME Analysis
# Covers: QC, Harmony Integration, Annotation, InferCNV, DEG, and GSEA
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(readxl)
library(harmony)
library(clustree)
library(infercnv)
library(RColorBrewer)
library(ggrepel)
library(presto)
library(msigdbr)
library(fgsea)
library(tibble)
library(scRNAtoolVis)

# Set working directory to project root
# setwd("./") 

# ------------------------------------------------------------------------------
# 1. Data Loading and Quality Control (Source: 1.1GSE182109.R)
# ------------------------------------------------------------------------------
cat(">>> [Step 1] Loading Data and Performing QC...\n")

# Define data directory
data_dir <- "./data/GSE182109_RAW_ALL/" 
file_list <- list.files(data_dir)

scRNAlist <- list()
for(i in 1:length(file_list)){
  # Gene column = 2 for 10X data
  counts <- Read10X(data.dir = file.path(data_dir, file_list[i]), gene.column = 2)
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = file_list[i], min.cells = 3, min.features = 200)
  print(paste("Loaded:", file_list[i]))
}

# QC for each sample
for(i in 1:length(scRNAlist)){
  scRNA <- scRNAlist[[i]]
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  # Filter based on features and mitochondrial percentage
  scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  scRNAlist[[i]] <- scRNA
}

# Merge into one Seurat object
scRNA <- merge(scRNAlist[[1]], y = scRNAlist[2:length(scRNAlist)])
print(paste("Total cells:", ncol(scRNA)))
print(paste("Total genes:", nrow(scRNA)))

# Load metadata
# Assuming info.xlsx has columns: GSM, sample_id
info <- read_xlsx("./data/info.xlsx", sheet = 2)
scRNA@meta.data$sample <- info[match(scRNA@meta.data$orig.ident, info$GSM), 2]
# Extract the sample_id column if it's a dataframe
if(is.data.frame(scRNA@meta.data$sample)) {
  scRNA@meta.data$sample <- scRNA@meta.data$sample[[1]]
}

# Save checkpoint
if(!dir.exists("./results")) dir.create("./results")
save(scRNA, file="./results/182109_qc_over.Rdata")

# ------------------------------------------------------------------------------
# 2. Normalization, Integration, and Clustering (Source: 1.1GSE182109.R)
# ------------------------------------------------------------------------------
cat(">>> [Step 2] Normalization, Harmony Integration, and Clustering...\n")

load("./results/182109_qc_over.Rdata")

# Standard workflow
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))

# Elbow Plot
ElbowPlot(scRNA)

# Determine optimal PCs automatically
pct <- scRNA[["pca"]]@stdev / sum(scRNA[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
print(paste("Selected PCs:", pcs))

# Initial Clustering
scRNA <- FindNeighbors(scRNA, dims = 1:16)
scRNA <- FindClusters(scRNA, resolution = seq(0.1, 1.0, 0.1))

# Clustree visualization
# p1 <- clustree(scRNA@meta.data, prefix='RNA_snn_res.') + coord_flip()
# print(p1)

# UMAP
scRNA <- RunUMAP(scRNA, dims = 1:16)

# Harmony Integration (Batch Correction)
Idents(scRNA) <- 'sample'
sce <- scRNA %>% RunHarmony("orig.ident", plot_convergence = TRUE)

# UMAP and Clustering after Harmony
sce <- sce %>% 
  RunUMAP(reduction = "harmony", dims = 1:16) %>%
  FindNeighbors(reduction = "harmony", dims = 1:16) %>%
  FindClusters(resolution = seq(0.1, 1.0, 0.1))

save(sce, file="./results/182109_harmony_over.Rdata")

# ------------------------------------------------------------------------------
# 3. Cell Type Annotation (Source: 1.2 & 1.3)
# ------------------------------------------------------------------------------
cat(">>> [Step 3] Cell Type Annotation...\n")

load("./results/182109_harmony_over.Rdata")

# Select resolution for annotation
Idents(sce) <- "RNA_snn_res.0.2" # Adjust based on your analysis

# Marker Visualization (DotPlot)
list_genes <- list(
  B_cell = c("MZB1","MS4A1","IGHG1","IGHG3","CD79A"),
  Endothelial = c("VWF","ITM2A","CLEC14A","A2M","CLDN5"),
  Glia_Neuron = c("CLU","C1orf61","SOX2","TUBA1A","PTPRZ1","FABP7"),
  Myeloid = c("C1QA","C1QB","C1QC","SPP1","FCER1G"),
  Oligodendrocyte = c("PLLP","MBP","CNP","CLDN11","PLP1"),
  Pericyte = c("RGS5","PDGFRB","NOTCH3","DCN","CD248"),
  T_cell = c("CD2","CD3D","CD3E","CD3G","CD52")
)

dot_data <- DotPlot(sce, features = list_genes)$data

# Plotting logic
p_dot <- DotPlot(sce, features = list_genes) + 
  RotatedAxis() +
  scale_color_gradient2(low = "darkgrey", mid = "white", high = "red") +
  theme_bw() +
  labs(x = "Features", y = "Cluster")
print(p_dot)

# Apply Annotations
# NOTE: Replace the vector below with your actual cluster IDs based on the DotPlot
# e.g., if Cluster 0 is T cell, Cluster 1 is Myeloid...
# new.cluster.ids <- c("T cell", "Myeloid", "Glia", ...)
# Ensure length(new.cluster.ids) matches levels(sce)

# Placeholder: Define your cluster mapping here
current_ids <- levels(sce)
# new.cluster.ids <- vector("character", length = length(current_ids)) 
# names(new.cluster.ids) <- current_ids
# sce <- RenameIdents(sce, new.cluster.ids)

# Add celltype to metadata
sce$celltype <- Idents(sce)

# UMAP Visualization
p_umap <- DimPlot(sce, reduction = 'umap', group.by = 'celltype', label = TRUE, pt.size = 1)
print(p_umap)

# Stacked Bar Plot (Ratio)
Cellratio <- prop.table(table(Idents(sce), sce$sample), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("CellType", "Sample", "Freq")

p_bar <- ggplot(Cellratio, aes(x = Sample, y = Freq, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.7, position = "fill") + 
  theme_classic() +
  labs(x = 'Sample', y = 'Ratio')
print(p_bar)

# Find Markers for verification
scRNA.markers <- FindAllMarkers(sce, logfc.threshold = 0.5, only.pos = TRUE)
scRNA.markers_used <- scRNA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC) %>%
  dplyr::filter(pct.1/pct.2 > 1.5)

# Prepare data for InferCNV (excluding Myeloid for reference usually, or keep specific types)
# Keeping non-myeloid cells for initial InferCNV setup as per original code
sce1 <- subset(sce, idents = c("B cell", "Endothelial cell", "Glia and neuronal cell", 
                               "Oligodendrocyte", "Pericyte", "T cell"))
save(sce1, file="./results/182109_annotation_no_myeloid.Rdata")

# Save metadata for InferCNV
meta_df <- data.frame(row.names = rownames(sce1@meta.data), 
                      celltype = sce1@meta.data$celltype)
write.table(meta_df, "./results/182109_metadata_no_myeloid.txt", sep = '\t', quote = F, col.names = F)

# ------------------------------------------------------------------------------
# 4. InferCNV Analysis (Source: 1.4inferCNV.R)
# ------------------------------------------------------------------------------
cat(">>> [Step 4] Running InferCNV...\n")

# Ensure output directory exists
if(!dir.exists("./results/inferCNV_output")) dir.create("./results/inferCNV_output", recursive = TRUE)

# Prepare Inputs
expFile <- as.data.frame(GetAssayData(sce1, slot = "counts"))
groupFiles <- "./results/182109_metadata_no_myeloid.txt"
geneFile <- "./data/hg38_gencode_v27.txt" # Ensure this file exists in ./data/

# Create InferCNV Object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expFile,
                                     annotations_file = groupFiles,
                                     delim = "\t",
                                     gene_order_file = geneFile,
                                     ref_group_names = c("T cell", "B cell"))

# Run InferCNV (Heavy computation)
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1, 
                              out_dir = "./results/inferCNV_output", 
                              cluster_by_groups = TRUE, 
                              denoise = TRUE,
                              HMM = TRUE,
                              num_threads = 4)

# ------------------------------------------------------------------------------
# 5. InferCNV Scoring & Tumor Classification (Source: 1.5inferCNV_score可视化.R)
# ------------------------------------------------------------------------------
cat(">>> [Step 5] InferCNV Scoring and Tumor/Normal Classification...\n")

# Load InferCNV object
infer_CNV_obj <- readRDS('./results/inferCNV_output/run.final.infercnv_obj')
expr <- infer_CNV_obj@expr.data
meta <- sce1@meta.data

# Calculate thresholds based on reference cells
tmp1 <- expr[, infer_CNV_obj@reference_grouped_cell_indices$`T cell`]
tmp2 <- expr[, infer_CNV_obj@reference_grouped_cell_indices$`B cell`]
tmp <- cbind(tmp1, tmp2)

down <- mean(rowMeans(tmp)) - 2 * mean(apply(tmp, 1, sd))
up <- mean(rowMeans(tmp)) + 2 * mean(apply(tmp, 1, sd))
oneCopy <- up - down

# Define scoring thresholds
a1 <- down - 2 * oneCopy
a2 <- down - 1 * oneCopy
a3 <- up + 1 * oneCopy
a4 <- up + 2 * oneCopy 

# Score Calculation
cnv_score_mat <- as.matrix(expr)
cnv_score_table_pts <- matrix(0, nrow = nrow(cnv_score_mat), ncol = ncol(cnv_score_mat))

cnv_score_table_pts[cnv_score_mat > 0 & cnv_score_mat < a2] <- 2
cnv_score_table_pts[cnv_score_mat >= a2 & cnv_score_mat < down] <- 1
cnv_score_table_pts[cnv_score_mat >= down & cnv_score_mat < up] <- 0
cnv_score_table_pts[cnv_score_mat >= up & cnv_score_mat <= a3] <- 1
cnv_score_table_pts[cnv_score_mat > a3 & cnv_score_mat <= a4] <- 2
cnv_score_table_pts[cnv_score_mat > a4] <- 2

# Total CNV Score
cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"
meta$totalCNV <- cell_scores_CNV[match(rownames(meta), rownames(cell_scores_CNV)), 1]

# Classify Tumor vs Normal based on CNV Score
mean_Value <- mean(meta$totalCNV, na.rm = TRUE)

tumor_meta <- meta[meta$celltype %in% c("Glia and neuronal cell", "Oligodendrocyte"), ]
tumor_cells <- tumor_meta[tumor_meta$totalCNV >= mean_Value, ]
normal_cells <- tumor_meta[tumor_meta$totalCNV < mean_Value, ]

tumor_cells$type <- "Tumor"
normal_cells$type <- "Normal"
immune_meta <- meta[meta$celltype %in% c("T cell", "B cell"), ]
immune_meta$type <- "Immune"

# Combine for downstream analysis
new_meta <- rbind(tumor_cells, normal_cells, immune_meta)
save(new_meta, file="./results/inferCNV_sce_Cluster_new.Rdata")

# Re-cluster Tumor Cells (Visualized by UMAP)
sce_tumor <- subset(sce1, cells = rownames(new_meta))
sce_tumor$type <- new_meta$type
Idents(sce_tumor) <- sce_tumor$type

sce_tumor <- RunPCA(sce_tumor)
sce_tumor <- RunUMAP(sce_tumor, dims = 1:16) # Adjust dims as needed
save(sce_tumor, file="./results/inferCNV_sce_NTI_new.Rdata")

# ------------------------------------------------------------------------------
# 6. Differential Expression (DEG) & Hallmark Analysis (Source: 1.6 & Volcano)
# ------------------------------------------------------------------------------
cat(">>> [Step 6] Differential Expression & GSEA...\n")

load("./results/inferCNV_sce_NTI_new.Rdata")

# Identify Tumor vs Normal
Idents(sce_tumor) <- sce_tumor$type
sce_comp <- subset(sce_tumor, idents = c("Tumor", "Normal"))

# Run Wilcoxon Analysis using Presto
deg_results <- wilcoxauc(sce_comp, 'type')

# Filter for Tumor group
tumor_res <- deg_results %>%
  dplyr::filter(group == "Tumor") %>%
  arrange(desc(auc))

# GSEA - Hallmark Pathways
m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

ranks <- deframe(tumor_res %>% dplyr::select(feature, auc))
fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
save(fgseaRes, file="./results/scRNA_TumorVSN_hallmark.Rdata")

# Volcano Plot
data_vol <- tumor_res
log2FC_cutoff <- 0.1 # Customize threshold
padj_cutoff <- 0.05

data_vol$change <- as.factor(ifelse(data_vol$pval < padj_cutoff & abs(data_vol$logFC) > log2FC_cutoff,
                                    ifelse(data_vol$logFC > log2FC_cutoff, 'UP', 'DOWN'), 'NS'))

p_vol <- ggplot(data_vol, aes(x = -log10(pval), y = logFC)) +
  geom_point(aes(size = -log10(padj), color = -log10(padj))) +
  scale_color_gradientn(colors = c("#4E99C7", "#D1E5F0", "#FCD5BF", "#B2182B")) +
  geom_hline(yintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = 4) +
  theme_bw() +
  labs(title = "Tumor vs Normal DEGs") +
  geom_text_repel(data = subset(data_vol, abs(logFC) >= log2FC_cutoff & padj < padj_cutoff),
                  aes(label = feature), max.overlaps = 15)

print(p_vol)
ggsave("./results/Tumor_vs_Normal_Volcano.pdf", p_vol, width = 8, height = 6)

cat(">>> Step 1 Integrated Analysis Completed.\n")