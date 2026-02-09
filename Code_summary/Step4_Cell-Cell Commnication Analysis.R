# ==============================================================================
# Step 4: Cell-Cell Communication Analysis (CellChat) & Clinical Validation
# Project: GBM TME Analysis
# Description: This script performs cell-cell communication analysis using 
#              CellChat, visualizes signaling networks/patterns, and validates 
#              ligand-receptor pairs using survival analysis (KM curves).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Libraries
# ------------------------------------------------------------------------------
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(Seurat)
library(dplyr)
library(future)
library(NMF)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(survival)
library(survminer)
library(jjAnno)

# Set working directory (adjust as needed)
# setwd("./") 

# Create directories
if(!dir.exists("./results/step4_cellchat")) dir.create("./results/step4_cellchat", recursive = TRUE)

# ==============================================================================
# 1. Data Preparation & CellChat Object Construction
# Source: 1.1Cellchat前期分析.R
# ==============================================================================
cat(">>> [Step 4.1] Initializing CellChat Analysis...\n")

# Load processed Seurat object from Step 1
# Ensure this file exists (e.g., from Step 1 results)
load("./results/182109_annotation_over.Rdata") # Variable: sce or sce_chat

# Prepare Seurat object (if not already prepared)
# Ensure 'celltype' metadata exists
if(!"celltype" %in% colnames(sce@meta.data)){
  sce$celltype <- Idents(sce)
}

# Create CellChat Object
cellchat <- createCellChat(object = sce, group.by = "celltype")

# Set CellChat Database (Human)
CellChatDB <- CellChatDB.human 
# showDatabaseCategory(CellChatDB)

# Use the full database or specific categories (e.g., Secreted Signaling)
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use

# Subsetting data (saving computation)
cellchat <- subsetData(cellchat)

# Parallel processing setup
# future::plan("multisession", workers = 4) 

# ==============================================================================
# 2. Inference of Cell-Cell Communication
# Source: 1.1Cellchat前期分析.R
# ==============================================================================
cat(">>> [Step 4.2] Inferring Communication Networks...\n")

# Identify over-expressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression data onto PPI (Optional but recommended)
# cellchat <- projectData(cellchat, PPI.human)

# Compute communication probability
# raw.use = FALSE if projectData was used
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter out communications with few cells
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer signaling pathways
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network
cellchat <- aggregateNet(cellchat)

# Save intermediate result
save(cellchat, file = "./results/step4_cellchat/Tumor_subtype_Cellchat.Rdata")

# ==============================================================================
# 3. Visualization: Network & Patterns
# Source: 1.1Cellchat前期分析.R & 1.3LR-dotplot.R
# ==============================================================================
cat(">>> [Step 4.3] Visualizing Networks and Patterns...\n")

# A. Aggregate Network (Circle Plot)
groupSize <- table(cellchat@idents) %>% as.numeric()
pdf("./results/step4_cellchat/NetVisual_Circle_Count_Weight.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction weights/strength")
dev.off()

# B. Specific Pathways (Hierarchy/Chord)
# Pathways of interest defined in your code
pathways.show <- c('MIF', 'PTN', 'CD99', 'MK', 'THBS', 'CDH','PDGF','NCAM','TNF', 'CD45')

# Loop to generate plots for each pathway
# Note: vertex.receiver definition depends on your specific clusters. 
# Adjust c(5,4,6,8,9) to match indices of T cells/Myeloid cells in your data
vertex.receiver <- c(1:3) # Example placeholder

for (pathway in pathways.show) {
  # Skip if pathway is not significant
  if(pathway %in% cellchat@netP$pathways) {
    pdf(file = paste0("./results/step4_cellchat/Pathway_", pathway, ".pdf"), width = 12, height = 10)
    tryCatch({
      print(netVisual_aggregate(cellchat, signaling = pathway, layout = 'chord'))
      print(netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = 'Reds'))
    }, error = function(e) { message(paste("Error plotting", pathway)) })
    dev.off()
  }
}

# C. Communication Patterns (NMF)
# Outgoing Patterns
# selectK(cellchat, pattern = "outgoing") # Run manually to check elbow point
nPatterns_out <- 4 # Based on your previous code
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns_out)

pdf("./results/step4_cellchat/NMF_Outgoing_Patterns.pdf", width = 10, height = 8)
print(netAnalysis_dot(cellchat, pattern = "outgoing"))
dev.off()

# Incoming Patterns
# selectK(cellchat, pattern = "incoming")
nPatterns_in <- 3 # Based on your previous code
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns_in)

pdf("./results/step4_cellchat/NMF_Incoming_Patterns.pdf", width = 10, height = 8)
print(netAnalysis_dot(cellchat, pattern = "incoming"))
dev.off()

# ==============================================================================
# 4. Customized Bubble Plot for L-R Pairs
# Source: 1.3LR-dotplot.R
# ==============================================================================
cat(">>> [Step 4.4] Creating Customized Bubble Plots...\n")

# Define specific L-R pairs of interest
TB_target <- data.frame(interaction_name = c("SPP1_CD44", "SPP1_ITGA4_ITGB1", "PTN_PTPRZ1", 
                                             "PTN_NCL", "MIF_CD74_CXCR4", "MIF_CD74_CD44"))
ICI <- data.frame(interaction_name = c("LGALS9_HAVCR2", "CD274_PDCD1", "PDCD1LG2_PDCD1"))
Cytokine <- data.frame(interaction_name = c("IL21_IL13RA2", "IL7_IL4R_IL13RA2"))

# Define colors
cols <- c("#D6251F","#3778AD","#B699C6","#4EA74A","#202D70","#E57FB0","#A15528","#8F4C9A")

# Identify cluster indices for sources and targets
# NOTE: You must verify these indices match your `levels(cellchat@idents)`
# Example logic:
# sources.use = Tumor clusters
# targets.use = T/B cell clusters
# Adjust indices (e.g., c(2,3,7,10)) based on your actual cluster order

# Plotting Bubble Plot
# Combine custom pairs if needed, or plot for all signaling
p_bubble <- netVisual_bubble(cellchat, 
                             sources.use = c(2,3,7,10), # Tumor cells (Adjust indices!)
                             targets.use = c(5,9),      # T/B cells (Adjust indices!)
                             remove.isolate = FALSE,
                             pairLR.use = ICI) +        # Using ICI pairs as example
  theme(legend.position = 'right',
        legend.key.width = unit(1,'cm'),
        plot.margin = unit(c(2.5, 2.5, 2.5, 2.5),'cm')) +
  coord_cartesian(clip = 'off')

# Add annotation segment (requires jjAnno)
# p_bubble <- annoSegment(object = p_bubble, annoPos = 'top', aesGroup = T, aesGroName = 'source', ...)

pdf("./results/step4_cellchat/Bubble_ICI_Pairs.pdf", width = 8, height = 6)
print(p_bubble)
dev.off()

# ==============================================================================
# 5. Functional Enrichment of LR Genes
# Source: 1.2cellchat富集+提取互作.R
# ==============================================================================
cat(">>> [Step 4.5] Performing Functional Enrichment...\n")

# Extract significant Ligands and Receptors
if(!is.null(cellchat@LR[["LRsig"]])){
  ligand <- unique(cellchat@LR[["LRsig"]][["ligand"]])
  receptor <- unique(cellchat@LR[["LRsig"]][["receptor"]])
  genes_lr <- c(ligand, receptor)
  
  # Convert to ENTREZID
  gene_entrez <- bitr(genes_lr, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  # Reactome Pathway Enrichment
  pathway_enrichment <- enrichPathway(gene_entrez$ENTREZID, organism = "human", readable = TRUE)
  
  # Visualization
  pdf("./results/step4_cellchat/LR_Reactome_Enrichment.pdf", width = 10, height = 8)
  print(barplot(pathway_enrichment, showCategory = 20))
  print(cnetplot(pathway_enrichment, layout = "kk", node_label = "all", 
                 color.params = list(category = "#E5C494", gene = "#B3B3B3")))
  dev.off()
}

# ==============================================================================
# 6. Clinical Validation: L-R Pair Survival Analysis
# Source: LR-KM曲线.R
# ==============================================================================
cat(">>> [Step 4.6] Validating L-R Pairs with Survival Analysis...\n")

# Load Bulk Expression Data (from Step 2/3)
load("./data/processed/TCGA_CGGA_input_exp.Rdata") # exp_new
# Or load specific dataset like GSE108474
if(file.exists("./data/GEO/GSE108474exp_over.Rdata")) load("./data/GEO/GSE108474exp_over.Rdata")

# Select Dataset for Validation (e.g., GSE108474 or TCGA)
# Ensure 'exp_final_gse1' or 'exp_new' is active
if(exists("exp_final_gse1")) {
  validation_data <- exp_final_gse1
  dataset_name <- "GSE108474"
} else {
  validation_data <- exp_new
  dataset_name <- "TCGA_CGGA"
}

# Pre-processing
colnames(validation_data) <- trimws(colnames(validation_data))
surv_info <- validation_data[, 1:2] # Assuming first two cols are OS.time and OS

# Define L-R genes to test (e.g., EFNB1 - OPHN1)
genes_to_test <- c("EFNB1", "OPHN1") 

# Check if genes exist in dataset
valid_genes <- genes_to_test[genes_to_test %in% colnames(validation_data)]

if(length(valid_genes) == 2) {
  exp_subset <- validation_data[, valid_genes]
  
  # Calculate Interaction Term (Gene A * Gene B) to represent L-R strength
  exp_subset$interaction <- exp_subset[,1] * exp_subset[,2]
  
  # Combine with survival info
  sur_exp <- cbind(surv_info, exp_subset)
  colnames(sur_exp)[1:2] <- c("OS.time", "OS")
  sur_exp$OS.time <- as.numeric(sur_exp$OS.time)
  sur_exp$OS <- as.numeric(sur_exp$OS)
  
  # Determine Cutoff using surv_cutpoint
  res.cut <- surv_cutpoint(sur_exp, time = "OS.time", event = "OS", variables = "interaction")
  cutpoint <- res.cut[["cutpoint"]][["cutpoint"]]
  
  # Grouping
  sur_exp$risk_group <- ifelse(sur_exp$interaction > cutpoint, "High", "Low")
  
  # Survival Fit
  fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = sur_exp)
  
  # Plot KM Curve
  p_km <- ggsurvplot(
    fit,
    data = sur_exp,
    title = paste0(dataset_name, ": ", genes_to_test[1], "-", genes_to_test[2]),
    palette = c("#E7B800", "#2E9FDF"),
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    legend.labs = c(paste0("High ", paste(genes_to_test, collapse="*")), 
                    paste0("Low ", paste(genes_to_test, collapse="*"))),
    ggtheme = theme_classic2()
  )
  
  pdf(paste0("./results/step4_cellchat/Survival_", genes_to_test[1], "_", genes_to_test[2], ".pdf"))
  print(p_km)
  dev.off()
  
} else {
  warning("Specified genes for survival analysis not found in the dataset.")
}

cat(">>> Step 4 Analysis Completed.\n")