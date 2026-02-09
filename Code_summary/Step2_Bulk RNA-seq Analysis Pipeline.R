# ==============================================================================
# Step 2: Integrated Bulk RNA-seq Analysis Pipeline
# Project: GBM TME Analysis
# Description: This script performs bulk RNA-seq data processing, integration 
#              (TCGA, CGGA, GTEx), batch correction (ComBat), immune-based 
#              consensus clustering, and downstream differential analysis.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Libraries
# ------------------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(stringr)
library(tinyarray)
library(AnnoProbe)
library(GenomicFeatures)
library(clusterProfiler)
library(org.Hs.eg.db)
library(sva)           # For ComBat batch correction
library(limma)         # For DEG analysis
library(GSVA)          # For ssGSEA
library(ConsensusClusterPlus)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(survival)
library(survminer)
library(readxl)
library(msigdbr)
library(fgsea)
library(estimate)      

# Set working directory (adjust as needed)
# setwd("./") 

# Create directories
if(!dir.exists("./results/step2_bulk")) dir.create("./results/step2_bulk", recursive = TRUE)
if(!dir.exists("./data/processed")) dir.create("./data/processed", recursive = TRUE)

# ==============================================================================
# 1. GTEx Data Preprocessing
# Source: 1.1GTEx data预处理.R
# ==============================================================================
cat(">>> [Step 2.1] Processing GTEx Normal Data...\n")

# Load GTEx expected count data
dat <- data.table::fread("./data/gtex_gene_expected_count.gz", data.table = F)
exp <- column_to_rownames(dat, "sample") %>% as.matrix()

# Format rownames (Remove version numbers from ENSEMBL IDs)
rownames(exp) <- str_split(rownames(exp), "\\.", simplify = T)[,1]

# Annotate probes
an <- annoGene(rownames(exp), ID_type = "ENSEMBL")
exp <- trans_array(exp, ids = an, from = "ENSEMBL", to = "SYMBOL")

# Load phenotype data
clinical <- data.table::fread("./data/GTEX_phenotype.gz")
# Filter for Brain samples
clinical <- clinical[clinical$`_primary_site` != "<not provided>", ]
colnames(clinical)[3] <- "site"
clinical_subset <- subset(clinical, site == "Brain")

# Match samples
s <- intersect(colnames(exp), clinical_subset$Sample)
clinical_subset <- clinical_subset[match(s, clinical_subset$Sample), ]
exp <- exp[, s]

# Convert log2(expected_count+1) back to raw counts
exp <- round(2^exp - 1, 4)

# Save processed GTEx data
write.table(data.frame(ID = rownames(exp), exp), 
            file = "./data/GTEX_Brain_Normal_sample.txt", 
            sep = "\t", quote = F, row.names = F)
write.table(clinical_subset, 
            file = "./data/GTEX_Brain_Normal_sample_phenotype.txt", 
            sep = "\t", quote = F, row.names = F)

# Filter for "Brain - Cortex" specific samples
gtex_info <- read.table("./data/GTEX_Brain_Normal_sample_phenotype.txt", sep = '\t', header = T)
gtex_info1 <- gtex_info[which(gtex_info$body_site_detail..SMTSD. == "Brain - Cortex"), ]

gtex_sample <- colnames(exp)
# Format sample IDs
rownames(gtex_info1) <- gsub('[-]', '.', gtex_info1$Sample)
gtex_sample <- gtex_sample[which(gtex_sample %in% rownames(gtex_info1))]

gtex_exp <- exp[, which(colnames(exp) %in% gtex_sample)]
save(gtex_exp, file = "./data/processed/GTEx_Normal105.Rdata")

# ==============================================================================
# 2. TPM Calculation and Merging (TCGA + CGGA + GTEx)
# Source: 1.2bulk_TPM 合并样本.R & 1.3bulk_TPMexp_log_merge.R
# ==============================================================================
cat(">>> [Step 2.2] Calculating TPM and Merging Datasets...\n")

# --- Helper Function: Counts to TPM ---
# Requires effective length from GTF
# txdb <- makeTxDbFromGFF("./data/gencode.v36.annotation.gtf", format="gtf")
# exons.list.per.gene <- exonsBy(txdb, by="gene")
# exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))
# gfe <- data.frame(gene_id=names(exonic.gene.sizes), length=exonic.gene.sizes)
# save(gfe, file = "./data/gfe_length.Rdata")

load("./data/gfe_length.Rdata")

Counts2TPM <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

# --- Process TCGA (Example for LGGGBM) ---
load("./data/TCGA-GBM-Origion_count_Tumor2.Rdata") # exp_new1
# Map Symbols to ENSEMBL
gene_map <- bitr(rownames(exp_new1), fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
exp_new1 <- exp_new1[which(rownames(exp_new1) %in% gene_map$SYMBOL), ]
gene_map <- gene_map[match(rownames(exp_new1), gene_map$SYMBOL), ]
rownames(exp_new1) <- gene_map$ENSEMBL

# Match with gene lengths
gfe$ENSG <- sub("\\..*", "", gfe$gene_id)
common_genes <- intersect(rownames(exp_new1), gfe$ENSG)
exp_TCGA <- exp_new1[common_genes, ]
gfe_new <- gfe[match(common_genes, gfe$ENSG), ]
effLen <- gfe_new$length

# Calculate TPM
TCGA_TPMs <- apply(exp_TCGA, 2, Counts2TPM, effLen = effLen)
TCGA_TPMs <- as.data.frame(TCGA_TPMs)
save(TCGA_TPMs, file = "./data/processed/TCGA_GBMexp_Tumor_TPM.Rdata")

# --- Process GTEx TPM ---
GTEx_TPMs <- apply(gtex_exp[common_genes, ], 2, Counts2TPM, effLen = effLen)
save(GTEx_TPMs, file = "./data/processed/GTEx_Normal105_TPM.Rdata")

# --- Process CGGA (FPKM to TPM) ---
exp_cgga <- read.table("./data/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", header = T, sep = "\t", row.names = 1)
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
CGGA_TPMs <- apply(exp_cgga, 2, FPKM2TPM)
save(CGGA_TPMs, file = "./data/processed/CGGA_693_Tumor_TPM.Rdata")

# --- Merge and Log Transform ---
cat(">>> Merging TCGA, CGGA, and GTEx...\n")

# Load all TPMs
load("./data/processed/TCGA_GBMexp_Tumor_TPM.Rdata") # LG_TPMs
load("./data/processed/CGGA_693_Tumor_TPM.Rdata")   # TPMS693
load("./data/processed/GTEx_Normal105_TPM.Rdata")   # gtex

# Log2 Transform
LG_TPM_log <- log2(LG_TPMs + 1)
CGGA_TPM_log <- log2(TPMS693 + 1)
GTEx_TPM_log <- log2(gtex + 1)

# Map ENSEMBL back to SYMBOL for merging
trans_to_symbol <- function(data) {
  # Remove version from ENSEMBL
  rownames(data) <- gsub("\\..*", "", rownames(data))
  gene_map <- bitr(rownames(data), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  data <- data[rownames(data) %in% gene_map$ENSEMBL, ]
  gene_map <- gene_map[match(rownames(data), gene_map$ENSEMBL), ]
  
  # Handle duplicates by taking max mean expression
  data$Symbol <- gene_map$SYMBOL
  data <- data %>% group_by(Symbol) %>% summarize_all(mean) %>% column_to_rownames("Symbol")
  return(data)
}

TCGA_merge <- trans_to_symbol(LG_TPM_log)
GTEx_merge <- trans_to_symbol(GTEx_TPM_log)
# CGGA is usually already in Symbol, if not apply similar logic

# Intersect genes
common_genes <- intersect(rownames(TCGA_merge), rownames(CGGA_TPM_log))
common_genes <- intersect(common_genes, rownames(GTEx_merge))

Tumor_merge <- cbind(TCGA_merge[common_genes, ], CGGA_TPM_log[common_genes, ])
Normal_merge <- GTEx_merge[common_genes, ] # Add TCGA normals if available

save(Tumor_merge, file = "./data/processed/GBM_Tumor_merge_over.Rdata")
save(Normal_merge, file = "./data/processed/GBM_Normal_merge_over.Rdata")

# ==============================================================================
# 3. Clinical Info Merging
# Source: 1.4合并info信息.R
# ==============================================================================
cat(">>> [Step 2.3] Merging Clinical Information...\n")

# Load clinical files
LG_info <- read.xlsx("./data/TCGA_LGGGBM_ALLinfo2.xls", sheetIndex = 1)
GBM_info <- read.xlsx("./data/TCGA_LGGGBM_ALLinfo_over.xlsx", sheetIndex = 1)
CGGA325_info <- read.xlsx("./data/CGGA325_GBM_info.xlsx", sheetIndex = 1)
CGGA693_info <- read.xlsx("./data/CGGA693_GBM_info.xlsx", sheetIndex = 1)

# Standardize columns
colnames(LG_info) <- c("Sample", "Age", "Gender", "OS", "OS.time", "Platform")
colnames(GBM_info) <- c("Sample", "Age", "Gender", "OS", "OS.time", "Platform")
# ... (Standardize CGGA columns similarly)

# Merge
Tumor_info <- bind_rows(LG_info, GBM_info, CGGA325_info, CGGA693_info)
Tumor_info$age1 <- ifelse(Tumor_info$Age < 60, "<60", ">=60")
Tumor_info$Gender <- str_to_title(Tumor_info$Gender)

write.xlsx(Tumor_info, file = "./results/step2_bulk/GBM-Tumor_info.xlsx")

# ==============================================================================
# 4. Batch Correction (ComBat) & Differential Analysis
# Source: 1.5Tumor-Normal差异表达分析.R
# ==============================================================================
cat(">>> [Step 2.4] Batch Correction and Tumor vs Normal Analysis...\n")

load("./data/processed/GBM_Tumor_merge_over.Rdata")
load("./data/processed/GBM_Normal_merge_over.Rdata")

# 1. ComBat for Tumor
# Define batch vector based on sample origin
batch_tumor <- c(rep("TCGA", ncol(TCGA_merge)), rep("CGGA", ncol(CGGA_TPM_log)))
tumor_combat <- ComBat(dat = as.matrix(Tumor_merge), batch = batch_tumor)
save(tumor_combat, file = "./data/processed/GBM_Tumor_merge_Combatover.Rdata")

# 2. ComBat for Normal (if combining TCGA normal and GTEx)
# batch_normal <- c(...)
# normal_combat <- ComBat(...)

# 3. Merge Tumor and Normal for DEG
# Ensure matching genes
common <- intersect(rownames(tumor_combat), rownames(Normal_merge))
exp_final <- cbind(tumor_combat[common, ], Normal_merge[common, ])
group_list <- c(rep("Tumor", ncol(tumor_combat)), rep("Normal", ncol(Normal_merge)))

# 4. Limma Analysis
design <- model.matrix(~0 + factor(group_list))
colnames(design) <- c("Normal", "Tumor")
rownames(design) <- colnames(exp_final)

contrast.matrix <- makeContrasts("Tumor-Normal", levels = design)
fit <- lmFit(exp_final, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

nrDEG <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(nrDEG)

write.table(nrDEG, "./results/step2_bulk/diff.ControlVSTumor.txt", sep = '\t', quote = F)
save(nrDEG, file = "./results/step2_bulk/nrDEGoutput.Rdata")

# 5. Volcano Plot
DEG <- nrDEG
logFC_cutoff <- 1
p_cutoff <- 0.05
DEG$change <- as.factor(ifelse(DEG$adj.P.Val < p_cutoff & abs(DEG$logFC) > logFC_cutoff,
                               ifelse(DEG$logFC > logFC_cutoff, 'UP', 'DOWN'), 'NOT'))

p_vol <- ggplot(data = DEG, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c('#71C6D4', 'grey', '#E31A1C')) +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10 adj.P.Val", title = "Tumor vs Normal") +
  geom_hline(yintercept = -log10(p_cutoff), lty = 4) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 4)

print(p_vol)
ggsave("./results/step2_bulk/Volcano_TVN.pdf", p_vol)

# ==============================================================================
# 5. Consensus Clustering (ssGSEA)
# Source: 1.6scMarkers_ssGSEA_select_samples.R
# ==============================================================================
cat(">>> [Step 2.5] ssGSEA and Consensus Clustering...\n")

# Load Immune Markers
im <- read.xlsx("./data/immune_markers_over.xlsx", sheetIndex = 1)
im_list <- as.list(im)
# Remove NAs from list
im_list <- lapply(im_list, function(x) x[!is.na(x)])

# Run ssGSEA
load("./data/processed/GBM_Tumor_merge_Combatover.Rdata")
gsvaP <- ssgseaParam(exprData = as.matrix(tumor_combat), geneSets = im_list, 
                     minSize = 1, maxSize = Inf, alpha = 0.25, normalize = T)
ssGSEA_matrix <- gsva(gsvaP)
save(ssGSEA_matrix, file = "./results/step2_bulk/TCGA_tumor_ssGSEA_immune.Rdata")

# Consensus Clustering
title_dir <- "./results/step2_bulk/ConsensusCluster/"
results <- ConsensusClusterPlus(as.matrix(ssGSEA_matrix), maxK = 9, reps = 1000, pItem = 0.8,
                                pFeature = 0.8, title = title_dir, clusterAlg = "hc", 
                                distance = "pearson", seed = 1262118388.71279, plot = "pdf")

# Extract Cluster 2 Results (Sample Clustering)
annCol <- data.frame(results = paste0("Cluster", results[[2]][['consensusClass']]),
                     row.names = colnames(ssGSEA_matrix))
write.xlsx(annCol, file = "./results/step2_bulk/tumorSample_cluster_results.xlsx")

# ==============================================================================
# 6. Cluster Characterization (Heatmap & Validation)
# Source: 1.7info预处理+heatmap.R & cibersort+estimate验证.R
# ==============================================================================
cat(">>> [Step 2.6] Cluster Visualization...\n")

# Heatmap
load("./results/step2_bulk/TCGA_tumor_ssGSEA_immune.Rdata")
annCol <- read.xlsx("./results/step2_bulk/tumorSample_cluster_results.xlsx", sheetIndex = 1)
rownames(annCol) <- annCol$row.names

# Load Clinical Info to Add Annotations
info <- read.xlsx("./results/step2_bulk/GBM-Tumor_info.xlsx", sheetIndex = 1)
info <- info[match(rownames(annCol), info$Sample), ]

annotation_col <- data.frame(
  Cluster = annCol$results,
  Gender = info$Gender,
  Status = ifelse(info$OS == 0, "Alive", "Dead"),
  Age = info$age1
)
rownames(annotation_col) <- rownames(annCol)

ann_colors <- list(
  Cluster = c("Cluster1" = "#EEA599", "Cluster2" = "#92B4C8"),
  Gender = c("Male" = "#D9BDD8", "Female" = "#458DA3"),
  Status = c("Alive" = "#93B237", "Dead" = "#F6C957"),
  Age = c("<60" = "#FFB77F", ">=60" = "#8086bb")
)

pheatmap(as.matrix(ssGSEA_matrix),
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#6892c5", "white", "#b34557"))(100),
         filename = "./results/step2_bulk/Cluster_Heatmap.pdf")

# ESTIMATE Validation
filterCommonGenes(input.f = "./data/processed/GBM_Tumor_merge_Combatover.Rdata", # Needs conversion to txt
                  output.f = "./results/step2_bulk/estimate_input.gct", 
                  id = "GeneSymbol")
estimateScore(input.ds = "./results/step2_bulk/estimate_input.gct",
              output.ds = "./results/step2_bulk/estimate_score.gct", 
              platform = "affymetrix")

# ==============================================================================
# 7. Subtype DEG & GSEA
# Source: 1.8TCGA_hot-cold DEG.R & 1.9bulk_GSEA.R
# ==============================================================================
cat(">>> [Step 2.7] Subtype DEG and GSEA...\n")

# DEG between Clusters
exp <- tumor_combat[, rownames(annotation_col)]
design <- model.matrix(~0 + factor(annotation_col$Cluster))
colnames(design) <- levels(factor(annotation_col$Cluster))
rownames(design) <- colnames(exp)

contrast.matrix <- makeContrasts("Cluster1-Cluster2", levels = design)
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
nrDEG_subtype <- topTable(fit2, coef = 1, n = Inf)
save(nrDEG_subtype, file = "./results/step2_bulk/nrDEG_Subtype.Rdata")

# GSEA Analysis
geneList <- nrDEG_subtype$logFC
names(geneList) <- rownames(nrDEG_subtype)
geneList <- sort(geneList, decreasing = TRUE)

m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes <- fgsea(fgsea_sets, stats = geneList, nperm = 1000)
save(fgseaRes, file = "./results/step2_bulk/Bulk_Subtype_GSEA.Rdata")

# GSEA Visualization
topPathways <- fgseaRes[head(order(pval), n = 10), pathway]
plotGseaTable(fgsea_sets[topPathways], geneList, fgseaRes, 
              gseaParam = 0.5)

cat(">>> Step 2 Bulk Analysis Completed.\n")