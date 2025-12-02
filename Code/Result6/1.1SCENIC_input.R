#download packages and reference txt
#######################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
BiocManager::install(c("AUCell", "RcisTarget"),ask = F,update = F) 
BiocManager::install(c("GENIE3"),ask = F,update = F)  # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"),ask = F,update = F)
BiocManager::install("rbokeh")
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"),ask = F,update = F)
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"),ask = F,update = F)
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools)
devtools::install_github("aertslab/SCENIC") 

#devtools::install_github("aertslab/SCENIC", ref="v1.1.0")

packageVersion("SCENIC")

#  https://resources.aertslab.org/cistarget/
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
#             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
#             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs
dbFiles
basename(dbFiles[1])
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  #  (1041.7 MB)
  # 
}

##1, For human:
dbFiles<-c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather",
"https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather")
###############################
library(SCENIC)
## Load data
loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
library(SCopeLoomR)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)
dim(exprMat)
exprMat[1:4,1:4]
head(cellInfo)
table(cellInfo$CellType)
### Initialize settings
library(SCENIC)
# 保证 cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)
### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

#######################################
rm(list = ls()) 
library(Seurat) 
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
load("D:/A-BS/CellChat/182109_sce_chat_umap_over.Rdata")
scRNA<-subset(x=sce,idents=c("Glia and neuronal cell(Tumor-associated)","Oligodendrocyte(Tumor-associated)","CD8_T_EX","CD8_T_EM"))
exprMat  <-  as.matrix(scRNA@assays$RNA@data)

# 重命名列
metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident,
                                       nUMI = nCount_RNA, 
                                       nGene = nFeature_RNA)

setwd("D:/A-BS/SCENIC")
save(exprMat,file="4Cell_expMatrix.Rdata")
####################
setwd("D:/A-BS/SCENIC")
load("4Cell_expMatrix.Rdata")
load("4Cell_InfoMatrix.Rdata")
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  scRNA@meta.data[,c(18,3,2)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
save(cellInfo,file = "4Cell_InfoMatrix.Rdata")
### Initialize settings
library(SCENIC)
# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
mydbDIR <- "./cisTarget_databases"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather","hg19-tss-centered-5kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "5kb")

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

scenicOptions <- initializeScenic(org="hgnc",dbs = mydbs,
                                  dbDir=mydbDIR, nCores=1) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
################
setwd("D:/A-BS/SCENIC")
scenicOptions=readRDS(file="int/scenicOptions.Rds")
load("4Cell_expMatrix.Rdata")
load("4Cell_InfoMatrix.Rdata")
#scenicOptions <- initializeScenic(org="hgnc", 
#                                  dbDir="cisTarget_databases", nCores=1) 
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["5kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

rm(list = ls()) 
library(Seurat) 
library(SCENIC)
library(doParallel)

scenicOptions=readRDS(file="int/scenicOptions.Rds")

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))
