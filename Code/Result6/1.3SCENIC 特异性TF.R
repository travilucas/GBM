library(SCENIC)

#### 2.加载SeuratData
library(SeuratData) #加载seurat数据集  
data("pbmc3k")  
seurat.data = pbmc3k
#seurat.data = readRDS("./pbmc3k.test.seurat.Rds")
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = F)

n.pcs = 30
seurat.data <- seurat.data %>% 
  RunUMAP(reduction = "pca", dims = 1:n.pcs, verbose = F) %>% 
  FindNeighbors(reduction = "pca", k.param = 10, dims = 1:n.pcs)

# 这里有自带的注释
seurat.data$seurat_annotations[is.na(seurat.data$seurat_annotations)] = "B"
Idents(seurat.data) <- "seurat_annotations"
DimPlot(seurat.data,reduction = "umap",label=T ) 


sub_regulonAUC <- regulonAUC[,match(rownames(cellInfo),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$seurat_annotations))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$seurat_annotations)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

#保存一下
save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')


#########################
setwd("D:/A-BS/SCENIC/SCENIC_result")
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

load("D:/A-BS/SCENIC/4Cell_InfoMatrix.Rdata")
sub_regulonAUC <- regulonAUC[,match(rownames(cellInfo),colnames(regulonAUC))]

### 4.1. TF活性均值
# 看看不同单细胞亚群的转录因子活性平均值
# Split the cells by cluster:
selectedResolution <- "CellType" # select resolution
cellsPerGroup <- split(rownames(cellInfo), 
                       cellInfo[,selectedResolution])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)

regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
new<-regulonActivity_byGroup_Scaled[which(rownames(regulonActivity_byGroup_Scaled) %in% rssPlot[["rowOrder"]]),]
Heatmap(new,
  #regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

### 4.2. rss查看特异TF
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellInfo[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss,zThreshold = 0.1,thr=0.1)
rssPlot

plotly::ggplotly(rssPlot$plot)

###其他方式
rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)
