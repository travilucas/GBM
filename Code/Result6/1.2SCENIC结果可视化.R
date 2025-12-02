setwd("D:/A-BS/SCENIC/SCENIC_result")
library(Seurat) 
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

load("D:/A-BS/SCENIC/4Cell_expMatrix.Rdata")
load("D:/A-BS/SCENIC/4Cell_InfoMatrix.Rdata")

logMat <- log2(exprMat+1) # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) # default t-SNE

aucellApp <- AUCell_createViewerApp(auc=cells_AUC, thresholds=selectedThresholds, 
                                    tSNE=cellsTsne, exprMat=exprMatrix, cellInfo=cellInfo)

savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# Show TF expression:
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("CREM", "IRF1", "JUNB","FOSB")],], plots="Expression")

# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

library(KernSmooth)

library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

#par(bg = "black")
par(mfrow=c(1,2))
regulonNames <- c( "Dlx5","Sox10")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(0, 10, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-20,-10, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
regulonNames <- list(red=c("Sox10", "Sox8"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
text(5, 15, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(5, 15-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(5, 15-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)

regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Dlx5", "Irf1")]
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
###############
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")

tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset)
######################
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, 
                   #fontsize_row=3,
                   cluster_cols =F,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
                   cluster_rows =T,#是否对行聚类
                   treeheight_row =150,
                   color=colorRampPalette(c("#6892c5","white","#9F0000"))(100),
                   breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10,
                   border_color=NA)
#############################
# filename="regulonActivity_byCellType.pdf", width=10, height=20)
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)

#######################################################
##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC

load("D:/A-BS/GSE182109/result/182109_Tcell_sce_annotation_over.Rdata")
#pbmc <-readRDS("/home/wucheng/jianshu/function/data/pbmc.rds")  #导入pbmc3k数据
#current.cluster.ids <- c(0:8)
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
pbmc@meta.data$celltype <- plyr::mapvalues(x = pbmc@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
head(pbmc@meta.data)


scRNAauc <- AddMetaData(object=sce, 
                        metadata=AUCmatrix,
                        col.name = "AUCell")
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

##利用Seurat可视化AUC
dir.create('scenic_seurat')
#FeaturePlot
p1 = FeaturePlot(scRNAauc, features='CEBPB_extended_2290g', label=T, reduction = 'tsne')
p2 = FeaturePlot(scRNAbin, features='CEBPB_extended_2290g', label=T, reduction = 'tsne')
p3 = DimPlot(scRNA, reduction = 'tsne', group.by = "celltype_Monaco", label=T)
plotc = p1|p2|p3
ggsave('scenic_seurat/CEBPB_extended_2290g.png', plotc, width=14 ,height=4)

#################################################
#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "CEBPB_extended_2290g", group.by="celltype_Monaco") + 
  theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features = "CEBPB_extended_2290g", pt.size = 0, group.by="celltype_Monaco") + 
  theme(legend.position='none')
plotc = p1 + p2
ggsave('scenic_seurat/Ridge-Vln_CEBPB_extended_2290g.png', plotc, width=10, height=8)

#
library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'CellType')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#挑选部分感兴趣的regulons

my.regulons<-c('JUNB_68g','SOX4_163g','YBX1_extended_639g','SOX2_extended_372g','EGR1_101g',
               'CEBPD_extend_74g','IRF1_18g','ETS1_514g','CREM_35g','ELF1_104g','STAT1_80g',
               'RUNX3_254g','STAT2_28g','IKZF1_115g'
               )

#example
my.regulons <- c('ETS1_2372g','ETV7_981g','IRF7_239g','XBP1_854g','ATF4_37g',
                 'KLF13_78g','ATF6_129g','CREB3L2_619g','TAGLN2_13g',
                 'STAT1_extended_1808g','CEBPB_extended_2290g','IRF5_extended_422g',
                 'SPI1_1606g','HMGA1_14g','SPIB_1866g','IRF8_348g','BCL11A_136g',
                 'EBF1_40g','MAF_45g','BATF_131g','FOXP3_55g','TBX21_388g',
                 'EOMES_extended_101g','TCF7_extended_31g','LEF1_extended_49g')
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
#使用regulon原始AUC值绘制热图
ann_colors<-list(
  CellType=c("Glia and neuronal cell(Tumor-associated)"="#EEA599","Oligodendrocyte(Tumor-associated)"="#D9BDD8","CD8_T_EX"="#93B237","CD8_T_EM"="#92B4C8")
  
)

p<-pheatmap(myAUCmatrix, show_colnames=F,
            cluster_cols =F,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
            cluster_rows =T,#是否对行聚类
            annotation_colors =ann_colors,
            annotation_col=celltype)

pdf(file = "D:/A-BS/scenic.pdf",width = 7,height = 8)
plot(p)
dev.off()

         #filename = 'scenic_seurat/myAUCmatrix_heatmap.png',
         #width = 6, height = 5)
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         filename = 'D:/A-BS/scenic.pdf',
         color = colorRampPalette(colors = c("white","black"))(100),
         width = 6, height = 5)

