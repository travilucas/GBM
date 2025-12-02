library(xlsx)
#top 50 marker
im<-read.xlsx("D:/A-BS/ssGSEA/immune_markers_over.xlsx",sheetIndex = 1)

im_list<-as.list(im)
im_list[["B.cell"]]<-im_list[["B.cell"]][1:27]

library(GSVA)
load("D:/A-BS/UCSC-GBM/GBM_Tumor_merge_Combatover.Rdata")

gsvaP<-ssgseaParam(exprData = as.matrix(data),geneSets =im_list,minSize = 1,maxSize = Inf,alpha = 0.25,normalize = T)
gsva_re<-gsva(gsvaP)

ssGSEA_matrix<-gsva_re
save(ssGSEA_matrix,file = "D:/A-BS/ConsensusCluster/TCGA_tumor_ssGSEA_immune.Rdata")

######################
library(tibble)
library(magrittr)
library(reshape2)
library(ggplot2)
load( "D:/A-BS/ConsensusCluster/Tumor_ssGSEA_matrix.Rdata")

setwd("D:/A-BS/ConsensusCluster/")
title="D:/A-BS/ConsensusCluster/"
library(ConsensusClusterPlus)

results = ConsensusClusterPlus(as.matrix(ssGSEA_matrix), maxK = 9, reps = 1000, pItem = 0.8,
                               pFeature = 0.8, title = title, clusterAlg = "hc", 
                               distance = "pearson", seed = 1262118388.71279, plot = "pdf")

save(ssGSEA_matrix2,file = "./ssGSEA_matrix2.Rdata")

annCol <- data.frame(results = paste0("Cluster",
                                      results[[8]][['consensusClass']]),
                     row.names = colnames(ssGSEA_matrix))

head(annCol)

adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

a<-adjustdata (annCol)

sample<-a[which(a$results==c("Cluster2","Cluster3")),]

ssGSEA_matrix2<-ssGSEA_matrix[,which(colnames(ssGSEA_matrix)%in% sample$`rownames(data)`)]

results = ConsensusClusterPlus(as.matrix(ssGSEA_matrix2), maxK = 9, reps = 1000, pItem = 0.8, 
                               pFeature = 0.8, title = title, clusterAlg = "hc", 
                               distance = "pearson", seed = 1262118388.71279, plot = "pdf")

save(results2,file = "./ssGSEA_matrix2_results.Rdata")

load("./ssGSEA_matrix2_results.Rdata")
# 绘制热图
library(RColorBrewer)
library(pheatmap)
annCol <- data.frame(results = paste0("Cluster",
                                      results[[2]][['consensusClass']]),
                     row.names = colnames(ssGSEA_matrix2))
mycol <- brewer.pal(7, "Set3")

annColors <- list(results = c("Cluster1" = mycol[1],
                              "Cluster2" = mycol[2]
))

heatdata <- results[[2]][["consensusMatrix"]]
dimnames(heatdata) <- list(colnames(ssGSEA_matrix2), colnames(ssGSEA_matrix2))
heatdata[1:3,1:3]

pheatmap(mat = heatdata,
         color = colorRampPalette((c("white","#6892c5")))(50),
         border_color = NA,
         annotation_col = annCol,
         annotation_colors = annColors,
         show_colnames = F,
         show_rownames = F,
)

############



load("./ssGSEA_matrix2.Rdata")
load("./ssGSEA_matrix2_results.Rdata")
annCol5 <- data.frame(results = paste0("Cluster",
                                      results[[5]][['consensusClass']]),
                     row.names = colnames(ssGSEA_matrix2))

a<-adjustdata (annCol5)

write.xlsx(a,file = "tumorSample186_clus2.xlsx")

select_data<-data[,which(colnames(data) %in% rownames(annCol))]
a1<-adjustdata(select_data)
write.table(a1,"tumorSample_select150.txt",quote = F,sep = "\t",row.names = F)

# 查看聚类结果
print(consensus_result)

# 可视化聚类结果
# plot函数可以生成多个图，以下仅展示部分
pdf("consensus_cluster_plot.pdf")
plot(consensus_result, main="Consensus Clustering Plot", sub="Sample Clustering")
dev.off()

pdf("D:/A-BS/16免疫细胞_箱线图.pdf",width = 12,height = 8)

ggplot(tmp1,aes(Celltype,Score,fill = Celltype)) + 
  geom_boxplot(outlier.shape = 15,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Score") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(16))

dev.off()
############################################
