library(Seurat)
library(tidyverse)
library(patchwork)

library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
load("D:/A-BS/GSE182109/result/182109_Tcell_sce_annotation_over.Rdata")

#expression_matrix <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
cell_metadata$cell<-row.names(cell_metadata)
info<-cell_metadata[,c(19,18)]

Exp <-  as.data.frame(sce@assays$RNA@data)
Exp<-as.data.frame(t(Exp))

T_col<-c("#B699C6","#3778AD","#4EA74A","#202D70")
gene<-c('KLRB1','CLEC2D','PDCD1','LAG3','HAVCR2')

gene <- as.vector(gene)
Exp_plot <- Exp[,colnames(Exp) %in% gene]#提取需要作图得基因表达信息

#加载样本信息
Exp_plot$sam=info$celltype
Exp_plot$sam <- factor(Exp_plot$sam,levels=c("Naive_T","CD8_T_EM","CD8_T_EX","Treg_T"))

plist2<-list()
for (i in 1:length(gene)){
  bar_tmp<-Exp_plot[,c(gene[i],"sam")]
  colnames(bar_tmp)<-c("Expression","sam")
  
  my_comparisons2 <- list(c("Naive_T", "CD8_T_EX"))
  my_comparisons4 <- list(c("CD8_T_EM", "CD8_T_EX"))
  my_comparisons6 <- list(c("CD8_T_EX", "Treg_T"))
  pb1<-ggboxplot(bar_tmp,
                 x="sam",
                 y="Expression",
                 color="sam",
                 fill=NULL,
                 add = "jitter",
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.01,
                 font.label = list(size=30), 
                 palette = T_col)+theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")#
  pb1<-pb1+stat_compare_means(method="t.test",hide.ns = F,
                              comparisons =c(
                                             my_comparisons2,
                                             my_comparisons4,
                                             my_comparisons6),
                              label="p.signif")
  plist2[[i]]<-pb1 
}
plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
          plist2[[4]],plist2[[5]])