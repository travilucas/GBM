#single cell DEG 可视化结果图
#load("D:/A-BS/ssGSEA/182109_sce_subtype_markers.Rdata")

library(Seurat)
library(tidyverse)
library(scRNAtoolVis)
library(RColorBrewer)

#devtools::install_github ("junjunlab/scRNAtoolVis")
library(scRNAtoolVis)

T_col<-c("#3778AD","#B699C6","#4EA74A","#202D70")
Myeloid_col<-c("#8F4C9A","#E57FB0","#A15528")
#"#E89118","#56A8D7",
cols<-c("#D6251F","#3778AD","#B699C6","#4EA74A","#202D70","#E57FB0","#A15528","#8F4C9A")

load("D:/A-BS/GSE182109/result/182109_sce_all_celltype_DEG.Rdata")
sce2.markers <- sce.markers%>%
  dplyr::filter(p_val_adj<0.05&pct.1>0.25) %>% 
  group_by(cluster)



type<-c("B cell","Naive_T","CD8_T_EM","CD8_T_EX","Treg_T","Mac","Neut","DC")
immune.markers<-sce2.markers[which(sce2.markers$cluster %in% type),]

list_genes1<-list(
  DC=c("HLA-DRA","HLA-DPB1","HLA-DPA1"),
  Mac=c("CSF1R","VSIG4","SIGLEC1"),
  Neut=c("IL1R2","ITGAM","FPR2")
  
)

list_genes2<-list(
  CD8_T_EX=c("CD3D","CD3E","CD8A","CD8B","LAG3","PDCD1","TOX"),
  Naive_T=c("IL7R","CCR7","LEF1"),
  T_EFF=c("CCR5","CXCR3"),
  Treg_T=c("TNFRSF4","BATF","TNFRSF18","FOXP3","IL2RA")
)
B_cell=c("MZB1","MS4A1","IGHG1","IGHG3","CD79A")


mygene<-unique(c(unlist(list_genes1),unlist(list_genes2),B_cell))

#jjVolcano(diffData = immune.markers,
#          myMarkers = mygene)



jjVolcano(diffData=immune.markers,
          tile.col=cols,myMarkers =mygene, 
          size=3,legend.position=c(0.8,0.2),
          aesCol=c("grey","#E31A1C"),
          cluster.order=unique(immune.markers$cluster),
          back.col='white',
          polar=T
          )

#############
#tumor vs normal
library(Seurat)
load("D:/A-BS/inferCNV/result/inferCNV_sce_NTI_new.Rdata")
Idents(scRNA)<-scRNA$type
#
scRNA1<-subset(x=scRNA,idents=c("Tumor","Normal"))
scRNA1$type<-Idents(scRNA1)

#library(devtools)
#devtools::install_github('immunogenomics/presto')
library(presto)
pbmc.genes <- wilcoxauc(scRNA1, 'type')

#自定义阈值
log2FC = 0.1
padj = 0.05 
new_genes<- pbmc.genes %>%
  dplyr::filter(group == "Tumor") %>%
  dplyr::filter(logFC!=0)%>%
  dplyr::filter(avgExpr!=0)
#####
data<-new_genes

data$change = as.factor(ifelse(data$pval < 0.05 & abs(data$logFC) > log2FC,
                              ifelse(data$logFC > logFC_cutoff ,'up','down'),'ns')
)
library(ggrepel)
ggplot(data = data,aes(x=-log10(pval),y=logFC))+
  geom_point(aes(size = -log10(padj),color=-log10(padj)))+
  scale_color_gradientn(values = seq(0,1,0.05),
                        colors = c("#4E99C7","#A7CFE4","#D1E5F0","#FCD5BF","#DD6E57","#B2182B"))+
  scale_size_continuous(range = c(0.8,5))+
  labs(y="logFC",x="-log10(Pvalue)")+
  geom_hline(yintercept = c(-log2FC,log2FC),lty=4,lwd=0.6,alpha=0.8)+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.8,linetype = "solid"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 12,color = "black"),
        axis.title.x = element_text(size = 14,color = "black"),
        axis.title.y = element_text(size = 14,color = "black"),
        axis.line = element_blank()
        )+
  geom_text_repel(data = subset(data,abs(logFC)>=log2FC & padj<0.05),
                  aes(label=feature),col="black",alpha=0.8
                  )+
  annotate("text",x=100,y=-1,label="Down-regulated",color="#1f78b4",size=5,lineheight=0.8,vjust=0)+
  annotate("text",x=100,y=1,label="Up-regulated",color="#e31a1c",size=5,lineheight=0.8,vjust=0)+
  scale_x_continuous(limits = c(0, 300))+
  scale_y_continuous(limits = c(-1, 1))

###########method1
DEG<-new_genes
logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change = as.factor(ifelse(DEG$pval < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)

g = ggplot(data=DEG, aes(x=logFC, y=-log10(pval), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))## corresponding to the levels(res$change)
 

print(g)
############method2
library(ggrepel)
DEG$row<-rownames(DEG)
ggplot(DEG,aes(x=logFC, y=-log10(pval)))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(-1,1),lty=4,lwd=0.6,alpha=0.8)+
  geom_point(aes(size=-log10(pval),color=-log10(pval)))+
  scale_size_continuous(range = c(1,4))+
  scale_y_continuous(expand=expansion(add = c(0,0)))+
  scale_color_gradientn(values = seq(0,1,0.1),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  theme(panel.border = element_rect(fill=NA,color = "black",linewidth =1,linetype = "solid"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position.inside  = c(0.1,0.65),
        axis.text = element_text(size = 12,color = "black"),
        axis.title.x = element_text(size = 14,color = "black"),
        axis.title.y = element_text(size = 14,color = "black"),
        axis.line = element_blank()
  )+
  geom_text_repel(data = subset(DEG,abs(logFC)>=0.8 & pval<0.05),
                  aes(label = row,colour =-log10(pval)),alpha=0.8)
