library(dplyr)
library(Seurat)
library(patchwork)

load("D:/A-BS/GSE182109/result/182109_Tcell_sce_not_annotation.Rdata")
Idents(scRNA1)<-scRNA1$RNA_snn_res.0.9

cells_to_pick <- WhichCells(scRNA1, idents = c(0,2,3,5,6,8,9,13))

# 获取这些细胞的数据
sce<- subset(scRNA1, cells = cells_to_pick)

new.cluster.ids <-read.delim("clipboard",header=FALSE)
new.cluster.ids<-as.vector(new.cluster.ids)
new.cluster.ids<-unlist(new.cluster.ids)
names(new.cluster.ids) <- levels(sce) 
sce <- RenameIdents(sce, new.cluster.ids)       #修改Idents

#在metadata中，添加Celltype信息
sce$celltype <- Idents(sce)

list_genes<-list(
  CD8_T_EX=c("CD4","CD3D","CD3E","CD8A","CD8B","LAG3","PDCD1","TOX"),
  CD8_T_EM=c("CXCR4","GZMB","CD44"),
  Naive_T=c("IL7R","CCR7","LEF1"),
  Treg_T=c("TNFRSF4","BATF","TNFRSF18","FOXP3","IL2RA")
)
library(ggrepel)

dot <- DotPlot(sce, features = list_genes)$data
dot$id <- factor(dot$id,levels = unique(as.character(dot$id))[order(unique(as.character(dot$id)))])
#dot <- dot[order(dot$id),]
####做注释文件
data.anno <- data.frame(
  features.plot = rev(unique(dot$features.plot)),
  label = rev(c(rep("CD8_T_EX",8),rep("CD8_T_EM",3),
                rep("Naive_T",3),rep("Treg_T",5)
  )
  )
)

df.plot <- plyr::join(dot,data.anno)
df.plot$features.plot <- factor(df.plot$features.plot,levels = rev(levels(df.plot$features.plot)))

p <- ggplot(df.plot,aes(x=features.plot,y =  id,size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("% detected", range = c(0,8)) + #调整绘图点的相对大小
  scale_color_gradient2(low = "grey",mid = "white",high = "red",midpoint = 0,
                        name ="Average\nexpression" 
  ) +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id))+coord_flip()+#一开始设定的x轴和y轴错误，现在调转
  facet_grid(~id, scales="free_x",space = "free")+#theme_classic() +
  theme(
    axis.text.x = element_text(size=8, angle=0, hjust=0.5, color="black",face="bold"),#x轴标签样式
    #axis.text.y = element_blank(),####去掉y轴刻度
    axis.title.x = element_text(size=20,colour = 'black',vjust = -0.8,hjust = 0.5),#坐标轴标题
    axis.ticks.y = element_blank(),#坐标轴刻度
    axis.text.y.right = element_blank(),#坐标轴标签隐藏
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',linewidth = 0.2), #坐标轴轴线样式
    
    panel.spacing=unit(1, "mm"), #分面图图块间距
    strip.text.x = element_text(size=8, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),#分面标签样式
    strip.background.x = element_rect(linewidth = 1)
  )+scale_y_discrete(breaks=NULL)

p
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))

cols<-c("#DF9E9B","#D8E7CA","#99BADF","#99CDCE","#999ACD")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
pdf(file = "D:/A-BS/TCell_dotplot.pdf",width = 7,height = 8)
plot(g)
dev.off()

save(sce,file = "D:/A-BS/GSE182109/result/182109_Tcell_sce_annotation_over.Rdata")

load("D:/A-BS/GSE182109/result/182109_Tcell_sce_annotation_over.Rdata")
