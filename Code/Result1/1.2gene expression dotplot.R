library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)
library(ggplot2)
######
load("D:/A-BS/GSE182109/result/182109_annotation_over.Rdata")

Idents(sce)<-"celltype"
#182109
list_genes<-list(
  B_cell=c("MZB1","MS4A1","IGHG1","IGHG3","CD79A"),
  Endo=c("VWF","ITM2A","CLEC14A","A2M","CLDN5"),
  Glia_neuron=c("CLU","C1orf61","SOX2","TUBA1A","PTPRZ1","FABP7"),
  Myeloid=c("C1QA","C1QB","C1QC","SPP1","FCER1G"),
  Oligo=c("PLLP","MBP","CNP","CLDN11","PLP1"),
  Pericyte=c("RGS5","PDGFRB","NOTCH3","DCN","CD248"),
  T_cell=c("CD2","CD3D","CD3E","CD3G","CD52")
)

dot <- DotPlot(sce, features = list_genes)$data
###调整marker的顺序
dot$id <- factor(dot$id,levels = unique(as.character(dot$id))[order(unique(as.character(dot$id)))])
####做注释文件
data.anno <- data.frame(
  features.plot = rev(unique(dot$features.plot)),
  label = rev(c(rep("B cell",5),
                rep("Glia and neuronal cancer cell",6),rep("Endothelial cell",5),
                rep("Oligodendrocyte cancer cell",5),rep("Myeloid cell",5),
                rep("Pericyte",5),
                rep("T cell",5)
                )
              )
)

df.plot <- plyr::join(dot,data.anno)
df.plot$features.plot <- factor(df.plot$features.plot,levels = rev(levels(df.plot$features.plot)))

p <- ggplot(df.plot,aes(x=features.plot,y =  id,size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("% detected", range = c(0,5)) + #调整绘图点的相对大小
  scale_color_gradient2(low = "darkgrey",mid = "white",high = "red",midpoint = 0,
                        name ="Average\nexpression" 
                        ) +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id))+coord_flip()+#一开始设定的x轴和y轴错误，现在调转
  facet_grid(~id, scales="free_x",space = "fixed")+#theme_classic() +
  theme(
    axis.text.x = element_text(size=8, angle=0, hjust=0.5, color="black",face="bold"),#x轴标签样式
    #axis.text.y = element_blank(),####去掉y轴刻度
    axis.title.x = element_text(size=18,colour = 'black',vjust = -0.8,hjust = 0.5),#坐标轴标题
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
####scale_fill_npg的前五个配色
cols=c("#3C5488","#FFA500","#9370DB","#00A087","#F08080","#1E90FF","#EC8D63","#87CEEB","#DD9E82","#7CFC00","#7B68EE","#800080")
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
pdf(file = "D:/A-BS/inferCNV_TumorCell_dotplot.pdf",width = 7,height = 8)
plot(g)
dev.off()
