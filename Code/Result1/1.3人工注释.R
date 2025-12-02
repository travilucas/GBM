library(Seurat)
library(tidyverse)
library(ggrepel)

#细胞注释
new.cluster.ids <-read.delim("clipboard",header=FALSE)
new.cluster.ids<-as.vector(new.cluster.ids)
new.cluster.ids<-unlist(new.cluster.ids)

names(new.cluster.ids) <- levels(scRNA_combined) 
scRNA_combined <- RenameIdents(scRNA_combined, new.cluster.ids)       #修改Idents

#在metadata中，添加Celltype信息
scRNA_combined$celltype <- Idents(scRNA_combined)

head(scRNA_combined@meta.data)
#降维图
Idents(scRNA)<-"celltype"
DimPlot(scRNA,                                #seurat对象(数据)
        reduction = 'umap',                  #降维类型，umap，tsne，pca
        group.by = 'orig.ident',               #分组名称
        pt.size = 1,                         #图中点的大小
        #split.by = 'seurat_clusters',       #拆分绘图的名称
        label = T                          #图中是否添加label信息
        
        )                           

######################################################
load("D:/A-BS/GSE182109/result/182109_annotation_over.Rdata")
Idents(sce)<-"celltype"
###############################################
umap = sce1@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = sce1@meta.data$celltype) # 注释后的label信息 ，改为cell_type
colnames(umap)<-c("umap_1","umap_2","cell_type")
#131928
umap = scRNA1@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = scRNA1@meta.data$celltype)

head(umap)
###################################
allcolour=c("#3C5488","#FFA500","#9370DB","#00A087","#F08080","#1E90FF","#EC8D63")

p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  scale_color_manual(values = allcolour)
p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))

p3 <- p2 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 

p4 <- p3 + 
  geom_segment(aes(x = min(umap$umap_1) , y = min(umap$umap_2) ,
                   xend = min(umap$umap_1) +3, yend = min(umap$umap_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$umap_1)  , y = min(umap$umap_2)  ,
                   xend = min(umap$umap_1) , yend = min(umap$umap_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$umap_1) +1.5, y = min(umap$umap_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$umap_1) -1, y = min(umap$umap_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p4

cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

p5<-p4 +
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none")

pdf(file = "D:/A-BS/umap_182109.pdf",width = 7,height = 8)
plot(p5)
dev.off()


###########################
###堆积柱状图
table(sce$orig.ident)
table(sce$sex)#查看各组细胞数
prop.table(table(Idents(sce)))
table(Idents(sce), sce$sample)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(sce), sce$sample), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1,levels=c("B cell","Endothelial cell",
                                              "Glia and neuronal cell","Myeloid cell",
                                              "Oligodendrocyte","Pericyte","T cell"))

allcolour=c("#3C5488","#FFA500","#9370DB","#00A087","#F08080","#1E90FF","#EC8D63","#87CEEB","#DD9E82","#7CFC00","#7B68EE","#800080")
##
library(ggplot2)
ggplot(Cellratio,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_bar(stat = "identity",width = 0.7,position = "fill")+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))





