library(dplyr)
library(Seurat)
library(patchwork)
library(irlba)

load("D:/A-BS/inferCNV/inferCNV/182109_annotation_no_myeloid.Rdata")
#####
infer_CNV_obj<-readRDS('D:/A-BS/inferCNV/result/run.final.infercnv_obj')
expr<-infer_CNV_obj@expr.data
expr[1:4,1:4]
data_cnv<-as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)

phe<-sce1@meta.data
dim(phe)
sce1=sce1[,colnames(sce1) %in% rownames(phe)]

identical(colnames(sce1), rownames(phe))
sce1$celltype = phe$celltype
table(sce1$celltype )
meta <- sce1@meta.data


tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`T cell`]
tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`B cell`]
tmp= cbind(tmp1,tmp2)
down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
oneCopy=up-down
oneCopy
a1= down- 2*oneCopy
a2= down- 1*oneCopy
down;up
a3= up +  1*oneCopy
a4= up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # #####
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  cnv_score_table_pts[1:4,1:4]
  str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 


head(cell_scores_CNV) 
score=cell_scores_CNV
head(score)
meta$totalCNV = score[match(colnames(sce.all.int),
                            rownames(score)),1] 


setwd("D:/A-BS/inferCNV/result/")
#load("./inferCNV_scoreTable_result.Rdata")
load("./inferCNV_score_result.Rdata")
#inferCNV
#load("./inferCNV_Cellscore_result.Rdata")

library(ggplot2)
cols=c("#3C5488","#FFA500","#9370DB","#F08080","#1E90FF","#EC8D63")

ch<-ggplot(meta, aes(x=celltype  , y=totalCNV, fill=celltype)) +
  geom_boxplot()
# Assign custom color
ch+scale_fill_manual(values=cols)

pic <- ggplot(data = meta, aes(x = celltype, y = totalCNV, color = celltype))+
  geom_boxplot(size = 0.8, width = 0.8, alpha = 0)+
  geom_jitter(position = position_jitter(0.4), alpha = 0.2, size=1.5)+
  stat_summary(fun = mean, geom = "line", aes(group = celltype), size = 1.5, color = "red") 
  #labs(title = "")
pic+theme_classic(values=cols)
pic+scale_fill_manual(values=cols)

library(ggpubr)
library(patchwork)#如果没有安装要先安装

cols=c("#9370DB","#F08080","#EC8D63","#3C5488","#1E90FF","#FFA500")
ggboxplot(meta, x = "celltype", y = "totalCNV",
bxp.errorbar=T,#显示误差条`
width = 0.6,#箱体的宽度`
color = "celltype", #分组`
palette=cols#使用杂志aaas的配色
)+geom_hline(yintercept = mean(meta$totalCNV))


#########################################umap
load("D:/A-BS/inferCNV/inferCNV/182109_annotation_no_myeloid.Rdata")
load("./inferCNV_score_result.Rdata")
##################
#select tumor cell
library(dplyr)
library(Seurat)
library(patchwork)
library(irlba)

mean_Value<-mean(meta$totalCNV)
tumor_meta<-meta[which(meta$celltype==c("Glia and neuronal cell","Oligodendrocyte")),]
tumor<-tumor_meta[which(tumor_meta$totalCNV >= mean_Value),]
normal<-tumor_meta[which(tumor_meta$totalCNV < mean_Value),]
immune_meta<-meta[which(meta$celltype==c("T cell","B cell")),]
tumor[which(tumor$celltype=="Glia and neuronal cell"),18]<-"Glia and neuronal cell(Tumor-associated)"
tumor[which(tumor$celltype=="Oligodendrocyte"),18]<-"Oligodendrocyte(Tumor-associated)"

tumor$type<-"Tumor"
normal$type<-"Normal"
immune_meta$type<-"Immune"

new_meta<-rbind(tumor,normal,immune_meta)

save(new_meta,file="./inferCNV_sce_Cluster_new.Rdata")
##############
setwd("D:/A-BS/inferCNV/result/")

load("./inferCNV_sce_Cluster_new.Rdata")
load("D:/A-BS/inferCNV/inferCNV/182109_annotation_no_myeloid.Rdata")

head(colnames(sce1))
sce_tumor<-subset(x=sce1,cells = rownames(new_meta))

Idents(sce_tumor)<-new_meta$celltype
sce_tumor$celltype<-Idents(sce_tumor)

##########umap 绘制
umap = sce_tumor@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = sce_tumor@meta.data$celltype) # 注释后的label信息 ，改为cell_type
colnames(umap)<-c("umap_1","umap_2","cell_type")
############
cols<-c("#E89118","#56A8D7","#9370DB","#F08080","#00A087","#3C5488")
p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  scale_color_manual(values = cols)
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


#################
load("D:/A-BS/GSE182109/result/GSE182109_infercnv_tumor_umap.Rdata")

scRNA<- RunPCA(sce_tumor,features = VariableFeatures(object = sce_tumor))
scRNA$celltype<-Idents(scRNA)
scRNA$type<-new_meta$type
#绘制肘部图可视化确定聚类个数
ElbowPlot(scRNA)

#准确确定线性降维的维数
#
pct<-scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu>90 & pct<5)[1]
co1
co2<-sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1),decreasing = T)[1]+1
co2
pcs<-min(co1,co2)
pcs
plot_df<-data.frame(pct=pct,cumu=cumu,rank=1:length(pct))

library(ggplot2)
ggplot(plot_df,aes(cumu,pct,label=rank,color=rank>pcs))+
  geom_text()+
  geom_vline(xintercept = 90,color='grey')+
  geom_hline(yintercept = min(pct[pct>5]),color="grey")+
  theme_bw()
#UMAP
scRNA<- RunUMAP(scRNA, dims = 1:21)
save(scRNA,file="./inferCNV_sce_NTI_new.Rdata")
load("./inferCNV_sce_NTI_new.Rdata")

FeaturePlot(scRNA, features = c("MIF", "PTN", "SPP1" ))

######
#umap
######
umap = scRNA@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = scRNA@meta.data$celltype) # 注释后的label信息 ，改为cell_type
colnames(umap)<-c("umap_1","umap_2","cell_type")

cols<-c("#DF9E9B","#99BADF","#D8E7CA","#99CDCE","#999ACD","#FFD0E9")

cols <- c("#20B2AA","#1E90FF","#9370DB","#F08080","#EC8D63","#3C5488","#FFA500","#9370DB","#F08080","#20B2AA","#98FB98","#1E90FF") %>% as.vector() 

library(ggplot2)
library(ggrepel)

p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +  
  geom_point(size =0.8 , alpha =1 )  +  scale_color_manual(values = cols)
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

pdf(file = "D:/A-BS/umap_182109_tumor.pdf",width = 7,height = 8)
plot(p5)
dev.off()

#head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 0 from clusters 0 and 3
Idents(scRNA)<-scRNA$type
NT_ALLmarkers <- FindMarkers(scRNA, ident.1 = "Tumor", ident.2 = "Normal")
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

a<-adjustdata (NT_ALLmarkers)
save(NT_ALLmarkers,file="D:/A-BS/inferCNV/recluster_cell/GSE182109_inferCNV_NTcell_DEG.Rdata")

write.table(a,file = "D:/A-BS/inferCNV/recluster_cell/GSE182109_inferCNV_NTcell_DEG.txt",sep = "\t",quote = FALSE,row.names = F)

