library(Seurat)
library(CARD)
library(ggplot2)
library(patchwork)

load("D:/A-BS/ST_GBM/data/GBM_sc_V4.Rdata")
Idents(sc.v4)<-sc.v4$celltype

# 修改特定细胞类型
sc.v4$ct <- as.character(sc.v4$celltype)
sc.v4$ct[sc.v4$ct %in% c("Treg_T","Naive_T","CD8_T_EX","CD8_T_EM")] <- "T cell"
sc.v4$ct[sc.v4$ct %in% c("Neut","Mac","DC")] <- "Myeloid cell"

# 将修改后的列设置为idents
Idents(sc.v4) <- sc.v4$ct
table(sc.v4$ct)

wj<-list.files("D:/A-BS/ST_GBM_new/data/")

i=19
wj[i]
st1<-readRDS(paste0("D:/A-BS/ST_GBM_new/data/",wj[i],collapse = ""))
id<-gsub(".*?(\\d+_\\w+).*", "\\1", wj[i])
id
#手动改id
spatial_location<-st1@images$zh916t1@coordinates
spatial_location<-spatial_location[,4:5]
colnames(spatial_location)<-c("x","y")

CARD_obj = createCARDObject(
  sc_count = sc.v4@assays$RNA@counts,
  sc_meta = sc.v4@meta.data,
  spatial_count = st1@assays$Spatial@counts,
  spatial_location = spatial_location,
  ct.varname = "ct",
  ct.select = unique(sc.v4@meta.data$ct),
  sample.varname = "sample",
  minCountGene = 100,
  minCountSpot = 5) 

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

#Refined Grid
CARD_all <- CARD.imputation(CARD_obj,NumGrids = 3000,ineibor = 10,exclude = NULL)

location_imputation <- cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_all@refined_prop),split="x"),"[",1)),
                                        y=as.numeric(sapply(strsplit(rownames(CARD_all@refined_prop),split="x"),"[",2)))
rownames(location_imputation) <- rownames(CARD_all@refined_prop)

save(CARD_all,location_imputation,
     file = paste0("D:/A-BS/ST_GBM_new/CARD/SpatialData",id,"_CARD.Rdata",collapse=""))
############
cell_colors <- c(
  "T cell"="#EC8D63", 
  'B cell'="#3C5488",
  'Glia and neuronal cell'="#9370DB",
  'Oligodendrocyte'="#F08080",
  "Myeloid cell"="#00A087",
  'Pericyte'="#1E90FF",
  'Endothelial cell'="#FFA500"
)
p1 <- CARD.visualize.pie(proportion = CARD_all@Proportion_CARD,
                         spatial_location = CARD_all@spatial_location, colors = cell_colors)
pdf(file = paste0("D:/A-BS/ST_GBM_new/CARD/picture/",id,"_proportion.pdf",collapse=""),width = 7,height = 8)
plot(p1)
dev.off()

############
ct.visualize = unique(as.vector(sc.v4$ct))
p5 <- CARD.visualize.prop(
  proportion = CARD_all@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("lightblue","lightyellow","red"),    
  NumCols = 4)                                  

pdf(file = paste0("D:/A-BS/ST_GBM_new/CARD/picture/",id,"_celltype.pdf",collapse=""),width = 7,height = 8)
plot(p5)
dev.off()
##############
###增强分辨率下的基因表达
TF<-c("SOX11","SOX4","CEBPD","EGR1","RUNX3","ETS1","ELF1")
GENE<-c("AEBP1","IL13RA2","DCC","PRPS1","HDAC5","ASF1A","OPHN1","PDCD1","KLRB1","CLEC2D","LAG3","HAVCR2")
g<-c(TF,GENE)
g1<-g[g %in% rownames(CARD_all@refined_expression)]
p6 <- CARD.visualize.gene(
  spatial_expression = CARD_all@refined_expression,
  spatial_location = location_imputation,
  gene.visualize = g1,
  colors = NULL,
  NumCols = 6)
pdf(file = paste0("D:/A-BS/ST_GBM_new/CARD/picture/",id,"_gene.pdf",collapse=""),width = 7,height = 8)
plot(p6)
dev.off()
################
refined_prop <- CARD_all@refined_prop
# 根据GBM肿瘤微环境特征定义功能区域
tumor_regions <- list(
  "Tumor_region" = c("Glia and neuronal cell","Oligodendrocyte"),  # 肿瘤相关细胞
  "Stromal_region" = c("Pericyte","Endothelial cell"),   # 间质细胞
  "Immune_region" = c("T cell","B cell","Myeloid cell")  # 免疫细胞
)

# 计算每个网格点的区域组成
region_proportion <- sapply(names(tumor_regions), function(region) {
  cell_types <- tumor_regions[[region]]
  if(length(cell_types) == 1) {
    refined_prop[, cell_types]
  } else {
    rowSums(refined_prop[, cell_types])
  }
})

# 添加空间坐标信息
result_df <- cbind(location_imputation, region_proportion)

# 确定主导区域
result_df$Dominant_Region <- names(tumor_regions)[apply(region_proportion, 1, which.max)]

# 主导区域离散分布图
p1 <- ggplot(result_df, aes(x, y, color=Dominant_Region)) +
  geom_point(size=1.5, alpha=0.8) +
  scale_color_brewer(palette="Set1") +
  ggtitle("Dominant Functional Regions") +
  theme_bw() +
  coord_fixed()

# 特定区域连续比例图（示例：肿瘤核心）
p2 <- ggplot(result_df, aes(x, y, color=Tumor_region)) +
  geom_point(size=1.5) +
  scale_color_gradientn(colors=c("lightgrey", "yellow", "red"),
                        limits=c(0,1)) +
  ggtitle("Tumor Core Proportion") +
  theme_bw() +
  coord_fixed()

# 组合可视化
#p1 + p2
pdf(file = paste0("D:/A-BS/ST_GBM_new/CARD/picture/",id,"_region.pdf",collapse=""),width = 7,height = 8)
plot(p1 + p2)
dev.off()
