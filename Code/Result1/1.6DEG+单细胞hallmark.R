library(Seurat)
load("D:/A-BS/inferCNV/result/inferCNV_sce_NTI_new.Rdata")
Idents(scRNA)<-scRNA$type
#

#免疫细胞和胶质细胞和神经元 全部认定为正常细胞
cells.use <- WhichCells(scRNA, idents = 'Immune')
scRNA1 <- SetIdent(scRNA, cells = cells.use, value = 'Normal')
scRNA1$type<-Idents(scRNA1)
#library(devtools)
#devtools::install_github('immunogenomics/presto')
library(presto)
pbmc.genes <- wilcoxauc(scRNA1, 'type')

# 我们拥有每个cluster的所有基因
dplyr::count(pbmc.genes, group)

library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(ggrepel)

m_df<- msigdbr(species = "Homo sapiens", category = "H")

new_genes<- pbmc.genes %>%
  dplyr::filter(group == "Tumor") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

library(tibble)
ranks<- deframe(new_genes)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

save(fgseaRes,file="D:/A-BS/GSEA/scRNA_TumorVSN+I_hallmark.Rdata")

#############
#火山图
#自定义阈值
log2FC = 0.1
padj = 0.05 
new_genes<- pbmc.genes %>%
  dplyr::filter(group == "Tumor") %>%
  dplyr::filter(logFC!=0)%>%
  dplyr::filter(avgExpr!=0)
#####
data<-new_genes
index=data$padj<0.05 & abs(data$logFC)>log2FC
data$label<-0
data$label[index & data$logFC >0]=1
data$label[index & data$logFC <0]=-1
data$label<-factor(data$label,levels = c(1,0,-1),labels=c('up','ns','down'))

ggplot(data = data,aes(x=-log10(pval),y=logFC))+
  geom_point(aes(size = -log10(padj),color=-log10(padj)))+
  scale_color_gradientn(values = seq(0,1,0.05),
                        colors = c("#4E99C7","#A7CFE4","#D1E5F0","#FCD5BF","#DD6E57","#B2182B"))+
  scale_size_continuous(range = c(0.8,4))+
  labs(y="logFC",x="-log10(pval)")+
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
  annotate("text",x=-10,y=-1,label="Down-regulated",color="#1f78b4",size=5,lineheight=0.8,vjust=0)+
  annotate("text",x=-10,y=1,label="Up-regulated",color="#e31a1c",size=5,lineheight=0.8,vjust=0)+
  scale_x_continuous(limits = c(0, 300))+
  scale_y_continuous(limits = c(-1, 1))










