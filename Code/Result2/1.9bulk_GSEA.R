library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
################################
load("D:/A-BS/UCSC-GBM/Tumor_Cluster2 result/DEGoutput_up-down.Rdata")
TumorG<-DEG

geneList <- as.matrix(TumorG$logFC)
names(geneList) <- toupper(rownames(TumorG))
geneList <- sort(geneList,decreasing = T)

m_df<- msigdbr(species = "Homo sapiens", category = "H")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes<- fgsea(fgsea_sets, stats = geneList, nperm = 1000)

ggplot(fgseaResTidy[1:9,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色

save(fgseaRes,file="D:/A-BS/GSEA/bulk_Tumor_hallmark.Rdata")
################

load("D:/A-BS/GSEA/bulk_Tumor_hallmark.Rdata")
######################
library(GSEABase)
library(limma) 
library(clusterProfiler)
library(enrichplot)
library(GseaVis)
library(RColorBrewer)
library(tidyverse)
setwd("D:/A-BS/GSEA/")
geneset <- read.gmt("./h.all.v2024.1.Hs.symbols.gmt")

TumorG<-read.xlsx("D:/A-BS/UCSC-GBM/Tumor_clus2 result/diff.Tumor_clus2_DEG.xlsx",sheetIndex = 1)
rownames(TumorG)<-TumorG[,1]
TumorG<-TumorG[,-1]

TumorG<-TumorG[which(TumorG$P.Value<0.00001),]
TumorG<-TumorG[which(TumorG$adj.P.Val<0.00001),]

geneList <- as.matrix(TumorG$logFC)
names(geneList) <- toupper(rownames(TumorG))
geneList <- sort(geneList,decreasing = T)

gsea_results <- GSEA(
  geneList = geneList,
  TERM2GENE = geneset,
  pvalueCutoff = 0.05
  
)
save(gsea_results,file = "./TCGA_tumor_Hallmark_gsea_result_over.Rdata")

load("./scNT_Hallmark_gsea_result.Rdata")

list<-gsea_results@result[["ID"]]
for (i in 1:length(list)) {
  p <- gseaplot2(x=gsea_results,geneSetID=list[i]) 
  d <- paste("./TCGA_tumor_hallmark_result/",list[i],".pdf",sep="")
  pdf(file=d,family = "Times",width=10,height = 6)
  print(p)
  dev.off()
}

p <- gseaplot2(gsea_results,
               gsea_results@result[["ID"]][1:5],
               pvalue_table = TRUE)
#gseaplot2(gsea_results, gsea_results@result[["ID"]][1:17], pvalue_table = TRUE)

write.table(list,file = "TCGA_tumor_clus2_hallmarkID.txt",quote = F,sep = "\t",row.names = F)
#####################
genesetID<-c('HALLMARK_GLYCOLYSIS',
             'HALLMARK_ADIPOGENESIS','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','HALLMARK_HYPOXIA'
             )

mygene<-c('S100A9','AQP9','CCL7','CXCL6','CREB5','PTPRS','ZNF311')

color<-c("#A6CEE3","#1F7BB4","#B2DF8A","#33A92C")

gseaNb(object=gsea_results,
       termWidth = 20,
       legend.position = c(0.85,0.8),
       geneSetID = genesetID,
       curveCol = color,
       subPlot = 2,
       addPval = T,
       addGene = gene,
       htCol = c("#78D3AC","#74C1F0"),
       pvalX = 0.71,pvalY = 1.0)
