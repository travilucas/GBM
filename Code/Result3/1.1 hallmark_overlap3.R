library(dplyr)
#scRNA
load("D:/A-BS/GSEA/scRNA_TumorVSN+I_hallmark.Rdata")
sc_hallmark<-fgseaRes
sc_hallmark <- sc_hallmark %>%
  as_tibble() %>%
  arrange(desc(NES))
#bulk tumor clus2
load("D:/A-BS/GSEA/bulk_Tumor_hallmark_coldVShot.Rdata")
bulk_hallmark<-fgseaRes
bulk_hallmark <- bulk_hallmark %>%
  as_tibble() %>%
  arrange(desc(NES))
#################
#method3
#先交后并
sc<-sc_hallmark[1:28,]
bulk<-bulk_hallmark[1:6,]
id<-intersect(sc$pathway,bulk$pathway)
write.table(id,"D:/A-BS/GSEA/一致性hallmark_id.txt",quote = F,sep = "\t")
###############
#sc level 富集到的gene
sc_new<-sc[sc$pathway %in% id,]
#bulk level富集到的基因
bulk_new<-bulk[bulk$pathway %in% id,]
gene1<-intersect(unlist(sc_new$leadingEdge),unlist(bulk_new$leadingEdge))
#######
#TCGA tumor&normal DEG
load("D:/A-BS/UCSC-GBM/NT_result/DEGoutput_up-down.Rdata")
UP<-DEG[which(DEG$change=="UP"),]
###########
gene2<-intersect(unlist(bulk$leadingEdge),rownames(UP))
gene3<-intersect(unlist(sc$leadingEdge),rownames(UP))

g1<-union(gene1,gene2)
g2<-union(g1,gene3)
save(g2,file = "D:/A-BS/mime1_new/汇报版/hallmark_gene349.Rdata")
####
#韦恩图
library(ggvenn)
x=list("Bulk enriched hallmark genes"=unlist(bulk$leadingEdge),
       "sc enriched hallmark genes"=unlist(sc$leadingEdge),
       "DEGs_UP"=rownames(UP))
ggvenn(x,
       c("Bulk enriched hallmark genes","sc enriched hallmark genes", "DEGs_UP"),
       show_percentage = F,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295"))
