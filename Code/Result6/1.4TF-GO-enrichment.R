library(clusterProfiler)
library("org.Hs.eg.db")
library(xlsx)
setwd("D:/A-BS/SCENIC/TF-cor-network-GO/")
net<-read.xlsx("RUNX3.xlsx",sheetIndex = 1)

gene<-net$gene
#gene symbol id to entrez id
gene1 = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene2 = gene1$ENTREZID

#GO分析
enrich_go <- enrichGO(gene=gene2, 
                      OrgDb = 'org.Hs.eg.db', 
                      keyType = 'ENTREZID', 
                      ont = 'BP', #BP,CC,MF
                      pvalueCutoff = 0.05, #pvalue cutoff
                      qvalueCutoff = 0.05, #qvalue cutoff
                      readable = TRUE #whether mapping gene ID to gene Name
)

write.csv(summary(enrich_go),"RUNX3_GO_BPenrich.csv",row.names =F)

RUNX3<-summary(enrich_go)

#SOX4<-read.csv("SOX4_GO_BPenrich.csv",header = T)

term<-intersect(ELF1$Description,RUNX3$Description)

write.table(term,"ELF1-RUNX3_intersect_BPid.txt",quote = F,sep = "\t",row.names = F)
