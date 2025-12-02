load("D:/A-BS/UCSC-GBM/GBM_Tumor_merge_Combatover.Rdata")
info<-read.xlsx("D:/A-BS/UCSC-GBM/GBM-Tumor_info_over.xlsx",sheetIndex = 2)
library(xlsx)

info$Sample<-gsub('[-]', '.', info$Sample)
colnames(tumor_combat)<-gsub('[-]', '.',colnames(tumor_combat))
tumor_combat<-as.data.frame(tumor_combat)

exp<-tumor_combat[,which(colnames(tumor_combat) %in% info$Sample)]
info_new<-info[match(colnames(exp),info$Sample),]

identical(colnames(exp),info_new$Sample)
type<-info_new$Cluster
#############
library(limma)
###############
design <- model.matrix(~0+factor(type))
colnames(design)=levels(factor(type))
rownames(design)=colnames(exp)
#???????ݡ??????ȽϾ???
contrast.matrix<-makeContrasts(paste0(c("Cluster1","Cluster2"),collapse = "-"),levels = design)#cold-hot
fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##??һ??????Ҫ?????ҿ??????п???Ч??
fit2 <- eBayes(fit2)
#?õ?????????????
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
save(nrDEG,file = "D:/A-BS/UCSC-GBM/Tumor_Cluster2 result/nrDEGoutput_coldVShot.Rdata")

load("D:/A-BS/UCSC-GBM/Tumor_Cluster2 result/nrDEGoutput_coldVShot.Rdata")
DEG=nrDEG
logFC_cutoff <-0.5
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
this_tile
###################
load("D:/A-BS/UCSC-GBM/Tumor_Cluster2 result/DEGoutput_up-down.Rdata")
library(ggrepel)
DEG$row<-rownames(DEG)
ggplot(DEG,aes(x=logFC, y=-log10(P.Value)))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(-1,1),lty=4,lwd=0.6,alpha=0.8)+
  geom_point(aes(size=-log10(P.Value),color=-log10(P.Value)))+
  scale_size_continuous(range = c(1,4))+
  scale_y_continuous(expand=expansion(add = c(0,0)))+
  scale_color_gradientn(values = seq(0,1,0.1),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  theme(panel.border = element_rect(fill=NA,color = "black",size=1,linetype = "solid"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.1,0.65),
        axis.text = element_text(size = 12,color = "black"),
        axis.title.x = element_text(size = 14,color = "black"),
        axis.title.y = element_text(size = 14,color = "black"),
        axis.line = element_blank()
  )+
  geom_text_repel(data = subset(DEG,abs(logFC)>=2 & P.Value<0.05),
                  aes(label = row,colour =-log10(P.Value)),alpha=0.8)


save(DEG,file = "D:/A-BS/UCSC-GBM/Tumor_Cluster2 result/DEGoutput_up-downcoldVShot.Rdata")
