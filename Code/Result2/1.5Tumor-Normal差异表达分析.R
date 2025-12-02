#Tumor Normal DEG
load("D:/A-BS/UCSC-GBM/GBM_Tumor_merge_over.Rdata")
plat<-c(rep("TCGA-LGGGBM",108),rep("TCGA-GBM",144),rep("CGGA325",71),rep("CGGA693",105))
library(sva)
tumor_combat<-ComBat(dat = Tumor_merge, batch = plat)
save(tumor_combat,file ="D:/A-BS/UCSC-GBM/GBM_Tumor_merge_Combatover.Rdata" )
###
load("D:/A-BS/UCSC-GBM/GBM_Normal_merge_over.Rdata")
plat<-c(rep("TCGA-GBM",5),rep("GTEX",105))
library(sva)
normal_combat<-ComBat(dat = Normal_merge, batch = plat)
#########
tumor_combat<-tumor_combat[which(rownames(tumor_combat) %in% rownames(normal_combat)),]
normal_combat<-normal_combat[which(rownames(normal_combat) %in% rownames(Tumor_merge)),]

tumor_combat<-tumor_combat[order(row.names(tumor_combat)), ]
normal_combat<- normal_combat[order(row.names(normal_combat)), ]

identical(rownames(tumor_combat),rownames(normal_combat))
GBM_merge<-cbind(tumor_combat,normal_combat)
exp<-as.data.frame(GBM_merge)
type<-c(rep("Tumor",428),rep("Normal",110))
################
library(limma) 
###############
design <- model.matrix(~0+factor(type))
colnames(design)=levels(factor(type))
rownames(design)=colnames(exp)
#???????ݡ??????ȽϾ???
contrast.matrix<-makeContrasts(paste0(c("Tumor","Normal"),collapse = "-"),levels = design)

fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##??һ??????Ҫ?????ҿ??????п???Ч??
fit2 <- eBayes(fit2)  
#default no trend !!! eBayes() with trend=TRUE
#??ͼ???? ??ִ?д???
results <- decideTests(fit2) 
vennDiagram(results) #????Τ??ͼ
#?õ?????????????
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#######################################
#??????????????????
write.table(nrDEG, "D:/A-BS/UCSC-GBM/NT_result/diff.ControlVSTumor.txt",sep = '\t',quote = F)
save(nrDEG,file = "D:/A-BS/UCSC-GBM/NT_result/nrDEGoutput.Rdata")
#
load(file = "D:/A-BS/UCSC-GBM/NT_result/nrDEGoutput.Rdata")
colnames(nrDEG)
plot(nrDEG$logFC,-log10(nrDEG$P.Value))
library(ggplot2)
DEG=nrDEG
logFC_cutoff<-1
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
this_tile
head(DEG)
library(ggplot2)
g = ggplot(data=DEG, 
           aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.6, size=1) +
  scale_colour_manual(values = c('#71C6D4','#ADB6B6FF','#71C6D4'))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),lty=4,lwd=0.6,alpha=0.8)
  
print(g)


g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)
save(DEG,file = "D:/A-BS/UCSC-GBM/NT_result/DEGoutput_up-down.Rdata")
