#pre process data
###
library(xlsx)
K2Result<-read.xlsx(file = "D:/A-BS/ConsensusCluster/tumorSample186_clus2.xlsx",sheetIndex = 1)
tumor<-tumor_combat[,which(colnames(tumor_combat) %in% K2Result$rownames.data.)]
tumor<-as.data.frame(tumor)

adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

a<-adjustdata (t(tumor))

a<-as.data.frame(a)

a<-a[match(K2Result$rownames.data.,a$V1),]
identical(a$V1,K2Result$rownames.data.)

exp<-as.data.frame(t(a[,-1]))
identical(colnames(exp),K2Result$rownames.data.)

a1<-adjustdata (exp)
#################
setwd("D:/A-BS/CIBERSORT/")
write.table(a1,file = "mergeTumorExp_clus2.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#######################estimate
es_marker<-read.xlsx("D:/A-BS/inferCNV/estimate.xlsx",sheetIndex  = 1)
immune_markers<-t(es_marker[3,3:143])
immune_markers<-as.vector(immune_markers)
immune_markers<-unlist(immune_markers)

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
setwd("D:/A-BS/CIBERSORT/")
library(estimate)
rawdata_file <- "mergeTumorExp_clus2.txt"
filterCommonGenes(input.f=rawdata_file,
                  output.f="my.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "my.gct",
              output.ds="estimate_score.gct",
              platform="affymetrix")

raw_data = read.table("estimate_score.gct",skip = 2,header = F,sep ='\t',stringsAsFactors=FALSE)
rownames(raw_data)=raw_data[,1]
raw_data= as.data.frame(t(raw_data[,3:ncol(raw_data)]))
row.names(raw_data)<-c(1:nrow(raw_data))
colnames(raw_data)[1] <- c("sample")
for(i in 2:5){
  raw_data[,i] <- as.numeric(as.character(raw_data[,i]))
}
library(xlsx)
library(dplyr)
phenotype_file <- read.xlsx(file = "D:/A-BS/ConsensusCluster/tumorSample186_clus2.xlsx",sheetIndex = 1)
phenotype_file <- phenotype_file[phenotype_file[,1] %in% raw_data[,"sample"],]
colnames(phenotype_file)<-c("sample","Group")
merge_matrix <- left_join(phenotype_file,raw_data,by="sample") #合并数据集
my_comparisons <- list( c("Cluster1", "Cluster2"))

max_labely <- max(merge_matrix[,"ImmuneScore"])
min_labely <- min(merge_matrix[,"ImmuneScore"])
my_labely <- c(max_labely, max_labely, max_labely, max_labely)
library(ggpubr)

ggboxplot(merge_matrix, x = "Group", y = "ImmuneScore",
          color = "Group",alpha = 0.8,size=.3,
          add = "jitter",palette = c("#E58579","#7e9bb7"),
          add.params = list(color = "Group",alpha = .5,size = 0.8)) +
  #stat_compare_means(comparisons = my_comparisons,method = "t.test",
  #                   label.y = my_labely,na.rm=TRUE)+
  stat_compare_means(aes(group = Group),
                     label = "p.signif",
                     method = "t.test",
                     label.y = my_labely,
                     hide.ns = T)+
  ggtitle("Immune Score") +theme_bw()+
  theme(legend.position="none",
        plot.title = element_text(size=15,face="bold",hjust = 0.5),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Immune Score") +
  xlab("Group") +
  ylim(min_labely, max_labely)
