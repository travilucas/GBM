library(xlsx)
load("D:/A-BS/ConsensusCluster/ssGSEA_matrix2.Rdata")
#clus2
K2Result<-read.xlsx(file = "D:/A-BS/ConsensusCluster/tumorSample186_clus2.xlsx",sheetIndex = 1)

info<-read.xlsx("D:/A-BS/UCSC-GBM/GBM-Tumor_info.xlsx",sheetIndex = 1)

info$Sample<-gsub('[-]', '.', info$Sample)
colnames(ssGSEA_matrix2)<-gsub('[-]', '.',colnames(ssGSEA_matrix2))

adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

a<-adjustdata (t(ssGSEA_matrix2))
a<-as.data.frame(a)

info_new<-info[match(K2Result$rownames.data.,info$Sample),]
info_new<-na.omit(info_new)
K2Result<-K2Result[which(K2Result$rownames.data. %in% info_new$Sample),]
identical(info_new$Sample,K2Result$rownames.data.)
info_new$Cluster<-K2Result$results
write.xlsx(info_new,"D:/A-BS/UCSC-GBM/GBM-Tumor_info_over.xlsx")
################ 
info<-read.xlsx("D:/A-BS/UCSC-GBM/GBM-Tumor_info_over.xlsx",sheetIndex = 2)
info$status <- ifelse((info$OS==0)==TRUE,"Alive","Dead")

a<-a[match(info$Sample,a$V1),]
identical(a$V1,info$Sample)
###
matrix<-a[,-1]
annotation_col <- data.frame(
  Cluster = as.vector(info$Cluster),
  Gender=as.vector(info$Gender),
  Status=as.vector(info$status),
  Age=as.vector(info$age1)
  
)
rownames(annotation_col)<-rownames(matrix)
ann_colors<-list(
  Cluster=c(Cluster1="#EEA599",Cluster2="#92B4C8"),
  #,Cluster3="#D9BDD8",Cluster4="#9180AC",Cluster5="#9DD0C7",Cluster6="#A9B6C1"),
  Gender=c("Male"="#D9BDD8","Female"="#458DA3","Female "="#458DA3"),#"Not Reported"="lightgrey"),
  Status=c("Alive"="#93B237","Dead"="#F6C957"),#"Not Reported"="lightgrey"),
  Age=c("<60"="#FFB77F",">=60"="#8086bb")
)

m1<-as.data.frame(lapply(matrix,as.numeric))
rownames(m1)<-rownames(matrix)
########
m1<-t(m1)
m1<-as.data.frame(m1)
###########
#1#数据缩放
result1<-scale(m1,center = F,scale = T) 

library(pheatmap)

gap<-c(135) #clus2

bk <- c(seq(-2,-1.9,by=0.01),
        seq(-1.8,-0.50,by=0.01),
        seq(-0.45,-0.3,by=0.01),
        seq(-0.2,-0.1,by=0.01),
        seq(0,2,by=0.01)
        #seq(0.3,2,by=0.01)
)

pheatmap(as.matrix(result1),
         scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         cluster_cols =F,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         cluster_rows =T,#是否对行聚类
         show_colnames = F,
         show_rownames = T,
         treeheight_row =80, #调节行 聚类树的长度
         clustering_method = "ward.D", #聚类方法
         color = c(
           colorRampPalette(colors=c("#6892c5","lightblue"))(length(bk)/6),
           colorRampPalette(colors=c("lightblue","white"))(length(bk)/4),
           colorRampPalette(colors=c("white","#EEA599"))(length(bk)/4),
           colorRampPalette(colors=c("#EEA599","#b34557"))(length(bk)/6),
           colorRampPalette(colors=c("#b34557","#9F0000"))(length(bk)/4)
           #colorRampPalette(colors=c("#E7B0A4","#b71c2c"))(length(bk)/2),
           #colorRampPalette(colors=c("#b71c2c","#990033"))(length(bk)/2),
           #colorRampPalette(colors=c("#9F0000","darkred"))(length(bk)/2),
           #colorRampPalette(colors=c("darkred","#660020"))(length(bk)/8)
         ),
         #color = colorRampPalette(c("#2464AB","#559BD4","lightblue","white","#E56F5E","#B22A2A","#660020"))(50),
         breaks = bk,
         legend_breaks = c(-2,0,1,2),
         gaps_col = gap,
         border_color = "black",
         legend = T,
         annotation_col = annotation_col,
         annotation_colors =ann_colors,
         fontsize_row = 10
)

####violin#data not normalize
m1_new<-as.data.frame(t(m1))
dat<-cbind(m1_new,info_new$Cluster)
dat$Sample<-rownames(dat)
colnames(dat)[9]<-"Cluster"
library(dplyr)
library(tidyr)
data <- dat %>% as.data.frame() %>%
  gather(key = Cell_type,value = Immune_core,-c(Sample,Cluster))

##########################################
library(ggplot2)
library(ggsci)
library(ggpubr)
library(scales)
#library(devtools)
#devtools::install_github("JanCoUnchained/ggunchained")
#install_github("JanCoUnchained/ggunchained")
library(ggunchained)
###自定义颜色
mypal=pal_simpsons(alpha = .6)(9)
mypal[c(1,7)]
show_col(mypal)
show_col(mypal[c(1,7)])
###自定义主题
mytheme <- theme(axis.text.x=element_text(size=12), 
                 axis.text.y=element_text(size=12), 
                 axis.title=element_text(size = 13), 
                 legend.text=element_text(size=12),
                 legend.title=element_text(size=12),
                 axis.line = element_line(size=0.7), 
                 panel.border = element_blank(),
                 panel.grid = element_blank())
##########################
###作图"#197EC099","#FED43999"
ggplot(data,aes(x=Cell_type,y = Immune_core,fill=Cluster))+
  geom_split_violin(trim = T,colour=NA)+
  geom_point(stat = 'summary',fun=mean,
             position = position_dodge(width = 0.2))+ #点的位置
  scale_fill_manual(values = c("#EEA599","#92B4C8"))+
  stat_summary(fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='black',
               width=0.01,size=0.5,
               position = position_dodge(width = 0.2))+#竖线的位置
  theme_bw()+
  stat_compare_means(aes(group = Cluster),
                     label = "p.signif",
                     method = "t.test",
                     label.y = max(data$Immune_core),
                     hide.ns = T)+
  mytheme+
  ylab("Immune_score")+xlab("Cell_Type")
