library(ggplot2)
library(reshape2)
#####################TCGA-CGGA exp
load("D:/A-BS/mime1/TCGA_CGGA_input_exp.Rdata")
exp<-as.data.frame(scale(t(exp_new[,4:14457])))
exp<-as.data.frame(t(exp))
#######
TF1<-unique(c("SOX11","CEBPD","SOX4","EGR1"))
ICI<-c("ETS1","ELF1","RUNX3")
gene<-ICI
#######
colnames(exp)<-trimws(colnames(exp))
exp_new<-exp[,which(colnames(exp) %in% gene)]
exprSet<-exp
y <- exprSet[,which(colnames(exprSet) %in% TF1)]
data<-cbind(exp_new,y)
#####计算相关性矩阵
cor_matrix <- cor(data)
#############
# 计算相关性 p 值矩阵-计算全部的
cor_pvalue <- function(mat) {
  n <- ncol(mat)
  p_mat <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      test <- cor.test(mat[, i], mat[, j])  # 计算每对基因的 p 值
      p_mat[i, j] <- test$p.value
    }
  }
  return(p_mat)
}
#############
#只计算下三角部分
cor_pvalue <- function(mat) {
  n <- ncol(mat)
  p_mat <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {  # 只计算下三角部分
        test <- cor.test(mat[, i], mat[, j])
        p_mat[i, j] <- test$p.value
      }
    }
  }
  return(p_mat)
}
###########
p_matrix <- cor_pvalue(data)
colnames(p_matrix)<-colnames(cor_matrix)
rownames(p_matrix)<-rownames(cor_matrix)
####
# 将相关性矩阵和 p 值矩阵转换为长格式
cor_melted <- melt(cor_matrix)
colnames(cor_melted) <- c("Gene1", "Gene2", "Correlation")

p_melted <- melt(p_matrix)
colnames(p_melted) <- c("Gene1", "Gene2", "Pvalue")

# 合并相关性和 p 值数据
merged_data <- merge(cor_melted, p_melted, by = c("Gene1", "Gene2"))

# 添加显著性标记
merged_data$Significance <- ifelse(merged_data$Pvalue < 0.05, "*", "")

#绘制下三角部分的热图
## 只保留下三角部分
merged_data1 <- merged_data[as.numeric(factor(merged_data$Gene2)) >= as.numeric(factor(merged_data$Gene1)), ]

# 绘制热图
ggplot(merged_data1, aes(x = Gene1, y = Gene2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label =  Significance), color = "black", size = 4, vjust = 0.8) +  # 显示相关性值和显著性标记
  scale_fill_gradient2(low = "#6892c5", mid = "white", high = "#9F0000", midpoint = 0) +
  labs(title = "Pairwise Variable Correlation",
       x = "Variable",
       y = "Variable",
       fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed()  # 保持单元格为正方形


