library(survival)
library(survminer)

#TF<-unique(c("SOX11","CEBPD","SOX4","EGR1","ETS1","RUNX3","ELF1"))

#TCGA-CGGA exp
load("D:/A-BS/mime1/TCGA_CGGA_input_exp.Rdata")

exp<-scale(t(exp_new[,4:14457]))
exp<-t(exp)
##tumor with T
TF1<-c("ETS1","SOX4")
TF2<-c("ETS1","EGR1")
TF3<-c("ETS1","SOX11")
TF4<-c("RUNX3","CEBPD")
##tumor
TF5<-c("CEBPD","EGR1")
TF6<-c("SOX11","SOX4")
#T
TF7<-c("ELF1","ETS1")
TF8<-c("ELF1","RUNX3")
TF9<-c("RUNX3","ETS1")

gene<-TF4
colnames(exp)<-trimws(colnames(exp))
sur<-exp_new[,2:3]

exp_new<-exp[,which(colnames(exp) %in% gene)]
exp_new<-as.data.frame(exp_new)
#colnames(exp_new)<-gene
########

exp_new$expr <- apply(exp_new, 1, sum)
cut_off<-mean(exp_new$expr)
risk_group <- ifelse(exp_new$expr > cut_off, "High", "Low")
#
sur$OS.time<-sur$OS.time/365
fit <- survfit(Surv(sur$OS.time, sur$OS) ~ risk_group, data = exp_new)
data<-cbind(sur,exp_new)

############
ggsurvplot(
  fit,
  data = data,
  size = 0.5,title="TCGA-CGGA",                 # 更改线条粗细
  # 配色方案，支持ggsci配色，自定义颜色，brewer palettes中的配色，等
  palette = c("#E7B800", "#2E9FDF"),
  conf.int = F,          # 可信区间
  pval = TRUE,              # log-rank P值，也可以提供一个数值
  pval.method = TRUE,       # 计算P值的方法，可参考https://rpkgs.datanovia.com/survminer/articles/Specifiying_weights_in_log-rank_comparisons.html
  log.rank.weights = "1",
  risk.table = T,        # 增加risk table
  risk.table.col = "strata",# risk table根据分组使用不同颜色
  legend.labs = c("High RUNX3+High CEBPD", "Low RUNX3+Low CEBPD"),    # 图例标签
  risk.table.height = 0.25, # risk table高度
  ggtheme = theme_classic2()      # 主题，支持ggplot2及其扩展包的主题
)



