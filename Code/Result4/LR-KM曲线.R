library(survminer)
library(survival)
##########
load("D:/A-BS/mime1/TCGA_CGGA_input_exp.Rdata")

exp<-scale(t(exp_new[,4:14457]))
exp<-t(exp)

#exp<-exp_new[,4:14457]
################
#42669
load("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE42669/GSE42669exp_over.Rdata")
#108474
load("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/GSE108474exp_over.Rdata")

#7696
load("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE7696/GSE7696exp_over.Rdata")
#################

gene<-c("EFNB1","OPHN1")
colnames(exp)<-trimws(colnames(exp))
sur<-exp[,1:2]

exp_new<-exp[,which(colnames(exp) %in% gene)]
exp_new<-as.data.frame(exp_new)

#method
# 考虑基因乘积效应
exp_new$interaction <- exp_new$EFNB1 * exp_new$OPHN1
sur_exp<-cbind(sur,exp_new)
res.cut <- surv_cutpoint(sur_exp, time="OS.time", event="OS", variables="interaction")
# 
sur_exp$risk_group <- ifelse(sur_exp$interaction > res.cut[["cutpoint"]][["cutpoint"]], "High", "Low")
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data=sur_exp)
##
ggsurvplot(
  fit,
  data = sur_exp,
  size = 0.5,title="GSE108474",                 # 更改线条粗细
  # 配色方案，支持ggsci配色，自定义颜色，brewer palettes中的配色，等
  palette = c("#E7B800", "#2E9FDF"),
  conf.int = F,          # 可信区间
  pval = TRUE,              # log-rank P值，也可以提供一个数值
  pval.method = TRUE,       # 计算P值的方法，可参考https://rpkgs.datanovia.com/survminer/articles/Specifiying_weights_in_log-rank_comparisons.html
  log.rank.weights = "1",
  risk.table = T,        # 增加risk table
  risk.table.col = "strata",# risk table根据分组使用不同颜色
  legend.labs = c("High OPHN1+High EFNB1", "Low OPHN1+Low EFNB1"),    # 图例标签
  risk.table.height = 0.25, # risk table高度
  ggtheme = theme_classic2()      # 主题，支持ggplot2及其扩展包的主题
)
