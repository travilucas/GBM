library(randomForestSRC)
direction = 'forward'
################
gene<-c("AEBP1","ASF1A","DCC","HDAC5",
        "IL13RA2","OPHN1","PRPS1")

library(survival)
library(survminer)
#42669
load("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE42669/GSE42669exp_over.Rdata")
#108474
load("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE108474/GSE108474exp_over.Rdata")
#16011
load("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE16011/GSE16011exp_over.Rdata")
#7696
load("D:/A-BS/UCSC-GBM/GEO_Tumorsample/GSE7696/GSE7696exp_over.Rdata")

colnames(exp)<-trimws(colnames(exp))
sur<-exp[,1:2]

exp_new<-exp[,which(colnames(exp) %in% gene)]
exp_new<-as.data.frame(exp_new)

#Stepcox=forward
fit <- step(coxph(Surv(time = sur$OS.time, event = sur$OS) ~ ., exp_new), direction = direction)
rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因
  
est_dd2 <- exp[, c("OS.time", "OS", rid)]
#RSF
fit2 <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
             ntree = 1000, nodesize =5, #
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed =5201314)
#p<-predict(fit2, newdata = Dataset1)
rs <- predict(fit2, newdata = est_dd2)$predicted
#d1<-summary(coxph(Surv(OS.time, OS) ~ RS,rs$Dataset1))

# 分组
cut_off <- quantile(rs, probs = 0.5) #中位数
risk_group <- ifelse(rs > cut_off, "High", "Low")
#
fit <- survfit(Surv(sur$OS.time, sur$OS) ~ risk_group, data = exp_new)
data<-cbind(sur,exp_new)

ggsurvplot(
  fit,
  data = data,
  size = 1,title="GSE7696",                 # 更改线条粗细
  # 配色方案，支持ggsci配色，自定义颜色，brewer palettes中的配色，等
  palette = c("#E7B800", "#2E9FDF"),
  conf.int = TRUE,          # 可信区间
  pval = TRUE,              # log-rank P值，也可以提供一个数值
  pval.method = TRUE,       # 计算P值的方法，可参考https://rpkgs.datanovia.com/survminer/articles/Specifiying_weights_in_log-rank_comparisons.html
  log.rank.weights = "1",
  risk.table = T,        # 增加risk table
  risk.table.col = "strata",# risk table根据分组使用不同颜色
  legend.labs = c("High", "Low"),    # 图例标签
  risk.table.height = 0.25, # risk table高度
  ggtheme = theme_classic2()      # 主题，支持ggplot2及其扩展包的主题
)
