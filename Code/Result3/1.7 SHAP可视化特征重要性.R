library(magrittr)
library(tidyverse)
library(fastshap)
library(shapviz)
library(randomForestSRC)
direction = 'forward'
################
gene<-c("AEBP1","ASF1A","DCC","HDAC5",
        "IL13RA2","OPHN1","PRPS1")

library(survival)
library(survminer)

load("D:/A-BS/mime1/TCGA_CGGA_input_exp.Rdata")

colnames(exp_new)<-trimws(colnames(exp_new))
sur<-exp_new[,2:3]
ID<-exp_new[,1]

exp1<-exp_new[,which(colnames(exp_new) %in% gene)]
exp1<-as.data.frame(exp1)

#Stepcox=forward
fit <- step(coxph(Surv(time = sur$OS.time, event = sur$OS) ~ ., exp1), 
            direction = direction)
rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因

exp2 <- exp_new[, c("OS.time", "OS", rid)]
#RSF
rfsrc_model <- rfsrc(Surv(OS.time,OS)~., data = exp2,
              ntree = 1000, nodesize =5, #
              splitrule = 'logrank',
              importance = T,
              proximity = T,
              forest = T,
              seed =5201314)

pred_wrapper <- function(model, newdata) {
  predict(model, newdata = newdata)$predicted
}

shap <- explain(rfsrc_model, X = exp1, pred_wrapper = pred_wrapper, nsim = 10, adjust = TRUE,
                shap_only = FALSE)

shv.global <- shapviz(shap)

sv_importance(shv.global)






