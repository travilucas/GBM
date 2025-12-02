load("D:/A-BS/mime1_new/汇报版/res.feature.all_result.Rdata")
#############
load("D:/A-BS/mime1_new/汇报版/hallmark_gene349.Rdata")

load("D:/A-BS/mime1/TCGA_CGGA_input_exp.Rdata")

Dataset1<-exp_new[129:427,]
Dataset2<-exp_new[1:128,]
#scale
d1<-scale(t(Dataset1[,4:14457]))
TCGA1<-cbind(Dataset1[,1:3],t(d1))

d2<-scale(t(Dataset2[,4:14457]))
TCGA2<-cbind(Dataset2[,1:3],t(d2))
list_train_vali_Data <- list(Train_set = TCGA1, Test_set = TCGA2)##Dataset1是训练数据集，其他数据集作为验证数据集
library(Mime1)

res1 <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Train_set,
                        list_train_vali_Data = list_train_vali_Data,
                        unicox.filter.for.candi = T,
                        unicox_p_cutoff = 0.05,
                        candidate_genes = g2,
                        mode = 'all',nodesize =5,seed = 5201314 )

save(res1,file = "D:/A-BS/mime1_new/汇报版/C_TCGA_CGGA_gene349.Rdata")
#########
load("D:/A-BS/mime1_new/汇报版/C_TCGA_CGGA_gene349.Rdata")

cindex_dis_all(res1,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35)
cindex_dis_select(res1,
                  model="StepCox[forward] + RSF",
                  order= names(list_train_vali_Data))#最有组合算法可以自由替换


res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Train_set,
                                              candidate_genes = g2,
                                              mode = "all",nodesize =5,seed = 5201314 )

save(res.feature.all,file = "D:/A-BS/mime1/res.feature326.all_result.Rdata")
#根据不同数据集中特定模型计算的风险评分绘制患者的生存曲线：
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res1, model_name = "StepCox[forward] + RSF",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)

all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res1,train_data = list_train_vali_Data[["Train_set"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res1,train_data = list_train_vali_Data[["Train_set"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res1,train_data = list_train_vali_Data[["Train_set"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")

auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)

roc_vis(all.auc.5y,
        model_name = "StepCox[forward] + RSF",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=5)

auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="StepCox[forward] + RSF",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))

unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res1,optimal.model = "StepCox[forward] + RSF",
                                   type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))

rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('GBM'),
                                         list_input_data = list_train_vali_Data)
#
HR_com(rs.glioma.lgg.gbm,
       res1,
       model_name="StepCox[forward] + RSF",
       dataset=names(list_train_vali_Data),
       type = "categorical")

cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('GBM'),
                                             list_input_data = list_train_vali_Data)

cindex_comp(cc.glioma.lgg.gbm,
            res1,
            model_name="StepCox[forward] + RSF",
            dataset=names(list_train_vali_Data))

auc.glioma.lgg.gbm.1<-cal_auc_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('GBM'),
                                           AUC_time=c(1),
                                           list_input_data = list_train_vali_Data)
auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name="StepCox[forward] + RSF",
         dataset=names(list_train_vali_Data))
####################
load("D:/A-BS/mime1/overlap3_geneset.Rdata")
gene<-gene3_over
##################
res.ici <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$Dataset1,
                                    list_train_vali_Data = list_train_vali_Data,
                                    candidate_genes = g2,
                                    methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                    seed = 5201314
                                    #cores_for_parallel = 60
)

auc_vis_category_all(res.ici,dataset = c("training","validation"),
                     order= c("training","validation"))

plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.ici,model_name = i,dataset = c("training","validation"),
                                   order= c("training","validation"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)

auc.other.pre <- cal_auc_previous_sig(list_train_vali_Data = list_train_vali_Data,seed = 5201314,
                                      train_data = list_train_vali_Data$Dataset1)
                                      #cores_for_parallel = 32)

load("./Example.cohort.Rdata")
load("./genelist.Rdata")
#核心特征选择
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Dataset1,
                                              candidate_genes = gene3,
                                              mode = "all",nodesize =5,seed = 5201314 )

core_feature_select(res.feature.all)
core_feature_rank(res.feature.all, top=20)
dataset_col<-c("#3182BDFF","#E6550DFF")
corplot <- list()
for (i in c(1:2)) {
  print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
                               dataset=names(list_train_vali_Data)[i],
                               color = dataset_col[i],
                               feature1="PSEN2",
                               feature2="WNT5B",
                               method="pearson"))
}
aplot::plot_list(gglist=corplot,ncol=2)

survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-core_feature_sur("HK2", 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                                        #color=c("blue","green"),
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)









