####################################################
load("D:/A-BS/CellChat/Tumor&subtype_Cellchat.Rdata")
TB_target <- as.data.frame(c("SPP1_CD44",
                              "SPP1_ITGA4_ITGB1",
                              "PTN_PTPRZ1",
                              "PTN_NCL",
                              "MIF_CD74_CXCR4",
                              "MIF_CD74_CD44"))


ICI<-as.data.frame(c("LGALS9_HAVCR2",
                     "CD274_PDCD1",
                     "PDCD1LG2_PDCD1"
                     ))

r<-as.data.frame(c("IL21_IL13RA2",#581
                   "IL7_IL4R_IL13RA2"#577 579
                   
))

#a<-as.data.frame(cellchat@net[["prob"]])


colnames(ICI) <- 'interaction_name'

MAC_source <- as.data.frame(c("SPP1_CD44",
                            
                              "MIF_CD74_CXCR4",
                       "MIF_CD74_CD44","IL1B_IL1R2","C3_ITGAX_ITGB2","C3_ITGAM_ITGB2",
                       "LGALS9_CD44","LGALS9_CD45"
  
))
colnames(MAC_source) <- 'interaction_name'

myeloid_source<-as.data.frame(c("SPP1_CD44",
                                "SPP1_ITGA4_ITGB1",
                                "PTN_PTPRZ1",
                                "PTN_NCL",
                                "MIF_CD74_CXCR4","THBS1_CD47",
                                "MIF_CD74_CD44","LGALS9_CD44","LGALS9_CD45"
  
))
colnames(myeloid_source) <- 'interaction_name'

B_source<-as.data.frame(c("SPP1_CD44",
                            "SPP1_ITGA4_ITGB1",
                            "PTN_NCL","MDK_NCL","MDK_ITGA4_ITGB1","CD99_CD99",
                            "MIF_CD74_CXCR4",
                            "MIF_CD74_CD44","APP_CD74"
                            
))
colnames(B_source) <- 'interaction_name'

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
TB1111 = netVisual_bubble(cellchat, sources.use =c(2,3,7,10), #tumor cell
                     targets.use = c(5,9), #T B cell 
                     remove.isolate = FALSE,
                     #signaling=c("MIF","PTN","SPP1","C3","IL1B"))+
                     pairLR.use = ICI)+
  theme(legend.position = 'right',
        legend.key.width = unit(1,'cm'),
        legend.margin = margin(0.5,0.5,0,0,'cm'),
        plot.margin=unit(c(2.5, 2.5, 2.5, 2.5),'cm'))+
  coord_cartesian(clip = 'off')

cols<-c("#D6251F","#3778AD","#B699C6","#4EA74A","#202D70","#E57FB0","#A15528","#8F4C9A")
T_col<-c("#3778AD","#B699C6","#4EA74A","#202D70")
Myeloid_col<-c("#8F4C9A","#E57FB0","#A15528")
col<-c("#E89118",'#56A8D7')

library(jjAnno)
p<-annoSegment(object = TB1111,
                  annoPos = 'top', 
                  aesGroup = T,
                  aesGroName = 'source',
                  yPosition =9.7,
                  segWidth = 0.8,
                  pCol="#D6251F",
                  addText=T,
                  textSize = 5,
                  textCol = c("black","black"))

ggsave("D:/A-BS/CellChat/LR-dotplot/myeloid_TB_LR_new.pdf", p, width = 6.2, height = 4.5) 
