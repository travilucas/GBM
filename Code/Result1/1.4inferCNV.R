setwd("D:/A-BS/inferCNV/inferCNV")
load("D:/A-BS/GSE182109/result/182109_annotation_no_myeloid.Rdata")

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(infercnv)
library(RColorBrewer)
#######################################read files
expFile<-as.data.frame(sce1[["RNA"]]@counts)
groupFiles<-read.table("./182109_metadata_no_myeloid.txt",sep = "\t",row.names=1,header = F)
geneFile<-read.table("./hg38_gencode_v27.txt",sep = "\t",row.names = 1,header = F)
##########
#inferCNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("T cell","B cell"))  ## 这个取决于自己的分组信息里面的

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="try", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads = 4)

infercnv::plot_cnv(infercnv_obj,
                   plot_chr_scale = T,
                   output_filename = "better_plot",output_format = "pdf",
                   custom_color_pal = color.palette(c("#8DD3C7","white","#BC80BD"),c(2,2))
                   )
