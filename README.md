# Dissecting Glioblastoma Risk Signatures in the Tumor Immune Microenvironment Based on Multi-dimensional Transcriptomics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)

## Research summary

Glioblastoma (GBM) is characterized by pronounced tumor heterogeneity and a complex immune microenvironment, contributing to poor patient survival outcomes. In this study, we comprehensively dissected the tumor microenvironment (TME) and uncovered potential molecular mechanisms by integrating **single-cell, bulk, and spatial transcriptomic data**. 

Hallmarks of malignancy and cell cycle regulatory pathways were consistently enriched across these modalities. We identified **seven hallmark-related prognostic signatures (HMsig)** using machine learning algorithms—namely **AEBP1, ASF1A, PRPS1, DCC, OPHN1, IL13RA2, and HDAC5**—whose importance in predicting patient outcomes was validated through SHAP algorithm analysis. 

Ligand-receptor (LR) interaction analysis revealed that interactions involving **OPHN1** were associated with poorer prognosis. Additionally, Immune checkpoint genes (ICG) **LAG3, PDCD1, and HAVCR2** were found to be substantially upregulated along the pseudotime trajectory of T-cell progression. Spatial transcriptomic analysis consistently demonstrated the existence of synergistic gene interactions, deciphering the immunomodulatory functions of GBM biomarkers in the TME.

## Repository Structure

The code and processed data are organized into sequential steps corresponding to the analysis workflow in the manuscript.

```text
GBM/
├── README.md                 # Project Overview & Instructions
├── LICENSE                   # MIT License
├── Code_summary/             # Analysis Source Code
│   ├── Step1_Single-cell RNA-seq Analysis Pipeline.R      # QC, Clustering, Annotation, InferCNV, Hallmarks of scRNA-seq
│   ├── Step2_Bulk RNA-seq Analysis Pipeline.R             # Batch Correction, Consensus Clustering, Hallmarks of bulk RNA-seq
│   ├── Step3_Model Construction and Validation.R          # ML Model (StepCox+RSF), SHAP, HMsig (7 genes)
│   ├── Step4_Cell-Cell Commnication Analysis.R             # L-R Interactions
│   ├── Step5_Trajectory Analysis.R                        # Monocle3 (T-cell trajectory)
│   ├── Step6_Regulon Analysis.R                           # Regulon Activity & Clinical Outcome
│   └── Step7_Spatial Transcriptomics Analysis.R           # Spatial Deconvolution (CARD) for Validation
└── Data/                     # Key Processed Data for Figures

Figure,Processed Data File(s),Description
Figure 1,"Data/Figure 1/sce_data.txt, scRNA_Tumor_hallmark.Rdata","scRNA-seq clustering, markers, and enrichment results."
Figure 2,"Data/Figure 2/consensus.Rdata, bulk_tumor_hallmark.Rdata",Consensus clustering (Cold/Hot) and bulk GSEA results.
Figure 3,Data/Figure 3/res.feature.all_result.Rdata,ML model features and HMsig risk scores.
Figure 4,"Data/Figure 4/Risk gene-LR.xlsx, string_interactions_short.tsv",LR communication probabilities and PPI network data.
Figure 5,"Data/Figure 5/pseudotime.txt, module_gene.txt",T-cell trajectory coordinates and pseudotime gene modules.
Figure 6,"Data/Figure 6/Step2_regulonTargetsInfo.tsv, 4Cell_InfoMatrix.Rdata",SCENIC regulon activity and TF-target interaction info.
Figure 7,Data/Figure 7/GBM5_CARD_celltype_prop.txt,Spatial deconvolution results (CARD proportions).
Supp. Figs,Data/Supplementary Figure 1-5/,"Validation data for sub-clusters, ML benchmarking, and spatial mapping."
