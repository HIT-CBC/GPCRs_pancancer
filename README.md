# DNA Methylation-Specific Analysis of G Protein-Coupled Receptor-Related Genes in Pan-Cancer
DNA methylation specific analysis of G-protein-coupled receptor-associated genes in pan-cancer. In this study, we analyzed the pan-cancer DNA methylation data from TCGA and GEO, and proposed a computational method to quantify the degree of specificity based on the level of DNA methylation of G protein-coupled receptor-related genes (GPCRs-related genes) and to identify specific GPCRs DNA methylation biomarkers (GRSDMs) in pan-cancer. Then, a ridge regression-based method was used to discover potential drugs through predicting the drug sensitivities of cancer samples.

![GA](https://user-images.githubusercontent.com/97509376/190980691-626ac5f7-4eef-4c47-8e28-47644f7ec8b6.png)



# Getting Started
R(https://www.r-project.org/) >=4.1.1 are required for the current codebase.
# R code
The codes used in this study are all in pan-cancer.R
+ ├─ Data preprocessing         
+ ├─ debatch		        //// There are three categories based on the normal sample size
+ ├─ Boruta				      //// Use of the Boruta algorithm
+ ├─ Specific gene classification model		   //// main code for look for specific genes in cancer
+ ├─ KEGG                    //// Differential gene pathway annotation
+ ├─ Predict Drugs                //// main code for ic50 prediction and drug target prediction

Please cite : Zhang M, Zhao J, Dong H, Xue W, Xing J, Liu T, Yu X, Gu Y, Sun B, Lu H, Zhang Y. DNA Methylation-Specific Analysis of G Protein-Coupled Receptor-Related Genes in Pan-Cancer. Genes (Basel). 2022 Jul 7;13(7):1213. doi: 10.3390/genes13071213. PMID: 35885996; PMCID: PMC9320183.
