#! /usr/bin/env Rscript

setwd("R")
source("0_loadLibraries.R")
loadpkg("dplyr")

## Load data from module 1

load("brca_rnaseq.RData")
load("sample_data.RData")

# Subsetting normal and tumourous samples

are.tumours <- as.numeric(substr(colnames(brca_rnaseq), 14,15)) < 10
brca_rnaseq.tumour <- brca_rnaseq[, are.tumours]
brca_rnaseq.normal <- brca_rnaseq[, !are.tumours]

colnames(brca_rnaseq.normal) <- substr(colnames(brca_rnaseq.normal), 1,12)
colnames(brca_rnaseq.tumour) <- substr(colnames(brca_rnaseq.tumour), 1,12)
brca_rnaseq.normal <- brca_rnaseq.normal[,!duplicated(colnames(brca_rnaseq.normal))]
brca_rnaseq.tumour <- brca_rnaseq.tumour[,!duplicated(colnames(brca_rnaseq.tumour))]

# Subsetting according to cancer subtype

tnbc_samples <- sample_data %>% dplyr::filter(`ER Status` == "Negative" & `PR Status` == "Negative" & `HER2 Final Status` == "Negative" & `PAM50 mRNA` != "Luminal A")
tnbc_barcodes <- tnbc_samples$`Complete TCGA ID`

luminalA_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Luminal A")
luminalA_barcodes <- luminalA_samples$`Complete TCGA ID`

luminalB_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Luminal B")
luminalB_barcodes <- luminalB_samples$`Complete TCGA ID`

basal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Basal-like")
basal_barcodes <- basal_samples$`Complete TCGA ID`

her2.enrich_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "HER2-enriched")
her2.enrich_barcodes <- her2.enrich_samples$`Complete TCGA ID`

brca_rnaseq.luminalA <- brca_rnaseq.tumour[,colnames(brca_rnaseq.tumour) %in% luminalA_barcodes]
brca_rnaseq.luminalB <- brca_rnaseq.tumour[,colnames(brca_rnaseq.tumour) %in% luminalB_samples]
brca_rnaseq.tnbc <- brca_rnaseq.tumour[,colnames(brca_rnaseq.tumour) %in% tnbc_barcodes]
brca_rnaseq.basal <- brca_rnaseq.tumour[,colnames(brca_rnaseq.tumour) %in% basal_barcodes]
brca_rnaseq.her2.enrich <- brca_rnaseq.tumour[,colnames(brca_rnaseq.tumour) %in% her2.enrich_barcodes]

save(brca_rnaseq.luminalA, file = "brca_rnaseq-luminalA.RData")
save(brca_rnaseq.luminalB, file = "brca_rnaseq-luminalB.RData")
save(brca_rnaseq.tnbc, file = "brca_rnaseq-tnbc.RData")
save(brca_rnaseq.basal, file = "brca_rnaseq-basal.RData")
save(brca_rnaseq.her2.enrich, file = "brca_rnaseq-her2-enrich.RData")
save(brca_rnaseq.normal, file = "brca_rnaseq-normal.RData")
save(brca_rnaseq.tumour, file = "brca_rnaseq-tumour.RData")
