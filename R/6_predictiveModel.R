#! /usr/bin/env Rscript

# Retrieve clinical data and build a model for fixed-time survival 
# predictions, integrating gene expression and clinical data.

setwd("R")
source("0_loadLibraries.R")
loadpkg("glmnet")
loadpkg("TCGAretriever")
loadpkg("dplyr")
loadpkg("pROC")
set.seed(1953)

load("diff.exp.genes.RData")
load("brca_rnaseq-tumour.RData")

months <- 24

# Get data
cdat.tcga <- suppressWarnings(get_clinical_data("brca_tcga_all"))
# Order data by barcode
cdat.tcga <- cdat.tcga[with(cdat.tcga, order(CASE_ID)),]
# Trim the barcode
cdat.tcga$CASE_ID <- substr(cdat.tcga$CASE_ID, 1,12)
# Remove duplicates
cdat.tcga <- cdat.tcga[!duplicated(cdat.tcga$CASE_ID),]
# Select those samples who are in target type (all tumour by default)
cdat.tcga <- cdat.tcga[cdat.tcga$CASE_ID %in% colnames(brca_rnaseq.tumour),]

# Select differentially expressed genes
filtered.exp.genes <- brca_rnaseq.tumour[rownames(brca_rnaseq.tumour) %in% diff.exp.genes$Gene,]
filtered.exp.genes <- log2(filtered.exp.genes + 1)

rm(brca_rnaseq.tumour)

# Clean, transform and select clinical-pathologic variables
cdat.tcga %>% mutate(
        HER2 = as.factor(ifelse("Positive" == ifelse(
            IHC_HER2 == "Equivocal" | IHC_HER2 == "Indeterminate" | IHC_HER2 == "",
            ifelse(HER2_FISH_STATUS == "", NA, HER2_FISH_STATUS), 
            IHC_HER2), 1, 0)),
        ER = as.factor(ifelse("Positive" == ifelse(
            ER_STATUS_BY_IHC == "Equivocal" | ER_STATUS_BY_IHC == "Indeterminate" | ER_STATUS_BY_IHC == "",
            NA, 
            ER_STATUS_BY_IHC), 1, 0)),
        PR = as.factor(ifelse("Positive" == ifelse(
            PR_STATUS_BY_IHC == "Equivocal" | PR_STATUS_BY_IHC == "Indeterminate" | PR_STATUS_BY_IHC == "",
            NA, 
            PR_STATUS_BY_IHC), 1, 0)),
        MENOPAUSE = as.factor(
            ifelse(grepl("^Post", MENOPAUSE_STATUS), 1,
            ifelse(grepl("^Pre|^Peri", MENOPAUSE_STATUS), 0, NA))),
        PT = as.factor(
            ifelse(grepl("^T1", AJCC_TUMOR_PATHOLOGIC_PT), 0, 1)),
        PN = as.factor(
            ifelse(grepl("^N1", AJCC_NODES_PATHOLOGIC_PN), 0, 1)),
        PM = as.factor(
            ifelse(grepl("^M1", AJCC_METASTASIS_PATHOLOGIC_PM), 0, 1)),
        AGE = as.numeric(AGE),
        OS_MONTHS = as.numeric(OS_MONTHS),
        DFS_MONTHS = as.numeric(DFS_MONTHS)
    ) %>% select(c(CASE_ID, OS_MONTHS, OS_STATUS, DFS_MONTHS, DFS_STATUS, AGE, 
                   RACE, HER2, ER, PR, PM, PN, PT, MENOPAUSE)) -> cdat.tcga

# Transform continuous survival data to fixed-time survival (`months` months)
# To predict living status
cdat.tcga %>% mutate(DEAD = ifelse(OS_STATUS == "DECEASED",
                                     ifelse(OS_MONTHS < months, 1, 0),
                                     ifelse(OS_MONTHS >= months, 0, NA))) %>% 
    filter(!is.na(DEAD)) -> clin.os

# To predict progression
cdat.tcga %>% filter(DFS_STATUS == "DiseaseFree" || DFS_STATUS == "Recurred/Progressed") %>% 
    mutate(RECURRED = as.factor(ifelse(DFS_STATUS == "Recurred/Progressed",
                               ifelse(DFS_MONTHS < months, 1, 0),
                               ifelse(DFS_MONTHS >= months, 0, NA)))) %>% 
    filter(!is.na(RECURRED)) -> clin.dfs

# Rowwise removal of missing values
clin.dfs <- na.omit(clin.dfs)
clin.os <- na.omit(clin.os)

# Merge expression data with clinical
clin.exp.os <- cbind(clin.os, t(filtered.exp.genes[,colnames(filtered.exp.genes) %in% clin.os$CASE_ID]))
clin.exp.dfs <- cbind(clin.dfs, t(filtered.exp.genes[,colnames(filtered.exp.genes) %in% clin.dfs$CASE_ID]))

# Without expression data
dfs.design.matrix <- model.matrix(RECURRED ~ AGE + RACE + HER2 + ER + PR + PM + 
                                  PN + PT + MENOPAUSE, data = clin.dfs)[,-1]
os.design.matrix <- model.matrix(DEAD ~ AGE + RACE + HER2 + ER + PR + PM + 
                                     PN + PT + MENOPAUSE, data = clin.os)[,-1]
# Without expression data
dfs.exp.design.matrix <- model.matrix(RECURRED ~AGE + RACE + HER2 + ER + PR + PM + 
                                          PN + PT + MENOPAUSE, data = clin.exp.dfs)[,-1]
os.exp.design.matrix <- model.matrix(DEAD ~AGE + RACE + HER2 + ER + PR + PM + 
                                         PN + PT + MENOPAUSE, data = clin.exp.os)[,-1]

nested.cv.glmnet <- function(x, y, outer.fold = 10, inner.fold = 10, alpha = 1, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    lvls <- levels(y)
    n <- nrow(x)
    prop <- n%/%outer.fold
    newseq <- rank(runif(n))
    k <- as.factor((newseq - 1)%/%prop + 1)
    vec.auc <- vector(length = outer.fold)
    vec.models <- vector(length = outer.fold)
    for (i in seq(outer.fold)) {
        fit <- cv.glmnet(x[k != i,], y[k != i], family = "binomial", 
                         type.measure = "auc", alpha = alpha, nfolds = inner.fold)
        vec.models[i] <- list(fit)
        prediction <- predict(fit$glmnet.fit, newx = x[k == i,], s = fit$lambda.1se, type = "class")
        print(prediction)
        if(length(levels(y[k == i])) == 1) {
            vec.auc[i] = NA
            next
        }
        vec.auc[i] <- pROC::auc(pROC::roc(y[k == i], as.numeric(prediction), levels = lvls))
    }
    avg.auc <- mean(vec.auc, na.rm = T)
    list(avgauc = avg.auc, aucs = vec.auc, models = vec.models)
}

res.dfs.lasso <- nested.cv.glmnet(dfs.design.matrix, as.factor(clin.dfs$RECURRED), 8, 8, 1, seed = 1234)
res.os.lasso <- nested.cv.glmnet(os.design.matrix, as.factor(clin.os$DEAD), 8, 8, seed = 1234)
res.exp.dfs.lasso <- nested.cv.glmnet(dfs.exp.design.matrix, as.factor(clin.exp.dfs$RECURRED), 8, 8, 1, seed = 1234)
res.exp.os.lasso <- nested.cv.glmnet(os.exp.design.matrix, as.factor(clin.exp.os$DEAD), 8, 8, 1, seed = 1234)

res.dfs.ridge <- nested.cv.glmnet(dfs.design.matrix, as.factor(clin.dfs$RECURRED), 8, 8, 0, seed = 1234)
res.os.ridge <- nested.cv.glmnet(os.design.matrix, as.factor(clin.os$DEAD), 8, 8, 0, seed = 1234)
res.exp.dfs.ridge <- nested.cv.glmnet(dfs.exp.design.matrix, as.factor(clin.exp.dfs$RECURRED), 8, 8, 0, seed = 1234)
res.exp.os.ridge <- nested.cv.glmnet(os.exp.design.matrix, as.factor(clin.exp.os$DEAD), 8, 8, 0, seed = 1234)
