#! /usr/bin/env Rscript

# Build a predictor by finding a gene signature to classify cancer subtype with gene expression data

# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
    stop("You must provide the name of the txt output file.", call.=FALSE)
}

output = args[1]

setwd("R")
source("0_loadLibraries.R")
loadpkg("glmnet")
loadpkg("dplyr")
loadpkg("pROC")
loadpkg("doParallel")
set.seed(1953)

# may not be compatible
doParallel::registerDoParallel(4)

load("diff.exp.genes.RData")
load("brca_rnaseq-tumour.RData")
load("sample_data.RData")

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
        fit <- cv.glmnet(x[k != i,], y[k != i], family = "multinomial", 
                         type.measure = "class", alpha = alpha, nfolds = inner.fold,
                         parallel = T)
        vec.models[i] <- list(fit)
        prediction <- predict(fit$glmnet.fit, newx = x[k == i,], s = fit$lambda.min, type = "class")
        if(length(levels(y[k == i])) == 1) {
            vec.auc[i] = NA
            next
        }
        mlt.roc <- pROC::multiclass.roc(as.numeric(factor(y[k == i])), as.numeric(factor(prediction)))
        vec.auc[i] <- pROC::auc(mlt.roc)
    }
    avg.auc <- mean(vec.auc, na.rm = T)
    list(avgauc = avg.auc, aucs = vec.auc, models = vec.models)
}

sample_data <- sample_data[order(sample_data$`Complete TCGA ID`),]

flt.cdata <- sample_data[
    sample_data$`PAM50 mRNA` != "NA" & sample_data$`PAM50 mRNA` != "Normal-like" &
        sample_data$`Complete TCGA ID` %in% colnames(brca_rnaseq.tumour),
    c("Complete TCGA ID", "PAM50 mRNA")]

flt.exp.genes <- log2(brca_rnaseq.tumour[
    rownames(brca_rnaseq.tumour) %in% diff.exp.genes$Gene,
    colnames(brca_rnaseq.tumour) %in% flt.cdata$`Complete TCGA ID`] + 1)

cmpl.data <- cbind(flt.cdata, t(flt.exp.genes))
rownames(cmpl.data) <- flt.cdata$`Complete TCGA ID`
cmpl.data$`Complete TCGA ID` <- NULL

design.mtx <- model.matrix(`PAM50 mRNA` ~ ., data = cmpl.data)[,-1]

message("Lasso fitting")
res.lasso.min <- nested.cv.glmnet(x = design.mtx, y = cmpl.data$`PAM50 mRNA`, 
                              outer.fold = 10, inner.fold = 10, alpha = 1, seed = 1953)
message("Ridge fitting")
res.ridge.min <- nested.cv.glmnet(x = design.mtx, y = cmpl.data$`PAM50 mRNA`, 
                              outer.fold = 10, inner.fold = 10, alpha = 0, seed = 1953)
message("Elastic-net fitting")
res.elasticnet.min <- nested.cv.glmnet(x = design.mtx, y = cmpl.data$`PAM50 mRNA`, 
                                   outer.fold = 10, inner.fold = 10, alpha = 0.5, seed = 1953)

message(paste("Average AUC for lasso:", res.lasso.min$avgauc))
message(paste("Average AUC for ridge:", res.ridge.min$avgauc))
message(paste("Average AUC for elastic net:", res.elasticnet.min$avgauc))

co <- coef(res.lasso.min$models[[1]], s=res.lasso.min$models[[1]]$lambda.min, exact=TRUE)

# sum(co$`Luminal A` != 0)
# sum(co$`Luminal B` != 0)
# sum(co$`HER2-enriched` != 0)
# sum(co$`Basal-like` != 0)

selected.genes.lumA <- rownames(co$`Luminal A`)[as.vector(co$`Luminal A` != 0)]
selected.genes.lumB <- rownames(co$`Luminal B`)[as.vector(co$`Luminal B` != 0)]
selected.genes.her2 <- rownames(co$`HER2-enriched`)[as.vector(co$`HER2-enriched` != 0)]
selected.genes.basal <- rownames(co$`Basal-like`)[as.vector(co$`Basal-like` != 0)]

gene.frequencies <- table(c(selected.genes.basal[-1], 
selected.genes.her2[-1], selected.genes.lumA[-1], selected.genes.lumB[-1]))

genes <- names(gene.frequencies)
write(genes, file = paste0("../", output, ".txt"))

save(genes, file = "selectedGenesLasso.RData")
