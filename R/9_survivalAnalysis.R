#! /usr/bin/env Rscript

# Perform the survival analysis for the found groups.

# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
    stop("You must provide the name of the output pdf file.", call.=FALSE)
}

output = args[1]

setwd("R")
source("0_loadLibraries.R")
loadpkg("survminer")
loadpkg("survival")
loadpkg("TCGAretriever")

load("groups.RData")

cdat.tcga <- suppressWarnings(get_clinical_data("brca_tcga_all"))
cdat.tcga$CASE_ID <- substr(cdat.tcga$CASE_ID, 1, 12)
cdat.tcga <- cdat.tcga[order(cdat.tcga$CASE_ID),]
cdat.tcga <- cdat.tcga[!duplicated(cdat.tcga$CASE_ID),]
survival.data <- cdat.tcga[cdat.tcga$CASE_ID %in% names(groups),]
survival.data$group <- factor(groups[names(groups) %in% survival.data$CASE_ID])

os.mapping <- c("LIVING"=0, "DECEASED"=1)
dfs.mapping <- c("DiseaseFree"=0, "Recurred/Progressed"=1)

survival.data <- survival.data %>%
    mutate(Survival_time = as.numeric(OS_MONTHS),
           Disease_free_time = as.numeric(DFS_MONTHS),
           Death = os.mapping[OS_STATUS],
           Recurrence = dfs.mapping[DFS_STATUS])

surv.fit.os <- survfit(formula = Surv(Survival_time, Death) ~ group, data = survival.data)
surv.fit.dfs <- survfit(formula = Surv(Disease_free_time, Recurrence) ~ group, data = survival.data)

pdf(file = paste0("../", output, ".pdf"), width = 9, height = 7.5)
ggsurvplot(
    surv.fit.os,
    conf.int = T, 
    censor = T, 
    cex.axis=3.0, 
    cex.lab=3.0, 
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    pval=T,
    title = "Overall survival",
    cumevents = T,
    xlab = "Time in months",
    break.time.by = 12)

ggsurvplot(
    surv.fit.dfs,
    conf.int = T, 
    censor = T, 
    cex.axis=3.0, 
    cex.lab=3.0, 
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    pval=T,
    title = "Disease free survival",
    cumevents = T,
    xlab = "Time in months",
    break.time.by = 12)
dev.off()