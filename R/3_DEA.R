#! /usr/bin/env Rscript
## DEA


# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)<5) {
    stop("Two first arguments must be cancer types ('Basal', 'TNBC', 'LumA', 'LumB', 'HER2+' 'Tumour' or 'Normal'), the third must be a valid p-value, the fourth must be a valid fold-change value and the fifth must be a valid output relative path for pdf and csv files.", call.=FALSE)
}

if (args[1] == args[2]) {
    stop("Cannot perform the DEA with the same subsets.")
}

cancertype1 = args[1]
cancertype2 = args[2]

# Load the expression data
setwd("R")

source("0_loadLibraries.R")
loadpkg("dplyr")
loadpkg("limma") # Differential gene expression analysis (DEA)
loadpkg("edgeR")
loadpkg("calibrate") # To label the volcano plot

# Pick the datasets to analyse

switch(cancertype1,
       "Basal"={
           load("brca_rnaseq-basal.RData")
           d1 = brca_rnaseq.basal
       },
       "TNBC"={
           load("brca_rnaseq-tnbc.RData")
           d1 = brca_rnaseq.tnbc
       },
       "LumA"={
           load("brca_rnaseq-luminalA.RData")
           d1 = brca_rnaseq.luminalA
       },
       "LumB"={
           load("brca_rnaseq-luminalB.RData")
           d1 = brca_rnaseq.luminalB
       },
       "Normal"={
           load("brca_rnaseq-normal.RData")
           d1 = brca_rnaseq.normal
       },
       "HER2+"={
           load("brca_rnaseq-her2-enrich.RData")
           d1 = brca_rnaseq.her2.enrich
       },
       "Tumour"={
           load("brca_rnaseq-tumour.RData")
           d1 = brca_rnaseq.tumour
       },
       stop(paste("Unrecognized cancer type:", cancertype1))
)

switch(cancertype2,
       "Basal"={
           load("brca_rnaseq-basal.RData")
           d2 = brca_rnaseq.basal
       },
       "TNBC"={
           load("brca_rnaseq-tnbc.RData")
           d2 = brca_rnaseq.tnbc
       },
       "LumA"={
           load("brca_rnaseq-luminalA.RData")
           d2 = brca_rnaseq.luminalA
       },
       "LumB"={
           load("brca_rnaseq-luminalB.RData")
           d2 = brca_rnaseq.luminalB
       },
       "Normal"={
           load("brca_rnaseq-normal.RData")
           d2 = brca_rnaseq.normal
       },
       "HER2+"={
           load("brca_rnaseq-her2-enrich.RData")
           d2 = brca_rnaseq.her2.enrich
       },
       "Tumour"={
           load("brca_rnaseq-tumour.RData")
           d2 = brca_rnaseq.tumour
       },
       stop(paste("Unrecognized cancer type:", cancertype2))
)

cat("Differential Expression analysis of breast cancer subtypes:", cancertype1, "vs", cancertype2, "\n")


# We Combine the two matrices for gene differential expression analysis (DEA). Further preprocessing included the removal of expression estimates with counts in less than 20% of cases.
rnaseq.for.de <- cbind(d1, d2)
counts = rnaseq.for.de[apply(rnaseq.for.de,1,function(x) sum(x==0))<ncol(rnaseq.for.de)*0.8,]

# Create a design matrix thar contains the RNA samples that are applied to each category
df <- rbind(
    data.frame(sample = colnames(d1), status = rep(1, length(colnames(d1)))),
    data.frame(sample = colnames(d2), status = rep(0, length(colnames(d2))))
)
design <- model.matrix(~ status, data = df)


dge <- DGEList(counts=counts)
A <- rowSums(dge$counts)
isexpr <- A > 100 # Keeping genes with total counts more than 100.
dge <- calcNormFactors(dge)
v <- voom(dge[isexpr,], design, plot=FALSE)

# find genes differentially expression between the two groups of samples combined above
fit <- lmFit(v, design)
fit <- eBayes(fit)

lfc = as.double(args[4])
pval = as.double(args[3])

diff.exp.df <- topTable(fit, coef = "status", n = Inf, sort = "p", p = pval) # Positive log-fold-changes mean higher expression in d1
diff.exp.df$gene.name <- rownames(diff.exp.df)
# With the code above, perform the DEA. First, we convert the read counts to log2-cpm, with associated weights, ready for linear modelling. As read counts follow a negative binomial distribution, which has a mathematical theory less tractable than that of the normal distribution, RNAseq data was normalised with the voom methodology. The voom method estimates the mean-variance of the log-counts and generates a precision weight for each observation. This way, a comparative analysis can be performed with all bioinformatic workflows originally developed for microarray analyses ([see this](https://bioconductor.org/packages/3.7/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html) and this paper Charity W Law et al. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. In: Genome biology 15.2 (Jan. 2014), R29–R29).

if (!nrow(diff.exp.df)) {
    stop(paste0("No genes filfill the given restrictions: pvalue=", pval, ", logfoldchange=", lfc, "."))
}

# Save the differential expression dataframe and types
save(diff.exp.df, cancertype1, cancertype2, file = "diff.exp.RData")

# Output, Volcano plot
tab = data.frame(logFC = diff.exp.df$logFC, negLogPval = -log10(diff.exp.df$adj.P.Val))
tab2 = data.frame(logFC = diff.exp.df$logFC, negLogPval = -log10(diff.exp.df$adj.P.Val), Gene=diff.exp.df$gene.name)
diff.exp.genes <- filter(tab2, abs(logFC) > lfc & negLogPval > -log10(pval))

save(diff.exp.genes, file = "diff.exp.genes.RData")
write.csv(filter(tab2, abs(logFC) > lfc & negLogPval > -log10(pval)), paste0("../", args[5], ".csv")) # write output

pdf(file = paste0("../", args[5], ".pdf"), width = 9, height = 4.5)
par(mar = c(5, 4, 4, 5))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
title(main = paste(cancertype1, "vs", cancertype2))
points(tab[(abs(tab$logFC) > lfc), ], pch = 16, cex = 0.8, col = "orange") 
points(tab[(tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "green") 
points(tab[(abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "red") 
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
with(subset(tab2, negLogPval > -log10(pval) & abs(logFC)>lfc), textxy(logFC, negLogPval, labs=Gene, cex=.4))
dev.off()