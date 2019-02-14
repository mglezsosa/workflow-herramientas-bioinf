#! /usr/bin/env Rscript

# Cluster tumour samples by gene expression

# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
    stop("You must provide the name of the output pdf file and the number of clusters k.", call.=FALSE)
}

output = args[1]
k <- as.numeric(args[2])

setwd("R")
source("0_loadLibraries.R")
loadpkg("gplots")
loadpkg("cluster")

load("brca_rnaseq-tumour.RData")
load("diff.exp.genes.RData")
diff.exp.genes <- diff.exp.genes$Gene

expression_of_deg <- as.matrix(brca_rnaseq.tumour[rownames(brca_rnaseq.tumour)%in%diff.exp.genes,])

custom.boxplot <- function(mymatrix) {
    boxplot(mymatrix[,sample(1:ncol(mymatrix), size = 100)], xaxt = "n", main = "Boxplot of 100 random samples.")
}

custom.hist <- function(mymatrix) {
    hist(mymatrix[,sample(1:ncol(mymatrix), size = 100)], breaks = seq(min(mymatrix), max(mymatrix), length.out = 50), main = "Distribution after normalisation.")
}

hclust.silhouette <- function(dist.matrix, k.min, k.max, hclust.method){
    
    data.hclust <- hclust(dist.matrix, method = hclust.method)
    
    sapply(k.min:k.max, function (k) {
        ss <- silhouette(cutree(data.hclust, k), dist.matrix)
        mean(ss[, 3])
    })
}

sil.plots <- function(dist.matrix, k.min = 2, k.max = 10, title = NULL) {
    par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
    for (method in c("complete", "average", "single", "centroid")) {
        sil <- hclust.silhouette(as.dist(dist.matrix), k.min, k.max, method)
        plot(k.min:k.max, sil, type = "b", pch = 19, xlab = "Number of clusters k", 
             ylab="Silhouette index", main = paste("Hclust with", method, "method"))
        abline(v = which.max(sil) + k.min - 1, lty = 2)
    }
    mtext(title, outer = TRUE, cex = 1.5)
}

log2.expr <- log2(expression_of_deg + 1)
sp.dist <- 1 - cor(log2.expr, method = "spearman")
ps.dist <- 1 - cor(log2.expr, method = "pearson")
dist.matrix <- dist(t(log2.expr))

pdf(file = paste0("../", output, ".pdf"), width = 9, height = 7.5)
par(mfrow=c(2, 1))
custom.boxplot(log2.expr)
custom.hist(log2.expr)
hm.obj <- heatmap.2(log2.expr, trace = "none",
                    xlab = "Samples", ylab = "Genes",
                    labRow = FALSE, labCol = FALSE, margins = c(2, 2),
                    main = "Gene expression clustering")
sil.plots(sp.dist, k.max = 20, title = "Spearman")
sil.plots(ps.dist, k.max = 20, title = "Pearson")
sil.plots(dist.matrix, k.max = 20, title = "Euclidean")
dev.off()

hc <- hclust(as.dist(sp.dist), method = "complete")
groups <- cutree(hc, k = k)
table(groups)

save(groups, file = "groups.RData")
