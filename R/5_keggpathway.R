#! /usr/bin/env Rscript


# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
        stop("provide a comma separated string with path ids.", call.=FALSE)
}

pathids <- unlist(strsplit(args[1], split = ","))

setwd("R")
source("0_loadLibraries.R")

loadpkg("dplyr")
loadpkg("pathview")

viewkeggpath <- function(dea.df, path, enrichment){
        # Get the genes in path from the enrichment summary
        genesymbols <- strsplit(as.character(enrichment[which(enrichment$ID == path),12]), split = "/")[[1]]
        geneids <- strsplit(as.character(enrichment[which(enrichment$ID ==path),11]), split = "/")[[1]]
        ll <- data.frame(geneid=geneids, gene.name=genesymbols)
        # dataframe with the genes and the log2FC
        kk <- filter(dea.df, gene.name %in% genesymbols) %>% dplyr::select(gene.name, logFC)
        kk <- kk %>%right_join(ll,by="gene.name")
        gene.vector <- kk$logFC
        names(gene.vector) <- kk$geneid
        pathview(
                gene.data  = gene.vector, 
                pathway.id = path, 
                species = "hsa", 
                limit = list(gene=max(abs(gene.vector)), cpd=1))
}

load("enrich_summary.RData")
load("diff.exp.RData")
setwd("..")

for (path in pathids) {
        viewkeggpath(dea.df = diff.exp.df, path=path, enrichment = enrich_summary)
}
