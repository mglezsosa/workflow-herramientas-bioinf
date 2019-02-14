#! /usr/bin/env Rscript

# Enrichment analysis with GO and KEGG of the genes selected in the previous LASSO 
# regularized model.

# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
    stop("You must provide the name of the output pdf file.", call.=FALSE)
}

output = args[1]

setwd("R")
source("0_loadLibraries.R")
loadpkg("dplyr")
loadpkg("clusterProfiler")
loadpkg("org.Hs.eg.db")
loadpkg("biomaRt")

load("selectedGenesLasso.RData")

orgdb <- "org.Hs.eg.db"
biomart_dataset <- "hsapiens_gene_ensembl"
keggname <- "hsa"
mart <- biomaRt::useMart(biomart = "ensembl", dataset = biomart_dataset)
entrezsymbol <- biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
entrezsymbol$entrezgene <- as.character(entrezsymbol$entrezgene)

summarize_cp <- function(res) {
    summaries <- data.frame()
    for (ont in names(res)) {
        ontsum <- as.data.frame(res[[ont]])
        ontsum$ont <- ont
        summaries <- rbind(summaries, ontsum)
    }
    return(summaries)
}

convert_enriched_ids = function(res, entrezsymbol) {
    res = res %>% mutate(geneID = strsplit(as.character(geneID), "/")) %>% 
        tidyr::unnest(geneID) %>% 
        left_join(entrezsymbol, by = c(geneID = "entrezgene")) %>% 
        group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, ont) %>% 
        summarise(geneID = paste(geneID, collapse = "/"), symbol = paste(hgnc_symbol, collapse = "/"))
    return(res)
}

gene.ids <- data.frame(genes, stringsAsFactors = F) %>% left_join(., entrezsymbol, by = c("genes" = "hgnc_symbol")) %>% filter(!is.na(entrezgene))
mf = enrichGO(gene.ids$entrezgene, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH",
                          qvalueCutoff = 1, pvalueCutoff = 1)
cc = enrichGO(gene.ids$entrezgene,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH",
                qvalueCutoff = 1, pvalueCutoff = 1)
bp = enrichGO(gene.ids$entrezgene,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH",
                qvalueCutoff = 1, pvalueCutoff = 1)
kg = enrichKEGG(gene = gene.ids$entrezgene, organism = keggname, pvalueCutoff = 1,
                qvalueCutoff = 1, pAdjustMethod = "BH")

all <- list(mf= mf, cc = cc, bp = bp, kg = kg, summary = summarize_cp(all))

rm(mf)
rm(cc)
rm(bp)
rm(kg)
gc()

pdf(file = paste0("../", output, ".pdf"), width = 9, height = 4.5)
dotplot(all$kg, x="count", showCategory=10, color="qvalue", title = "KEGG enrichment")
dotplot(all$mf, x="count", showCategory=10, color="qvalue", title = "GO Molecular function")
dotplot(all$cc, x="count", showCategory=10, color="qvalue", title = "GO Cellular component")
dotplot(all$bp, x="count", showCategory=10, color="qvalue", title = "GO Biological process")
dev.off()

enrich_summary <- all$summary %>% arrange(p.adjust)
# enrich_summary <- convert_enriched_ids(enrich_summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
# write.csv(enrich_summary, paste0("../", output, ".csv"))
save(enrich_summary, file = "enrich_summary.RData")
