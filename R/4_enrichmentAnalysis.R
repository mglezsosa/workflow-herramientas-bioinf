#! /usr/bin/env Rscript


# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
        stop("pvalue, log fold change and output filename for the pdf and csv files must be passed as arguments.", call.=FALSE)
}

pval = as.double(args[1])
lfc = as.double(args[2])
enrichmentfilename = args[3]

setwd("R")
source("0_loadLibraries.R")

loadpkg("biomaRt")
loadpkg("dplyr")
loadpkg("org.Hs.eg.db")
loadpkg("clusterProfiler")

orgdb = "org.Hs.eg.db"
biomart_dataset = "hsapiens_gene_ensembl"
keggname = "hsa"
mart = biomaRt::useMart(biomart = "ensembl", dataset = biomart_dataset)
# entrez = biomaRt::getBM(attributes = c("refseq_mrna", "entrezgene"), mart = mart)
# entrez$entrezgene = as.character(entrez$entrezgene)
entrezsymbol = biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
entrezsymbol$entrezgene = as.character(entrezsymbol$entrezgene)


summarize_cp = function(res, comparison) {
    summaries = data.frame()
    for (ont in names(res)) {
        ontsum = as.data.frame(res[[ont]])
        ontsum$ont = ont
        summaries = rbind(summaries, ontsum)
    }
    summaries$comparison = comparison
    return(summaries)
}

enrich_cp = function(res, comparison, pval, lfc, type="over") {
    res = res %>% data.frame()  %>% left_join(entrezsymbol, by = "hgnc_symbol") %>% filter(!is.na(entrezgene))
    # universe = brcaData@GISTIC@AllByGene$Gene.Symbol
    if(type=="all"){
        res <- res %>% filter(abs(logFC) > lfc & adj.P.Val < pval) # lfc and pval threshold defined above in the volcano plot
        genes = res$entrezgene
        
        mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        cc = enrichGO(genes,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        bp = enrichGO(genes,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1,
                        qvalueCutoff = 1, pAdjustMethod = "BH")
        all = list(mf = mf, cc = cc, bp = bp, kg = kg)
        all[["summary"]] = summarize_cp(all, comparison)
        return(all)
    }
    if(type=="over"){
        res.over <- res %>% filter(logFC > lfc  & adj.P.Val < pval)
        genes = res.over$entrezgene
        mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        cc = enrichGO(genes,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        bp = enrichGO(genes,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1,
                        qvalueCutoff = 1, pAdjustMethod = "BH")
        all = list(mf = mf, cc = cc, bp = bp, kg = kg)
        all[["summary"]] = summarize_cp(all, comparison)
        return(all)
    }
    
    if(type=="under"){
        res.under <- res %>% filter(logFC < -lfc & adj.P.Val < pval)
        genes = res.under$entrezgene
        mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        cc = enrichGO(genes,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        bp = enrichGO(genes,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH",
                      qvalueCutoff = 1, pvalueCutoff = 1)
        kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1,
                        qvalueCutoff = 1, pAdjustMethod = "BH")
        all = list(mf = mf, cc = cc, bp = bp, kg = kg)
        all[["summary"]] = summarize_cp(all, comparison)
        return(all)
    }
}

convert_enriched_ids = function(res, entrezsymbol) {
    res = res %>% mutate(geneID = strsplit(as.character(geneID), "/")) %>% tidyr::unnest(geneID) %>% 
        left_join(entrezsymbol, by = c(geneID = "entrezgene")) %>% group_by(ID, 
                                                                            Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, ont, 
                                                                            comparison) %>% summarise(geneID = paste(geneID, collapse = "/"), symbol = paste(hgnc_symbol, 
                                                                                                                                                             collapse = "/"))
    return(res)
}

load("diff.exp.RData")

res = diff.exp.df
names(res) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "hgnc_symbol")
enrich_rs = enrich_cp(res, paste0(cancertype1, "/", cancertype2), pval, lfc, type="all")

pdf(file = paste0("../", enrichmentfilename, ".pdf"), width = 9, height = 4.5)
dotplot(enrich_rs$kg, x="count", showCategory=10, color="qvalue")
dev.off()

enrich_summary = enrich_rs$summary %>% arrange(p.adjust)
enrich_summary = convert_enriched_ids(enrich_summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write.csv(enrich_summary, paste0("../", enrichmentfilename, ".csv"))
save(enrich_summary, file = "enrich_summary.RData")