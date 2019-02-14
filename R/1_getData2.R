#! /usr/bin/env Rscript

## Get the RNASeq and the clinical data

options(stringsAsFactors = FALSE)

setwd("R")
source("0_loadLibraries.R")

loadpkg("readxl")

tbl <- read.csv("./raw_counts.tsv", 
                header = T, 
                row.names = 1, 
                sep = "\t", 
                stringsAsFactors = F,
                check.names = F)

tbl <- tbl[!startsWith(rownames(tbl), '?'),]

new.names <- sapply(rownames(tbl), function(r) {
    strsplit(r, split = '|', fixed = T)[[1]][1]
})

are.duplicated <- duplicated(new.names)

brca_rnaseq <- tbl[!are.duplicated,]

rownames(brca_rnaseq) <- new.names[!are.duplicated]

## Download Supplementary data and import Table 1
## [Supplementary table 1 from Comprehensive molecular portraits of human breast tumours](http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html?foxtrotcallback=true#supplementary-information) has TCGA breast cancer samples classified by PAM50 subtypes and ER Status, PR Status and HER2 Final Status.

temp <- tempfile()
download.file("https://media.nature.com/original/nature-assets/nature/journal/v490/n7418/extref/nature11412-s2.zip",temp)
unzip(temp)
sample_data <- read_excel("nature11412-s2/Supplementary Tables 1-4.xls", sheet = 1, skip = 1)
rm(temp)

save(sample_data, file="sample_data.RData")
save(brca_rnaseq, file="brca_rnaseq.RData")
