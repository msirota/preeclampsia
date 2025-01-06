#setwd('') #to run this code set the directory in which this script is places as working directory

#import libraries and functions
library(MetaIntegrator)
library(GSVA)
source("./functions_to_process_RNAseq_Data.R")
source("./function.R")

#import datasets
d1 <- readRDS("datasets/metaobj_GSE25906.rds")
d2 <- readRDS("datasets/metaobj_GSE74341.rds")
d3 <- readRDS("datasets/metaobj_GSE75010.rds")

d4 <- readRDS("datasets/metaobj_GSE54618.rds")
d5 <- readRDS("datasets/metaobj_GSE114691.rds")
d6 <- readRDS("datasets/metaobj_GSE218039.rds")


obj <- list()
obj$originalData <- list()

obj$originalData$GSE25906 <- d1
obj$originalData$GSE74341 <- d2
obj$originalData$GSE75010 <- d3
obj$originalData$GSE54618 <- d4
obj$originalData$GSE114691 <- d5
obj$originalData$GSE218039 <- d6

if(!checkDataObject(obj, "Meta", "Pre-Analysis")) {
  stop()
}

#import genesets
hallmark_gs <- c(readRDS("genesets/hallmark_gs.RData"),
                 readRDS("genesets/kegg_gs.RData"),
                 readRDS("genesets/wikipathways_gs.RData")
                 )

for (gse in names(obj$originalData)) {
  expr <- obj$originalData[[gse]]$expr
  pas_mat <- gsva(expr = as.matrix(expr), gset.idx.list = hallmark_gs, method = 'ssgsea', ssgsea.norm = F, min.sz=5, parallel.sz=8)
  pas_mat <- scale(pas_mat)
  obj$originalData[[gse]]$expr <- pas_mat
  
  
  obj$originalData[[gse]]$keys <- rownames(pas_mat)
  names(obj$originalData[[gse]]$keys)<- rownames(pas_mat)
}

out <- runMetaAnalysis(obj)
res <- as.data.frame(out$metaAnalysis$pooledResults)


saveRDS(out, file = 'metanalyses_results/metanalysis_results_Pathways.rds')
