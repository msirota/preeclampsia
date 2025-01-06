#setwd('') #to run this code set the directory in which this script is places as working directory

#import libraries and functions
library(MetaIntegrator)
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

out <- runMetaAnalysis(obj)


deconv <- immunoStatesMeta(out)
deconv_meta <- runMetaAnalysis(deconv) 
res_deconv <- deconv_meta$metaAnalysis$pooledResults

saveRDS(deconv_meta, file = 'metanalyses_results/metanalysis_results_CellProp.rds')
