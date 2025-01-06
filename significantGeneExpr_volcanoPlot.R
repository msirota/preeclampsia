#setwd('') #to run this code set the directory in which this script is places as working directory

library(ggplot2)
library(ggrepel)

mint_genes <- readLines("mint/mint_signature.txt") 

res <- readRDS("metanalyses_results/metanalysis_results_GeneExpr.rds")
res <- res$metaAnalysis$pooledResults
res <- res[res$numStudies >=2,] #select genes derived from the integration of at least two datasets
res$X <- rownames(res)

res$col_label_text <- ifelse(res$X %in% mint_genes, 'black', 'black')

res$col_points <- ifelse(
  (res$effectSize < 0 & res$effectSizeFDR < 0.05), 'blue', ifelse((res$effectSize > 0 & res$effectSizeFDR < 0.05), 'red', 'black'))


res$X <- unlist(lapply(res$X, function(x) {
  if (x %in% mint_genes) {
    x <-  gsub(x, paste0("bold('", x, "')"), x)
  }
  return(x)
}

))


res$directionality <- ifelse(res$effectSize > 0, 'Upregulated', 'Downregulated')

es_threshold <- 1.5


jpeg(file = 'volcanoPlot_1010geneExpr.jpeg', res = 600, units = 'cm', width = 30, height = 27)
ggplot(res, aes(x = effectSize, y = -log10(effectSizeFDR))) +
  geom_point(aes(col = col_points)) +
  geom_vline(xintercept=c(-1, 1), col="black", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 'dashed') +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 18)) +
  xlim(c(-3, 3)) +
  geom_label_repel(
    data= subset(res, abs(effectSize) > es_threshold & effectSizeFDR < 0.05),
    aes(
      label = subset(res, abs(effectSize) > es_threshold & effectSizeFDR < 0.05)$X,
      color= factor(col_label_text), 
      segment.color="black"
    ),
    size = 5,
    box.padding = 0.3,  
    point.padding = 0.3,
    max.overlaps = Inf,
    parse = T
  ) + scale_color_identity()
dev.off()

