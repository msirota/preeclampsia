#setwd('') #to run this code set the directory in which this script is places as working directory

library(ggplot2)
library(ggrepel)

res <- readRDS("metanalyses_results/metanalysis_results_Pathways.rds")
res <- res$metaAnalysis$pooledResults
res <- res[res$numStudies >=2,] #select pathways derived from the integration of at least two datasets
res <- res[res$effectSizeFDR < 0.05, ] #select pathways with fdr < 0.05
res <- res[abs(res$effectSize) > 1,] #select pathways with abs(effect size) > 1
res$X <- rownames(res)
res$X <- gsub("WP_|KEGG_|HALLMARK_","",  res$X)

jpeg(file = 'barPlot_43pathways.jpeg', res = 600, units = 'cm', width = 38, height = 26)
ggplot(res) +
  geom_bar(
    aes(x = effectSize, y = reorder(X, effectSize)) , fill = ifelse(res$effectSize > 0, '#E6B111', 'darkgreen'),
    stat = 'identity') +
  theme_classic() +
  ylab('') +
  theme(axis.text.y = element_text(face="bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16)
        )

dev.off()
