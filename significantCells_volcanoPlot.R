#setwd('') #to run this code set the directory in which this script is places as working directory

library(ggplot2)
library(ggrepel)

redundant_cells <- readLines("redundant_cell_types.txt")

res <- readRDS("metanalyses_results/metanalysis_results_CellProp.rds")
res <- res$metaAnalysis$pooledResults
res$X <- rownames(res)
res <- res[!res$X %in% redundant_cells,]
res <- res[res$numStudies >=2,] #select cells derived from the integration of at least two datasets

es_threshold <- 0.5
col_lab <- ifelse((abs(res$effectSize) > es_threshold & res$effectSizeFDR < 0.05), 'Significant', 'Non-significant')

res$col_label <- col_lab

res$X <- c("Hematopoietic progenitors", "M2 macrophages", "MAST cells", "Naive B cells", "Memory B cells", "Neutrophils", "CD14+ monocytes", "M0 macrophages",
           "CD56bright NK cells", "CD8+ T cells", "CD56dim NK cells", "CD16+ monocytes", "Eosinophils", "Plasma cells",
           "Basophils", "CD4+ T cells", "M1 macrophages", "Plasmacytoid dendritic cells", "gamma-delta T cells",  "Myeloid dendritic cells") #formatting names

jpeg(file = 'volcanoPlot_6Cells.jpeg', res = 600, units = 'cm', width = 30, height = 26)
ggplot(res, aes(x = effectSize, y = -log10(effectSizeFDR))) +
  geom_point(aes(col=col_label), size =5) +
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 'dashed') +
  scale_color_manual(values=c("black", "red"))+
  theme_classic() +
  theme(axis.title = element_text(size = 19))+
  xlim(c(-3, 3)) +
  geom_label_repel(data=subset(res, abs(effectSize) > es_threshold & effectSizeFDR < 0.05),
                   aes(
                     label=subset(res, abs(effectSize) > es_threshold & effectSizeFDR < 0.05)$X
                     ),
                   max.overlaps = Inf,
                   size = 7,
                   box.padding = 0.3,  # Padding around the label
                   point.padding = 0.3,  # Padding between points and text
                   segment.color = 'grey50'
                  )+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 17)   # Increase font size of legend labels
  )

dev.off()