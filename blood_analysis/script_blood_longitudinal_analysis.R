setwd('~/Desktop/github_repos/preeclampsia/blood_analysis')
library(gridExtra)
library(labelled)   
library(rstatix)    
library(ggpubr)     
library(GGally)     
library(car)       
library(Epi)        
library(lme4)       # linear mixed-effects models
library(lmerTest)   
library(emmeans)    
library(multcomp)   
library(geepack)    
library(ggeffects)  
library(ggnewscale)
library(dplyr)
library(tidyverse)  
theme_set(theme_minimal() + theme(legend.position = "bottom")) 

mint_genes <- readLines("../mint/mint_signature.txt")

data <- readRDS("metaobj_GSE149437.rds")
expr <- as.data.frame(data$expr)

write.table(data$pheno, file = 'blood_metadata.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

wide_to_long <- function(dataframe, name_row_column, name_column_column, name_value_column) {
  dataframe <- as.data.frame(dataframe)
  out_long <- dataframe %>%
    rownames_to_column(var = name_row_column)
  out_long <- out_long %>%
    pivot_longer(cols = -name_row_column, 
                 names_to = name_column_column, 
                 values_to = name_value_column)
  
  return(out_long)
}

expr_long <- wide_to_long(expr, name_row_column = 'genes', name_column_column = 'samples' ,name_value_column = 'values') 


meta <- data$pheno
meta$ga_blood_draw <- as.numeric(meta$`gestational age:ch1`)
meta$class <- ifelse(grepl("Ear", meta$`group:ch1`), "case", "control")
meta$class <- factor(meta$class, levels = c('control', 'case'))
meta$patientid <- meta$`individual:ch1`

#identical(colnames(expr), meta$geo_accession)

all_pvals <- c()
all_plots <- list()
gene_names <- c()
c <- 0
for (gene in mint_genes) {
  
  subexpr <- expr[rownames(expr) == gene, , drop = F]
  meta <- meta[, c('geo_accession', 'ga_blood_draw', 'class', 'patientid')]
  
  if (!identical(colnames(subexpr), meta$geo_accession)) {
    stop("Non identical samples")
  }
  
  subexpr <- as.data.frame(t(subexpr))
  colnames(subexpr) <- 'expr_value'
  subexpr$geo_accession <- rownames(subexpr)
  merged_final <- merge(meta, subexpr, by = 'geo_accession')
  

  model <- lmer(expr_value ~ class * poly(ga_blood_draw, degree = 2) + (1 | patientid), data = merged_final)
  predictions <- ggpredict(model, terms = c("ga_blood_draw [all]", "class"))
  
  pval <- coef(summary(model))["classcase:poly(ga_blood_draw, degree = 2)1", "Pr(>|t|)"]
  all_pvals <- append(all_pvals, pval)
  gene_names <- append(gene_names, gene)
  
  if (pval < 0.05) {
    c <- c+1
    p2 <-  plot(predictions, show_data = T, dot_size = 3, line_size = 2) +
      theme_classic() +
      ggtitle(paste0(gene, " | pval: ", round(pval, 4))) +
      theme(axis.title = element_text(size = 20),
            plot.title = element_text(size = 27),
            legend.position = 'none') +
      scale_color_manual(values = c("#308ed3", '#f88933')) +
      scale_fill_manual(values = c("#308ed3", '#f88933'))
    
    all_plots[[c]] <- p2
    
  }

}


jpeg(file = 'longitudinal_plots_for_mint_gene_signature.jpeg', res = 300, width = 50, height = 80, units = 'cm')
do.call(grid.arrange,c(all_plots, ncol = 2))
dev.off()
