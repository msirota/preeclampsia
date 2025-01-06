#setwd('') #to run this code set the directory in which this script is places as working directory
source("./function_drug_rep.R")
#source("../functions_to_process_RNAseq_Data.R")

library(AnnotationDbi)
library(ggplot2)
library(qvalue)
library(parallel)
ncores = 2 #set number of cores

#Load in cmap drug profiles
load('~/Desktop/SirotaLab/PSPG_245B_drug_repurposing_workshop_2021/data/cmap_signatures.RData')
gene_list <- subset(cmap_signatures,select=1) #split first column (containing gene ids) from the rest of the dataset
cmap_signatures <- cmap_signatures[,2:ncol(cmap_signatures)]  #nb cmap signatures are ranks ??

#obtaining gene signature for preeclampsia 
disease_name <- 'preeclampsia_vs_PretermControls'
file <- "../metanalyses_results/metanalysis_results_GeneExpr.rds"
disease <- readRDS(file)
disease <- disease$metaAnalysis$pooledResults
disease <- disease[disease$effectSizeFDR <0.05,] #select significant adj pvalues
disease <- disease[disease$numStudies >=2,] 
disease <- disease[abs(disease$effectSize) > 1,] 

disease$gene <- rownames(disease)
disease <- disease[, c('gene', 'effectSize')]
colnames(disease) <- c('Gene_name', 'logFC')
disease$Directionality <- ifelse(disease$logFC > 0, 'Upregulated', 'Downregulated')


#map gene symbols to entrez ids
mapping_tbl <- read.csv("~/Desktop/general_scripts_R/gene_id_conversion_table.tsv", sep = '\t') #gene mapping table
mapping_tbl <- na.omit(mapping_tbl[, c('Gene_name', 'entrezID')]) #remove not mapped
mapping_tbl <- mapping_tbl[!duplicated(mapping_tbl), ]
merged <- merge(disease, mapping_tbl, by = 'Gene_name')
colnames(merged) <- c('X', 'logFC', "Directionality", "GeneID")
merged <- merged[order(merged$logFC, decreasing = T),]

dz_signature <- merged #disease signature, colnames = X (gene symbol) logFC Directionality GeneID
dz_signature <- dz_signature[order(dz_signature$logFC),] #sort by fold change in ascending order


####################
# Subset lists of up-regulated genes and down-regulated genes
dz_genes_up <- subset(dz_signature,Directionality=="Upregulated",select="GeneID") #list of upreg genes
dz_genes_down <- subset(dz_signature,Directionality=="Downregulated",select="GeneID") #list of downreg genes


# Intersection of genes in CMap drug profiles and dz_signature
dz_cmap_common_genes <- merge(gene_list, 
                              dz_signature,
                              by.x = "V1", 
                              by.y = "GeneID",
                              all.x = FALSE, 
                              all.y = FALSE)


saveRDS(dz_cmap_common_genes, file = 'dz_cmap_common_genes.rds')

write.table(sum(grepl("Up", dz_cmap_common_genes$Directionality)), file = paste0('upreg_genes_', disease_name,"_tissue.txt"), quote = F, col.names = F, row.names = F)
write.table(sum(grepl("Dow", dz_cmap_common_genes$Directionality)), file = paste0('downreg_genes_', disease_name,"_tissue.txt"), quote = F, col.names = F, row.names = F)

# Calculate distribution of scores using random genes (same number of up- and down-regulated genes)
# against random drugs
N_PERMUTATIONS <- 100000

rand_cmap_scores <- mclapply(sample(1:ncol(cmap_signatures), N_PERMUTATIONS, replace=T), function(exp_id) { 
  cmap_exp_signature <- cbind(gene_list, subset(cmap_signatures, select=exp_id))
  colnames(cmap_exp_signature) <- c("ids","rank")
  random_input_signature_genes <- sample(gene_list[,1], (nrow(dz_genes_up)+nrow(dz_genes_down)))
  rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
  rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
  cmap_score(rand_dz_gene_up, rand_dz_gene_down, cmap_exp_signature) 
}, mc.cores = ncores)

save(rand_cmap_scores,file=paste0("cmap_random_scores_100000_", disease_name, "_tissue.RData"))

#############################
# Compute reversal score for each drug profile in CMap using dz_signature

dz_cmap_scores <- as.numeric(mclapply(1:ncol(cmap_signatures),function(exp_id) {
  cmap_exp_signature <- cbind(gene_list,subset(cmap_signatures,select=exp_id))
  colnames(cmap_exp_signature) <- c("ids","rank")
  cmap_score(dz_genes_up,dz_genes_down,cmap_exp_signature)
}, mc.cores = ncores)) 

# Compute the significance against the random scores
random_scores <- unlist(rand_cmap_scores)
# Frequency-based p-value using absolute scores from sampling distribution to approximate two-tailed p-value
p_values <- sapply(dz_cmap_scores,function(score) {
  length(which(abs(random_scores) >= abs(score))) / length(random_scores)
})

q_values <- qvalue(p_values)$qvalues

subset_comparison_id <- disease_name
analysis_id <- "cmap"

drugs <- data.frame(exp_id = seq(1:length(dz_cmap_scores)), cmap_score = dz_cmap_scores, p = p_values, q = q_values,
                    subset_comparison_id, analysis_id)
results <- list(drugs, dz_signature)
save(results, file = paste0("cmap_predictions_", disease_name, ".RData"))

print('DONE')