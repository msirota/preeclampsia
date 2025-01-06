#setwd('') #to run this code set the directory in which this script is places as working directory
#source(../functions_to_process_RNAseq_Data.R')
library(pheatmap)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(ggtext)

n_upreg_genes <- as.numeric(readLines("upreg_genes_preeclampsia_vs_PretermControls_tissue.txt"))
n_downreg_genes <- as.numeric(readLines("downreg_genes_preeclampsia_vs_PretermControls_tissue.txt"))
tot_dis_signature <- n_upreg_genes + n_downreg_genes


#load cmap_signatures
load('./cmap_signatures.RData')   

disease_name <- 'preeclampsia_vs_PreTermControls'
#load pipeline results
load("cmap_predictions_preeclampsia_vs_PretermControls.RData")  # here you still don't have the drugs names

cmap_experiments <- read.csv("./cmap_drug_experiments_new.csv", stringsAsFactors =  F) #cmap profiles metadata
valid_instances <- read.csv("./cmap_valid_instances.csv", stringsAsFactors = F)

drug_preds <- results[[1]]
dz_sig <- results[[2]]

cmap_experiments_valid <- merge(cmap_experiments, valid_instances, by="id")
#keep valid profiles; drugs (as long as matched to Drugbank)
cmap_experiments_valid <- cmap_experiments_valid[cmap_experiments_valid$valid == 1 & cmap_experiments_valid$DrugBank.ID != "NULL", ]


drug_instances_all <- merge(drug_preds, cmap_experiments_valid, by.x="exp_id", by.y="id") #map drug to cmap score
write.csv(drug_instances_all, file =paste0("cmap_results_filtered_", disease_name, "_tissue.csv"))
write.csv(drug_instances_all, file =paste0("cmap_results_", disease_name, "_tissue.csv"))


#plot scores
drug_instances_plot <- subset(drug_instances_all, q < 0.05)
drug_instances_plot <- drug_instances_plot %>% 
  group_by(name) %>% 
  dplyr::slice(which.min(cmap_score))

drug_instances_plot <- drug_instances_plot[order(drug_instances_plot$cmap_score), ]
drug_instances_plot$name <- factor(drug_instances_plot$name, levels = drug_instances_plot$name)
drug_instances_plot <- drug_instances_plot %>% 
  mutate(mycolor = ifelse(cmap_score < 0, "red", "darkgreen"))

write.table(drug_instances_plot, file = paste0('drugs_sign_', disease_name, '_tissue.tsv'), col.names = T, row.names = F, quote = F, sep = '\t')

#plot all significant drugs
ggplot(drug_instances_plot, aes(x= name, y = cmap_score)) +
  geom_segment(aes(x=name, xend=name, y=0, yend=cmap_score,color=mycolor))+
  scale_color_identity() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#plot all significant drugs with camp_score < 0
drug_instances_plot_neg <- drug_instances_plot[drug_instances_plot$cmap_score<0,]


#create supp table 4
# supp_table4 <- drug_instances_plot_neg[, c('name', 'cmap_score', 'q', 'concentration', 'duration',
#                                            'cell_line', 'DrugBank.ID')]
# supp_table4 <- rename_single_column(supp_table4, 'name', 'Drug name')
# supp_table4 <- rename_single_column(supp_table4, 'q', 'q-value')
# supp_table4 <- rename_single_column(supp_table4, 'DrugBank.ID', 'DrugBank ID')
# write.table(supp_table4, file = 'suppTable4.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

plot <- ggplot(drug_instances_plot_neg, aes(x=reorder(name, -cmap_score), y=cmap_score)) +
  geom_point(size = 0) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cmap_score), linewidth=1) +
  coord_flip() +
  scale_y_reverse() +
  theme_classic() +
  theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_markdown(size =15),
        axis.title.x = element_text(size=14)
        ) +
  xlab('') +
  scale_x_discrete(labels = c(
    "lansoprazole" = "<span style='color:red; font-size:20px;'><strong>lansoprazole</strong></strong></span>")
    )

jpeg(filename = 'line_plot_list_drugs.jpeg', res = 800, width =9 , height =12 , units = 'in')
plot
dev.off()

write.table(drug_instances_plot_neg$name, file = 'drug_names.txt', sep = '\n', quote = F, col.names = F, row.names = F)
saveRDS(drug_instances_plot_neg, "drug_instances_plot_neg.rds")


###Plot heatmap


#Apply thresholds for significant hits: here, we apply FDR < 0.05 and keep only the reversed profiles (cmap_score < 0)
drug_instances <- subset(drug_instances_all, q < 0.05 & cmap_score < 0)

#Since drugs in cmap have been tested multiple times, we keep the most negative score to be as inclusive as possible
#Alternatively, you could aggregate the scores in a different manner (e.g. averages, or based on metadata)
drug_instances <- drug_instances %>% 
  group_by(name) %>% 
  dplyr::slice(which.min(cmap_score))

drug_instances <- drug_instances[order(drug_instances$cmap_score), ]
drug_instances_id <- c(drug_instances$exp_id) + 1 
#get candidate drugs
drug_signatures <- cmap_signatures[,c(1, drug_instances_id)] 

drug_dz_signature <- merge(dz_sig[, c("GeneID", "logFC")], drug_signatures, by.x = "GeneID", by.y="V1") 
colnames(drug_dz_signature)[2] <- "value"
drug_dz_signature <- drug_dz_signature[order(drug_dz_signature$value),] 

#Convert disease and drug values to ranks from 1:numgenes
drug_dz_signature[,2] <- -drug_dz_signature[,2] # higher rank corresponds to more overexpressed, so we need to reverse order of disease sig


for (i in 2:ncol(drug_dz_signature)){ #for each column turn fold change value into rank
  drug_dz_signature[,i] <- rank(drug_dz_signature[,i])
}
drug_dz_signature <- drug_dz_signature[order(drug_dz_signature[,2]),] #order by disease expression

gene_ids <- drug_dz_signature[,1] 
drug_dz_signature <- drug_dz_signature[, -1]
gene_names <- mget(x=as.character(gene_ids), envir=org.Hs.egSYMBOL)

drug_names <- sapply(2:ncol(drug_dz_signature), function(id){
  #need to minus 1 as in cmap_signatures, V1 is gene id
  new_id <- strtoi(paste(unlist(strsplit(as.character(colnames(drug_dz_signature)[id]),""))[-1], collapse="")) - 1 
  cmap_experiments_valid$name[cmap_experiments_valid$id == new_id]
})
colnames(drug_dz_signature)[-1] <- drug_names

#The disease signature is imbalanced (more downregulated than upregulated genes), so for the following heatmap
#visualization, we will rescale the disease signature rankings so that, for the disease signature, every
#blue square indicates downregulation and every red square indicates upregulation. For the drug profiles,
#we use a straight ranking scale (i.e. even numbers of red and blue).

#702 shared genes, 504 downregulated, 198 upregulated
reassigned_ranks <- c(
  seq(1,702/2,length.out=198), seq(702/2, 702, length.out=504)
) 


drug_dz_signature$value <- reassigned_ranks

write.csv(drug_dz_signature, "drug_dz_signature_all_cmap_hits.csv")

drug_names <- drug_names[drug_names %in% colnames(drug_dz_signature)]

colnames(drug_dz_signature) <- c("Preeclampsia Meta-Analysis", colnames(drug_dz_signature)[2:ncol(drug_dz_signature)])


#FIGURE: all hits from cmap, using a red/blue color scheme
pdf("heatmap_ALL_cmap_hits.pdf", width = 14, height = 6)
layout(matrix(1))
par(mar=c(10, 7, 4, 4))
colPal <- redblue(100)
image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.015),
     labels = c( "Preeclampsia Signature", drug_names), srt = 50, pos=2, offset=-0.2, xpd = TRUE, cex=1, 
     col = "black")
dev.off()


jpeg("heatmap_ALL_cmap_hits.jpeg", width = 8000, height = 4000, res =600)
layout(matrix(1))
par(mar=c(10, 7, 4, 4))
colPal <- redblue(100)
image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.015),
     labels = c( "Preeclampsia Signature", drug_names), srt = 50, pos=2, offset=-0.2, xpd = TRUE, cex=1, 
     col = "black")
dev.off()
