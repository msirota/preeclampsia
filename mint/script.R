#setwd('') #to run this code set the directory in which this script is places as working directory
library(mixOmics)
library(ComplexHeatmap)
library(circlize)

#load datasets
tissue_gses<-c("GSE74341", "GSE75010", "GSE25906",
               "GSE114691", "GSE218039", "GSE54618") 

files <- list.files("../datasets", recursive = T, pattern = '^metaobj.*rds', full.names = T)

meta_object<-list()
for (gse in tissue_gses){ #for each study
  file_expr_matrix<- files[grepl(gse, files)]
  expr_matrix<-readRDS(file_expr_matrix)$expr
  print(nrow(expr_matrix))
  meta_object$originalData[[gse]]$formattedName <- gse
  meta_object$originalData[[gse]]$expr <- expr_matrix
}
d<-meta_object #objects that contains all the gene expr matrices


common_genes<-c() #intersection of gene in all studies
for (i in seq(1, length(d$originalData))){
  if (i==1){
    common_genes<-rownames(d$originalData[[i]]$expr)
  }else{
    common_genes<-intersect(common_genes,rownames(d$originalData[[i]]$expr))
  }
}

for (i in seq(1, length(d$originalData))){ 
  d$originalData[[i]]$expr <- d$originalData[[i]]$expr[rownames(d$originalData[[i]]$expr) %in% common_genes,]
  d$originalData[[i]]$expr <- d$originalData[[i]]$expr[ order(row.names(d$originalData[[i]]$expr)), ]
}


data<-data.frame() #merge all studies in just one matrix, gene x samples
for (i in seq(1, length(d$originalData))){
  if (i==1){
    data<-d$originalData[[i]]$expr
  }else{
    data<-cbind(data, d$originalData[[i]]$expr)
  }
}

data<-t(data) #samples x genes


dataset_belonging<-c() #for each sample, vector that provides belonging dataset
for (i in seq(1, length(d$originalData))){
  dataset_belonging<-append(dataset_belonging, rep(d$originalData[[i]]$formattedName,ncol(d$originalData[[i]]$expr) ))
}
dataset_belonging<-as.factor(dataset_belonging)
summary(dataset_belonging)

gses<-c() #list of studies
for (i in seq(1, length(d$originalData))){
  gses<-append(gses, d$originalData[[i]]$formattedName)
  
}

obj_GSE25906 <- readRDS("../datasets/metaobj_GSE25906.rds")
pheno1 <- obj_GSE25906$pheno
class1 <- obj_GSE25906$class
pheno1 <- pheno1[pheno1$geo_accession %in% names(class1), ]
if ( !identical(pheno1$geo_accession, names(class1)) ) {
  stop("samples don't correspond")
}

obj_GSE74341 <- readRDS("../datasets/metaobj_GSE74341.rds")
pheno2 <- obj_GSE74341$pheno
class2 <- obj_GSE74341$class
pheno2 <- pheno2[pheno2$geo_accession %in% names(class2), ]
if ( !identical(pheno2$geo_accession, names(class2)) ) {
  stop("samples don't correspond")
}


obj_GSE75010 <- readRDS("../datasets/metaobj_GSE75010.rds")
pheno3 <- obj_GSE75010$pheno
class3 <- obj_GSE75010$class
pheno3 <- pheno3[pheno3$geo_accession %in% names(class3), ]
if ( !identical(pheno3$geo_accession, names(class3)) ) {
  stop("samples don't correspond")
}

obj_GSE114691 <- readRDS("../datasets/metaobj_GSE114691.rds")
pheno4 <- obj_GSE114691$pheno
class4 <- obj_GSE114691$class
pheno4 <- pheno4[pheno4$geo_accession %in% names(class4), ]
if ( !identical(pheno4$geo_accession, names(class4)) ) {
  stop("samples don't correspond")
}

obj_GSE218039 <- readRDS("../datasets/metaobj_GSE218039.rds")
pheno5 <- obj_GSE218039$pheno
class5 <- obj_GSE218039$class
pheno5 <- pheno5[pheno5$geo_accession %in% names(class5), ]
if ( !identical(pheno5$geo_accession, names(class5)) ) {
  stop("samples don't correspond")
}

obj_GSE54618 <- readRDS("../datasets/metaobj_GSE54618.rds")
pheno6 <- obj_GSE54618$pheno
class6 <- obj_GSE54618$class
pheno6 <- pheno6[pheno6$geo_accession %in% names(class6), ]
if ( !identical(pheno6$geo_accession, names(class6)) ) {
  stop("samples don't correspond")
}


pheno1$series_id <- rep('GSE25906', nrow(pheno1))
pheno2$series_id <- rep('GSE74341', nrow(pheno2))
pheno3$series_id <- rep('GSE75010', nrow(pheno3))
pheno4$series_id <- rep('GSE114691', nrow(pheno4))
pheno5$series_id <- rep('GSE218039', nrow(pheno5))
pheno6$series_id <- rep('GSE54618', nrow(pheno6))

pheno1$class <- class1
pheno2$class <- class2
pheno3$class <- class3
pheno4$class <- class4
pheno5$class <- class5
pheno6$class <- class6



pheno1 <- pheno1[, c('geo_accession', 'class', 'series_id')]
pheno2 <- pheno2[, c('geo_accession', 'class', 'series_id')]
pheno3 <- pheno3[, c('geo_accession', 'class', 'series_id')]
pheno4 <- pheno4[, c('geo_accession', 'class', 'series_id')]
pheno5 <- pheno5[, c('geo_accession', 'class', 'series_id')]
pheno6 <- pheno6[, c('geo_accession', 'class', 'series_id')]


metadata <- rbind(pheno1, pheno2, pheno3, pheno4, pheno5, pheno6) #create all metadata matrix

metadata$disease <- ifelse(grepl("0", metadata$class, ignore.case = T), 'Control', 'Preeclampsia')
metadata$class <- NULL

metadata<-metadata[match(rownames(data), metadata$geo_accession),]

disease_state<-metadata$disease
disease_state<-as.factor(disease_state)

# summary(disease_state)
# table(disease_state, dataset_belonging)


tune.mint<- tune(
  method = 'mint.splsda',
  X = data,
  Y = disease_state,
  ncomp = 3,
  study= dataset_belonging,
  test.keepX = seq(1, 100, 1),
  test.keepY = NULL,
  nrepeat = 10,
  validation = "Mfold",
  folds = 10,
  dist = "max.dist",
  auc = TRUE,
  progressBar = TRUE,
  center = TRUE,
  scale = TRUE,
  max.iter = 100,
  tol = 1e-09,
  light.output = TRUE)


saveRDS(tune.mint, file='tuned_mint_preeclampsia_vs_preterm.rds')
#or
#tune.mint<-readRDS('tuned_mint_preeclampsia_vs_preterm.rds')


mint.splsda.res = mint.splsda(X = data, Y = disease_state, study = dataset_belonging,
                              ncomp = 2,keepX = tune.mint$choice.keepX)


jpeg(file = 'res.splsda_dim_reduction_single_studies.jpeg', res = 600, width = 5000, height = 5000)
plotIndiv(mint.splsda.res, 
          study = 'all.partial', 
          title = 'MINT sPLS-DA | Preeclampsia single studies',
          subtitle = c("GSE114691", "GSE218039", "GSE25906",
                       "GSE54618", "GSE74341", "GSE75010"),legend = T,
          X.label = "X-variate 1", 
          Y.label = "X-variate 2", 
          size.legend = 18, 
          legend.title = 'Condition', 
          size.legend.title = 18
)
dev.off()


jpeg(file = 'res.splsda_dim_reduction_all_studies.jpeg', res = 600, width = 6000, height = 5000)
plotIndiv(mint.splsda.res, 
          study = 'global', 
          legend = TRUE, 
          title = 'MINT sPLS-DA', 
          subtitle = 'All studies', 
          ellipse=T,
          X.label = "X-variate 1", 
          Y.label = "X-variate 2", 
          size.legend = 18, 
          legend.title = 'Condition', 
          size.legend.title = 18)
dev.off()


ht<- cim(mint.splsda.res, 
         comp = 1, 
         row.sideColors = color.mixo(as.numeric(disease_state)), 
         row.names = T, title = "MINT sPLS-DA, first component", save='jpeg', name.save='heatmap_preeclampsia_30062024', scale = T)


mat <- ht$mat
mat <- mat[sort(ht$row.names), sort(ht$col.names)]

metadata <- metadata[rownames(mat), ]
#identical(metadata$geo_accession, rownames(mat))


palette_for_anno <- list(
  Condition = c("Control" = '#388ecd', 'Preeclampsia' = '#f78b33'),
  Study = c("GSE75010" = '#ef476f', 'GSE114691' = '#26547c', 'GSE74341' ='#ffd166' , 
              'GSE25906' = '#06d6a0', "GSE218039" = 'purple', 'GSE54618' = 'blue')
  
)

anno <- rowAnnotation(
  Condition = metadata$disease,
  Study = metadata$series_id,
  col = palette_for_anno,
  annotation_name_gp = gpar(fontsize = 18),
  annotation_legend_param = list(title_gp = gpar(fontsize = 16))
)


col_fun<-colorRamp2(c(-4.45,0,4.45), c("#2CB8A8","#fffec0","#EB3100"))

h <- Heatmap(mat, 
             right_annotation = anno,
             col = col_fun,
             #cluster_rows = ht$ddr,
             cluster_columns = ht$ddc,
             row_order = ht$row.names,
             column_names_gp = gpar(fontsize = 20),
             heatmap_width = unit(1, "cm")*ncol(mat),
             heatmap_height = unit(0.1, "cm")*nrow(mat),
             heatmap_legend_param = list(title = 'GeneExpr', title_gp = gpar(fontsize=16)),
             show_row_names = F, column_names_rot = 45
)

h <- draw(h, heatmap_legend_side = "top", annotation_legend_side = "top")


write.table(colnames(mat), file = 'mint_signature.txt',  sep = '\n', row.names = F, col.names = F, quote = F)

jpeg(file = 'mint_signature_heatmap.jpeg', res = 600, units = 'cm', width = 30, height = 30)
h
dev.off()
