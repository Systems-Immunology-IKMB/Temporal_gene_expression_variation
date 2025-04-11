## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(WGCNA)
library(DESeq2)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(reshape2)
library(lme4)
library(ICC)
library(circlize)
library(ComplexHeatmap)
library(rstatix)
library(ppcor)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


##Loading data ------------------------------------------------------------

gene_list <- read.csv("results/filtered_genes_samples25_fpm5.txt", sep = '\t')
gene_list <- gene_list %>%
  subset(.$gene_type == "protein_coding") %>%
  subset(! .$Chr %in% c("chrX", "chrY"))

count_data <- read.csv("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)
col_clinical_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')

col_clinical_data$Annual_season <- ifelse(col_clinical_data$Month %in% c("January", "February", "March"), "Winter", 
                                          ifelse(col_clinical_data$Month %in% c("July", "August", "September"), "Summer", 
                                                 ifelse(col_clinical_data$Month %in% c("April", "May"), "Spring", "Autumn")))

count_data <- count_data[as.character(gene_list$Gene_ID), as.character(col_clinical_data$Sample_ID)]
colnames(count_data) <- col_clinical_data$Description



## Run WGCNA ----------------

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ Individual_ID)
dds_counts <- estimateSizeFactors(dds_counts)
vst_counts <- vst(dds_counts)
vst_counts <- assay(vst_counts)
vst_counts <- t(vst_counts)

sample_tree <- hclust(dist(vst_counts), method = "average")

pdf("results/Coexpression/output/sample_clustering_to_detect_outliers.pdf", width=25)
par(cex =(0.2));
par(mar = c(0,4,2,0))
plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 0.025,cex.axis = 1.5, cex.main = 2)
abline(h = 95, col = "red")
dev.off()

clust <- cutreeStatic(sample_tree, cutHeight = 95, minSize = 10)
table(clust)

#determine parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(vst_counts, powerVector = powers, verbose = 5)
sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.9;

#Scale-free topology fit index as a function of soft-threshold power: Scale Independence Plot
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90, cool="red")

#Mean connectivity Plot
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Make gene tree
adj_mat <- adjacency(vst_counts, power=9)
TOM = TOMsimilarity(adj_mat)
dissTOM = 1-TOM
TOM_gene_tree <- hclust(as.dist(dissTOM), method = "average")

#Define and merge modules
module_labels <- cutreeDynamicTree(dendro=TOM_gene_tree, minModuleSize = 15, deepSplit = TRUE)
module_colors <- labels2colors(module_labels)
mergedColor <- mergeCloseModules(vst_counts, module_colors, cutHeight = .2)$color

#Save co-expression results
gene_module <- data.frame(Gene = colnames(vst_counts), ModuleLabel = module_labels, ModuleColor = mergedColor)
write.table(gene_module, "results/Coexpression/output/gene_module.txt", quote = FALSE, sep = '\t', row.names = FALSE)

MEsO <- moduleEigengenes(vst_counts, mergedColor)$eigengenes
MEs <- orderMEs(MEsO)
rownames(MEs) <- rownames(vst_counts)
write.table(MEs, "results/Coexpression/output/module_eigengenes.txt", quote = FALSE, sep = '\t')

module_size <- as.data.frame(table(gene_module$ModuleColor))
colnames(module_size) <- c("Module_color", "Number_of_genes")
module_size <- module_size[module_size$Module_color != "grey", ]
module_size <- module_size[order(module_size$Number_of_genes, decreasing = TRUE), ]
module_size$Module_name <- paste0("M", 1:nrow(module_size))
write.table(module_size, "results/Coexpression/output/Module_list.txt", sep = '\t', quote = FALSE)




## Compute partial correlation between module eigengenes and clinical parameters --------------------

#prepare clinical data
clinical_traits <- col_clinical_data[, c("Neutrophils...µL.", "Lymphocytes...µL.", "Creatinine..mg.dL.", "Annual_season", "Monocytes...µL.", "Hemoglobin..g.dL.", "Age..years.", "Triglycerides..mg.dL.", "Eosinophils...µL.", "Gender", "Thrombocytes..x1000.µL.", "Time_of_sampling", "Albumin..g.L.", "BMI..kg.m2.", "Uric.acid..mg.dL.", "hsCRP..mg.L.")]
clinical_traits$Gender <- ifelse(clinical_traits$Gender == "Man", 0, 1)
clinical_traits$Annual_season <- ifelse(clinical_traits$Annual_season == "Summer", -1, ifelse(clinical_traits$Annual_season == "Winter", 1, 0))
rownames(clinical_traits) <- col_clinical_data$Description
clinical_traits <- clinical_traits[rownames(vst_counts), ]

#compute spearman partial correlation
eigengenes <- as.data.frame(MEs)
filtered_eigengenes <- eigengenes[, !(colnames(eigengenes)=="MEgrey")]
module_trait_pcor <- pcor(na.omit(cbind(filtered_eigengenes, clinical_traits)), method = "spearman")
module_trait_pcor_matrix <- module_trait_pcor$estimate[colnames(filtered_eigengenes), colnames(clinical_traits)]
module_trait_pcor_pvalue <- module_trait_pcor$p.value[colnames(filtered_eigengenes), colnames(clinical_traits)]
module_trait_pcor_adj_pvalue <- p.adjust(module_trait_pcor_pvalue, method = 'BH')
module_trait_pcor_adj_pvalue <- matrix(module_trait_pcor_adj_pvalue, nrow = nrow(module_trait_pcor_pvalue), ncol = ncol(module_trait_pcor_pvalue))
rownames(module_trait_pcor_matrix) <- gsub("ME", "", rownames(module_trait_pcor_matrix))
rownames(module_trait_pcor_matrix) <- module_size[match(rownames(module_trait_pcor_matrix), module_size$Module_color), "Module_name"]
colnames(module_trait_pcor_matrix) <- c("Neutrophils", "Lymphocytes", "Creatinine", "Annual Season", "Monocytes", "Hemoglobin", "Age", "Triglycerides", "Eosinophils", "Sex", "Thrombocytes", "Time of sampling", "Albumin", "BMI", "Uric acid", "hsCRP")
colnames(module_trait_pcor_adj_pvalue) <- colnames(module_trait_pcor_matrix)
rownames(module_trait_pcor_adj_pvalue) <- rownames(module_trait_pcor_matrix)

module_trait_pcor_matrix <- module_trait_pcor_matrix[as.character(module_size$Module_name), ]
module_trait_pcor_adj_pvalue <- module_trait_pcor_adj_pvalue[as.character(module_size$Module_name), ]

write.table(module_trait_pcor_matrix, "results/Coexpression/protein_coding_genes/output/module_trait_partial_correlation.txt", sep = '\t', quote = FALSE)
write.table(module_trait_pcor_pvalue, "results/Coexpression/output/Clinical_traits/module_trait_partial_correlation_pval.txt", sep = '\t', quote = FALSE)
write.table(module_trait_pcor_adj_pvalue, "results/Coexpression/output/Clinical_traits/module_trait_partial_correlation_adj_pval.txt", sep = '\t', quote = FALSE)





## Plot correlation heatmap as in Figure 5b ------------------------------

#import gene co-expression results
res_module <- read.table("results/Coexpression/output/gene_module.txt", header = TRUE) %>%
  merge(annotation, by.x = "Gene", by.y = "Gene_ID")
df_eigengene <- read.table("results/Coexpression/output/module_eigengenes.txt", header = TRUE)
df_module_list <- read.table("results/Coexpression/output/module_list.txt", header = TRUE)

##compute ICC for clinical variables
var_names <- c("Neutrophils...µL." = "Neutrophils", 
               "Lymphocytes...µL." = "Lymphocytes", 
               "Monocytes...µL." = "Monocytes", 
               "Creatinine..mg.dL." = "Creatinine", 
               "Hemoglobin..g.dL." = "Hemoglobin", 
               "Age..years." = "Age", 
               "Triglycerides..mg.dL." = "Triglycerides", 
               "Eosinophils...µL." = "Eosinophils", 
               "Gender" = "Gender", 
               "Thrombocytes..x1000.µL." = "Thrombocytes", 
               "Albumin..g.L." = "Albumin", 
               "BMI..kg.m2." = "BMI", 
               "hsCRP..mg.L." = "hsCRP", 
               "Uric.acid..mg.dL." = "Uric acid", 
               "Annual_season" = "Annual Season", 
               "Time_of_sampling"  = "Time of sampling")

df_icc <- col_clinical_data[, c("Individual_ID", "Sample_ID", "Neutrophils...µL.", "Lymphocytes...µL.", "Creatinine..mg.dL.", "Monocytes...µL.", "Hemoglobin..g.dL.", "Triglycerides..mg.dL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL.", "Albumin..g.L.", "Uric.acid..mg.dL.", "hsCRP..mg.L.", "Time_of_sampling")] %>%
  melt(id.vars = c("Sample_ID", "Individual_ID")) %>%
  na.omit() %>%
  split(~ .$variable) %>%
  lapply(function(x){res <- ICCest(data = x, Individual_ID, value); return(data.frame("ICC" = res$ICC,
                                                                                      "lower_ICC" = res$LowerCI,
                                                                                      "upper_ICC" = res$UpperCI,
                                                                                      "var_within" = res$varw,
                                                                                      "var_between" = res$vara))}) %>%
  do.call(rbind, .) %>%
  add_column("Variable" = rownames(.), .before = 1)

df_icc$Name <- as_labeller(var_names)(df_icc$Variable) %>% unlist()

write.table(df_icc, "results/Coexpression/output/ICC_clinical_parameters.txt", quote = FALSE, sep = "\t", row.names = FALSE)


##compute ICC for module eigengenes
df_module_list <- df_module_list %>%
  add_column("Variable" = paste0("ME", .$Module_color))

df_icc_modules <- df_eigengene %>%
  add_column("Description" = rownames(.)) %>%
  melt(id.vars = "Description") %>%
  merge(col_clinical_data, by = "Description") %>%
  .[, c("Sample_ID", "Individual_ID", "variable", "value")] %>%
  split(~ .$variable) %>%
  lapply(function(x){res <- ICCest(data = x, Individual_ID, value); return(data.frame("ICC" = res$ICC,
                                                                                      "lower_ICC" = res$LowerCI,
                                                                                      "upper_ICC" = res$UpperCI,
                                                                                      "var_within" = res$varw,
                                                                                      "var_between" = res$vara))}) %>%
  do.call(rbind, .) %>%
  add_column("Variable" = rownames(.), .before = 1) %>%
  merge(df_module_list, by ="Variable")

write.table(df_icc_modules, "results/Coexpression/output/ICC_module_eigengenes.txt", quote = FALSE, sep = "\t", row.names = FALSE)


##plot partial correlation heatmap with ICC for modules and clinical parameters
df_icc_clin_param <- read.table("results/Coexpression/output/ICC_clinical_parameters.txt", sep = "\t", header = TRUE)
df_icc_modules <- read.table("results/Coexpression/output/ICC_module_eigengenes.txt", sep = "\t", header = TRUE)

res_partial_corr <- read.table("results/Coexpression/output/module_trait_partial_correlation.txt", header = TRUE, sep = "\t", check.names = FALSE)
res_adjp <- read.table("Nehas Analysis/Coexpression/output/module_trait_partial_correlation_adj_pval.txt", header = TRUE, sep = "\t", check.names = FALSE)

#convert to binned values to have a discrete scale
my_breaks <- seq(from = -1, to = 1, by = 2/15) %>% round(digits = 2)

my_level <- res_partial_corr$Neutrophils %>% 
  cut(breaks = my_breaks,include.lowest = TRUE, right = TRUE) %>%
  levels()
my_ht_colors <- structure(viridis::viridis(15), names = my_level)

res_partial_corr <- res_partial_corr %>%
  apply(cut, breaks = my_breaks, MARGIN = 1:2, include.lowest = TRUE, right = TRUE)

#add significance to the heatmap
df_adjp <- res_adjp %>%
  t() %>%
  as.data.frame() %>%
  apply(function(x){case_when(x >= 0.05 ~ "",
                              x >= 0.01 ~ "*",
                              x >= 0.0001 ~ "**",
                              .default = "***")}, MARGIN = 1:2)

add_pval <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%s", df_adjp[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
}

#add ICC of clinical parameters and module eigengenes
my_icc_param <- df_icc_clin_param %>%
  add_column("my_color" = ifelse(.$ICC > 0.5, "#009900", "#CCCC00")) %>%
  .[match(res_partial_corr %>% colnames(), .$Name), ] 

row_ha <- rowAnnotation("ICC" = anno_barplot(my_icc_param$ICC, 
                                             gp = gpar(fill = my_icc_param$my_color),
                                             width = unit(1, "cm"),
                                             border = TRUE,
                                             ylim = c(0, 1),
                                             axis_param = list(
                                               side = "bottom",
                                               gp = gpar(fontsize = 5),
                                               at = c(0, 0.5, 1), 
                                               labels = c(" 0.0", "0.5", "1.0"),
                                               labels_rot = 270)),
                        annotation_name_rot = 0,
                        annotation_name_side = "bottom",
                        annotation_name_gp = gpar(fontsize = 9),
                        gp = gpar(col = "black"))

my_icc_module <- df_icc_modules %>%
  add_column("my_color" = ifelse(.$ICC > 0.5, "#009900", "#CCCC00")) %>%
  .[match(res_partial_corr %>% rownames(), .$Module_name), ] 

col_ha <- columnAnnotation("ICC" = anno_barplot(my_icc_module$ICC, 
                                                gp = gpar(fill = my_icc_module$my_color),
                                                border = TRUE,
                                                ylim = c(0, 1),
                                                axis_param = list(
                                                  side = "right",
                                                  gp = gpar(fontsize = 5),
                                                  at = c(0, 0.5, 1), 
                                                  labels = c("0.0", "0.5", "1.0"),
                                                  labels_rot = 270)),
                           annotation_name_rot = 0,
                           annotation_name_side = "right",
                           annotation_name_gp = gpar(fontsize = 9),
                           gp = gpar(col = "black"))

lgd_corplot <- Legend(at = my_level %>% rev(), 
                      title = "Partial\ncorrelation",
                      legend_gp = gpar(fill = viridis::viridis(15) %>% rev()),
                      border =  gpar(col = "black"),
                      labels_gp = gpar(col = "black"),
                      title_gp = gpar(col = "black"))

p <- Heatmap(t(res_partial_corr), 
             name = "Partial\ncorrelation",
             col = my_ht_colors,
             show_heatmap_legend = FALSE,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             show_row_names = TRUE,
             row_names_side = "left",
             show_column_names = TRUE,
             column_names_side = "bottom",
             cell_fun = add_pval,
             top_annotation = col_ha, 
             right_annotation = row_ha,
             row_names_gp = gpar(fontsize = 9),
             column_names_gp = gpar(fontsize = 9),
             row_gap = unit(1, "mm"), 
             border = TRUE, 
             width = nrow(res_partial_corr)*unit(5, "mm"), 
             height = ncol(res_partial_corr)*unit(5, "mm"),
             row_title_gp = gpar(fontsize = 14),
             row_title_rot = 270,
             row_title_side = "right",
             cluster_row_slices = FALSE,
             column_title = NULL,
             show_column_dend = TRUE,
             cluster_column_slices = FALSE)


