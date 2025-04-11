## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(ggplot2)
library(ggpubr)
library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(tidyverse)
library(ggrepel)
library(ICC)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


## Compute variance of gene expression after randomly selecting one sample per individual ------------------------------------------------------------

#load data
count_data <- read.csv("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)

gene_list <- read.csv("info_files/filtered_genes_samples25_fpm5.txt", sep = '\t')
gene_list <- subset(gene_list, !(gene_list$Chr %in% c("chrX", "chrY")))
rownames(gene_list) <- gene_list$Gene_ID

col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')

#randomly select one sample per individual
set.seed(123)
randomly_selected_col_data <- col_data %>%
  group_by(Individual_ID) %>%
  sample_n(1) %>%
  ungroup()

count_data <- count_data[as.character(gene_list$Gene_ID), as.character(randomly_selected_col_data$Sample_ID)]
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = randomly_selected_col_data, design = ~ 1)
dds_counts <- estimateSizeFactors(dds_counts)
quantvst <- vst(dds_counts)
quantvst <- assay(quantvst)

gene_icc_data <- data.frame()
for (i in 1:nrow(quantvst)) {
  test_data <- data.frame(counts = quantvst[i, ], 
                          individual = randomly_selected_col_data$Individual_ID, 
                          timepoint = randomly_selected_col_data$Timepoint)
  total_var <- var(test_data$counts)
  res <- ICCest(test_data$individual, test_data$counts)
  gene_icc_data <- rbind(gene_icc_data, data.frame(Gene_ID = rownames(quantvst)[i], 
                                                   ICC = res$ICC, 
                                                   CI_lower = res$LowerCI, 
                                                   CI_upper = res$UpperCI, 
                                                   Var_within = res$varw, 
                                                   Var_between = res$vara, 
                                                   Total_var = total_var, 
                                                   Mean_expression = mean(test_data$counts)))
}

gene_icc_data <- left_join(gene_icc_data, gene_list, by = "Gene_ID")
write.table(gene_icc_data, "results/outliers_genes_analysis/output/Variation_with_one_sample_per_individual.txt", sep = '\t', row.names = FALSE, quote = FALSE)



## Dotplot for outlier genes as in Figure 6a ------------------------------

#load ICC results
one_sample_icc <- read.csv("results/outliers_genes_analysis/output/Variation_with_one_sample_per_individual.txt", sep = '\t')
icc_results <- read.csv("results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t')
icc_results <- inner_join(icc_results, one_sample_icc[, c("Gene_ID", "Total_var", "Mean_expression.x")], by = "Gene_ID", suffix = c("", ".onesample"))
icc_results$Var_diff <- icc_results$Total_var.onesample - icc_results$Var_between
icc_results$Var_diff_outlier <- ifelse(icc_results$Var_diff >= 2 * sd(icc_results$Var_between), "Outlier genes", "")

write.table(icc_results, "results/outliers_genes_analysis/Outlier_genes.txt", sep = '\t', quote = FALSE)


label_data <- icc_results[icc_results$Gene_ID %in% rownames(outlier_comb_mat[rowSums(outlier_comb_mat) > 1, ]), ]
label_data <- label_data[order(label_data$Var_diff, decreasing = TRUE), ]
label_data1 <- label_data[1:30,]

p <- ggplot(icc_results, aes(x = Var_between, y = Total_var.onesample, col = Var_diff_outlier)) + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 2 * sd(icc_results$Var_between), lty = 2) +
  geom_abline(slope = 1, intercept = -2 * sd(icc_results$Var_between), lty = 2) +
  coord_fixed(ratio = 1) +
  xlab("Inter-individual variance\n(Multiple samples per individual)") + 
  ylab("Inter-individual variance\n(One sample per individual)") +
  geom_text_repel(data = label_data1, aes(label = gene_name), color = "black", max.overlaps = 100, size = 3) +
  scale_color_manual(values = c("darkgray", "skyblue2"), guide = 'none') +
  theme_pubr() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12)) 




## Barplot for overlap of outlier genes with external DEGs as in Figure 7b ------------------------------------------------------------

#load ICC results
one_sample_icc <- read.csv("results/outliers_genes_analysis/output/Variation_with_one_sample_per_individual.txt", sep = '\t')
icc_results <- read.csv("results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t')
icc_results <- inner_join(icc_results, one_sample_icc[, c("Gene_ID", "Total_var", "Mean_expression.x")], by = "Gene_ID", suffix = c("", ".onesample"))
icc_results$Var_diff <- icc_results$Total_var.onesample - icc_results$Var_between
icc_results$Var_diff_outlier <- ifelse(icc_results$Var_diff >= 2 * sd(icc_results$Var_between), "Outlier genes", "")

#load DEGs from publications on various diseases
disease_gene_list <- read.csv("results/outliers_genes_analysis/output/Disease_genes_list.txt", sep = '\t')
disease_gene_list <- subset(disease_gene_list, disease_gene_list$Gene_ID %in% icc_results$Gene_ID)
disease_gene_list$Disease_category <- ifelse(disease_gene_list$Disease %in% c("Alzheimer's disease", "Huntington's disease", "Parkinson's disease"), "Neurodegenrative diseases", 
                                             ifelse(disease_gene_list$Disease %in% c("COVID-19", "Sepsis"), "Infectious diseases", 
                                                    ifelse(disease_gene_list$Disease == "Coronary Artery disease", "Cardiovascular diseases", "Inflammatory diseases")))

write.table(disease_gene_list, "results/outliers_genes_analysis/Disease_genes.txt", sep = '\t', quote = FALSE)

#prepare table for plot
input_list <- split(disease_gene_list$Gene_ID, disease_gene_list$Disease)
input_list$'Outlier genes' <- icc_results[icc_results$Var_diff_outlier == "Outlier genes", "Gene_ID"]

comb_mat <- as.data.frame(list_to_matrix(input_list))
outlier_comb_mat <- subset(comb_mat, comb_mat$`Outlier genes` == 1)
individual_disease_outliers <- colSums(outlier_comb_mat) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  `colnames<-`(c("Disease", "Intersections")) 
individual_disease_outliers$Group <- "Disease"
individual_disease_outliers$Total_genes <- colSums(comb_mat)

category_input_list <- split(disease_gene_list$Gene_ID, disease_gene_list$Disease_category)
category_input_list$'Outlier genes' <- icc_results[icc_results$Var_diff_outlier == "Outlier genes", "Gene_ID"]

category_comb_mat <- as.data.frame(list_to_matrix(category_input_list))
category_outlier_comb_mat <- subset(category_comb_mat, category_comb_mat$`Outlier genes` == 1)
disease_category_outliers <- colSums(category_outlier_comb_mat) %>% as.data.frame() %>% rownames_to_column() %>% `colnames<-`(c("Disease", "Intersections")) 
disease_category_outliers$Group <- "Category"
disease_category_outliers$Total_genes <- colSums(category_comb_mat)

disease_category_overlap <- data.frame(Disease = c("Cardiovascular diseases", "Infectious diseases", "Inflammatory diseases", "Neurodegenrative diseases"),
                                       Intersections = c(nrow(category_outlier_comb_mat[category_outlier_comb_mat$`Cardiovascular diseases` == 1 & rowSums(category_outlier_comb_mat) > 2, ]), 
                                                         nrow(category_outlier_comb_mat[category_outlier_comb_mat$`Infectious diseases` == 1 & rowSums(category_outlier_comb_mat) > 2, ]), 
                                                         nrow(category_outlier_comb_mat[category_outlier_comb_mat$`Inflammatory diseases` == 1 & rowSums(category_outlier_comb_mat) > 2, ]), 
                                                         nrow(category_outlier_comb_mat[category_outlier_comb_mat$`Neurodegenrative diseases` == 1 & rowSums(category_outlier_comb_mat) > 2, ])), 
                                       Group = "Category", 
                                       Total_genes = NA, 
                                       Category = c("Cardiovascular diseases1", "Infectious diseases1", "Inflammatory diseases1", "Neurodegenrative diseases1"), 
                                       width = 1, 
                                       label = NA)

plot_data <- rbind(individual_disease_outliers, disease_category_outliers)
plot_data$Category <- ifelse(plot_data$Group == "Disease", disease_gene_list[match(plot_data$Disease, disease_gene_list$Disease), "Disease_category"], plot_data$Disease) 
plot_data <- plot_data[!plot_data$Disease == "Outlier genes", ]
plot_data$width <- ifelse(plot_data$Group == "Disease", 0.9, 1)
plot_data$Disease <- factor(plot_data$Disease, levels = c("Coronary Artery disease", "Cardiovascular diseases", "Alzheimer's disease", "Huntington's disease", "Parkinson's disease", "Neurodegenrative diseases", "Systemic sclerosis", "Inflammatory bowel disease", "Juvenile idiopathic arthritis", "Inflammatory diseases", "Sepsis", "COVID-19", "Infectious diseases"))
plot_data$label <- paste0(plot_data$Intersections, "/", plot_data$Total_genes)


plot_data1 <- rbind(plot_data, disease_category_overlap)
plot_data1[10:13, "Intersections"] <- plot_data1[10:13, "Intersections"] - disease_category_overlap$Intersections

p1 <- ggplot(plot_data1, aes(x = Disease, y = Intersections, fill = Category, label = label)) + 
  geom_bar(aes(alpha = width), stat = "identity", col = "black") +
  geom_text(nudge_y = 20, angle = 60, vjust = 1, size = 2.5) +
  scale_alpha(range = c(0.7,1), guide = 'none') +
  scale_fill_manual(values = c("#D8A47F", "#B8E1FF","#272932", "#B8E1FF","#0F7173","#B8E1FF", "#F05D5E", "#B8E1FF"), name = "") +
  theme_pubr() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.5, "line")) 




