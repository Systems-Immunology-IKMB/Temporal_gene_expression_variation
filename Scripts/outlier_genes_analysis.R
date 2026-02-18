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



## Dotplot for outlier genes as in Figure 6e ------------------------------

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




## Barplot for overlap of outlier genes with external DEGs as in Figure 6f ------------------------------------------------------------

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

input_list <- split(disease_gene_list$Gene_ID, disease_gene_list$Disease)
input_list$'Outlier genes' <- icc_results[icc_results$Var_diff_outlier == "Outlier genes", "Gene_ID"]

comb_mat <- as.data.frame(list_to_matrix(input_list))
outlier_comb_mat <- subset(comb_mat, comb_mat$`Outlier genes` == 1)
individual_disease_outliers <- colSums(outlier_comb_mat) %>% as.data.frame() %>% rownames_to_column() %>% `colnames<-`(c("Disease", "Intersections")) 
individual_disease_outliers$Group <- "Disease"
individual_disease_outliers$Total_genes <- colSums(comb_mat)

# vector of disease columns (everything except Outlier genes)
disease_cols <- setdiff(colnames(comb_mat), "Outlier genes")

# total universe size
U <- nrow(icc_results)

results <- lapply(disease_cols, function(disease) {
  
  # genes in disease set
  disease_genes <- rownames(comb_mat)[comb_mat[, disease] == 1]
  
  # genes in outlier set
  outlier_genes <- rownames(comb_mat)[comb_mat[, "Outlier genes"] == 1]
  
  # counts
  a <- length(intersect(disease_genes, outlier_genes))  # overlap
  b <- length(outlier_genes) - a                         # outlier only
  c <- length(disease_genes) - a                         # disease only
  d <- U - a - b - c                                     # neither
  
  # Fisher test (enrichment)
  ft <- fisher.test(
    matrix(c(a, b, c, d), nrow = 2),
    alternative = "greater"
  )
  
  data.frame(
    Disease = disease,
    Disease_genes = length(disease_genes),
    Outlier_genes = length(outlier_genes),
    Overlap = a,
    Odds_ratio = unname(ft$estimate),
    P_value = ft$p.value
  )
})

results_df <- do.call(rbind, results)

# Multiple testing correction
results_df$FDR <- p.adjust(results_df$P_value, method = "BH")

results_df

write.table(results_df, "results/Disease_genes_analysis/Overlap_enrichment.txt", sep = '\t', quote = FALSE, row.names = FALSE)

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
plot_data$Disease <- factor(plot_data$Disease, levels = c("Coronary Artery disease", "Cardiovascular diseases", 
                                                          "Alzheimer's disease", "Huntington's disease", "Parkinson's disease", "Neurodegenrative diseases", 
                                                          "Systemic sclerosis", "Inflammatory bowel disease", "Juvenile idiopathic arthritis", "Inflammatory diseases", 
                                                          "Sepsis", "COVID-19", "Infectious diseases"))
plot_data$label <- paste0(plot_data$Intersections, "/", plot_data$Total_genes)

plot_data1 <- rbind(plot_data, disease_category_overlap)
plot_data1[10:13, "Intersections"] <- plot_data1[10:13, "Intersections"] - disease_category_overlap$Intersections

plot_data1 <- left_join(plot_data1, results_df, by = "Disease")
plot_data1$Sig <- ifelse(plot_data1$FDR < 0.0001, "***", 
                         ifelse(plot_data1$FDR < 0.01, "**", 
                                ifelse(plot_data1$FDR < 0.05, "*", "ns")))

plot_df <- plot_data1 %>%
  mutate(Category_base = sub("1$", "", as.character(Category)),
         is_common     = grepl("1$", as.character(Category)), 
         label1        = ifelse(is.na(Sig), label, paste(label, Sig, sep = '\n')))

text_df <- plot_df %>%
  group_by(Disease) %>%
  summarise(y = sum(Intersections, na.rm = TRUE),
            label = first(na.omit(label1)),   # take the non-NA label (e.g., 23/863)
            .groups = "drop")

base_cols <- c(
  "Cardiovascular diseases"   = "#e2b99b",
  "Infectious diseases"       = "#67676c",
  "Inflammatory diseases"     = "#6a9a9d",
  "Neurodegenrative diseases" = "#f57b89"
)

bar_width <- 0.9

disease_genes_plot <- ggplot(plot_df, aes(x = Disease, y = Intersections, fill = Category_base)) +
  geom_col(aes(alpha = is_common), color = "black", width = 0.9, position = "stack") +
  geom_text(data = text_df, aes(x = Disease, y = y+35, label = label), inherit.aes = FALSE, nudge_y = 15, angle = 60, vjust = 1, size = my_axis_text_size/.pt,) +
  scale_fill_manual(values = base_cols, name = "", breaks = c("Cardiovascular diseases", "Neurodegenrative diseases", "Inflammatory diseases", "Infectious diseases")) +
  scale_alpha_manual( values = c(`FALSE` = 0.7, `TRUE` = 1),  # common segment looks darker (less transparent)
                      breaks = "TRUE",
                      labels = "Common genes across different disease categories",
                      name = "") +
  guides( fill  = guide_legend(order = 1, nrow = 2),
          alpha = guide_legend(order = 2, override.aes = list(fill = "grey40", color = "black"))) 

x_labels <- plot_data1 %>%
  distinct(Disease, Group) %>%
  mutate( axis_label = ifelse(Group == "Category",
                              paste0("<b>", Disease, "</b>"),
                              as.character(Disease)))

disease_genes_plot_final <- disease_genes_plot +
  scale_x_discrete(labels = setNames(x_labels$axis_label, x_labels$Disease)) +
  theme_pubr() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.5, "line")) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_blank(),
        panel.grid = element_blank())


## Frequent DEG analysis ------------------------------------------

frequent_degs_table <- read.csv("Frequent_DEGs/pnas.1802973116.sd02.txt", sep = '\t')


## Frequent DEG sets at different DE prior thresholds 

freq_deg_999 <- frequent_degs_table %>%
  filter(DE_Prior_Rank > 0.99) %>%
  pull(Gene_Name) %>% unique()

freq_deg_995 <- frequent_degs_table %>%
  filter(DE_Prior_Rank > 0.95) %>%
  pull(Gene_Name) %>% unique()

freq_deg_90  <- frequent_degs_table %>%
  filter(DE_Prior_Rank > 0.90) %>%
  pull(Gene_Name) %>% unique()

## Outlier genes 

outlier_genes <- icc_results %>%
  filter(Var_diff_outlier == "Outlier genes") %>%
  pull(gene_name) %>% unique()

## Universe of genes (intersection of both datasets) 

universe_genes <- intersect(
  unique(outlier_genes_table$gene_name),
  unique(frequent_degs_table$Gene_Name)
)

length(universe_genes)  # total genes in universe

## Gene classification 

gene_class <- data.frame(
  gene       = universe_genes,
  is_outlier = universe_genes %in% outlier_genes,
  freq_99    = universe_genes %in% freq_deg_999,
  freq_95    = universe_genes %in% freq_deg_995,
  freq_90    = universe_genes %in% freq_deg_90
)




## Helper to compute enrichment stats for one cutoff ----------------------------

compute_enrichment <- function(flag_freq, cutoff_label) {
  tab <- table(
    Outlier  = gene_class$is_outlier,
    Frequent = flag_freq
  )
  ft  <- fisher.test(tab)
  
  # counts
  n_total    <- sum(tab)
  n_outliers <- sum(gene_class$is_outlier)
  n_freq     <- sum(flag_freq)
  n_overlap  <- sum(gene_class$is_outlier & flag_freq)
  
  # proportions
  prop_freq_in_outliers     <- n_overlap / n_outliers
  prop_freq_in_non_outliers <- (n_freq - n_overlap) / (n_total - n_outliers)
  prop_outliers_in_freq     <- n_overlap / n_freq
  prop_outliers_in_nonfreq  <- (n_outliers - n_overlap) / (n_total - n_freq)
  
  data.frame(
    DE_Prior_cutoff              = cutoff_label,
    Total_genes                  = n_total,
    Outlier_genes                = n_outliers,
    Frequent_genes               = n_freq,
    Overlap_genes                = n_overlap,
    Prop_frequent_in_outliers    = prop_freq_in_outliers,
    Prop_frequent_in_nonoutliers = prop_freq_in_non_outliers,
    Prop_outliers_in_frequent    = prop_outliers_in_freq,
    Prop_outliers_in_nonfreq     = prop_outliers_in_nonfreq,
    Odds_ratio                   = as.numeric(ft$estimate),
    CI_lower                     = ft$conf.int[1],
    CI_upper                     = ft$conf.int[2],
    P_value                      = ft$p.value
  )
}

## Summary table for 0.99, 0.95, 0.90 ------------------------------------------

summary_table <- bind_rows(
  compute_enrichment(gene_class$freq_99, "0.99"),
  compute_enrichment(gene_class$freq_95, "0.95"),
  compute_enrichment(gene_class$freq_90, "0.90")
)

summary_table


write.csv(summary_table,
          file = "results/Frequent_DEGs_analysis/frequentDEG_outlier_overlap.csv",
          row.names = FALSE)

summary_table$DE_Prior_cutoff <- factor(
  summary_table$DE_Prior_cutoff,
  levels = c("0.99", "0.95", "0.90")
)

##Enrichment plot as in Fig. 6g -------------------

gg_or <- ggplot(summary_table,
                aes(x = DE_Prior_cutoff,
                    y = Odds_ratio)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "DE prior cutoff",
       y = "Odds ratio (log scale)",
       title = "Enrichment of frequent DEGs among outlier genes") +
  theme_bw(base_size = 12)

##Upset plot as in Fig. Extended Data Fig. 9g ---------------------

# Convert to binary list format for UpSetR
upset_input <- list(
  Outliers      = gene_class$gene[gene_class$is_outlier],
  Frequent_99   = gene_class$gene[gene_class$freq_99],
  Frequent_95   = gene_class$gene[gene_class$freq_95],
  Frequent_90   = gene_class$gene[gene_class$freq_90]
)

# Convert to UpSet matrix
upset_data <- fromList(upset_input)

print(upset(upset_data,
            nsets = 4,
            order.by = "freq",
            mainbar.y.label = "Intersection size",
            sets.x.label = "Set size",
            text.scale = 1.3))

##GO dot plot as in Extended Data Fig. 9d -----------------

go_plot_data <- read.csv("10.input/Outlier_frequent_GO_plot_data.txt", sep = '\t')

go_plot_data$Term <- str_wrap(go_plot_data$Term, width = 45) 
go_plot_data$Group <- factor(go_plot_data$Group, levels = c("Outlier", "Frequent"))

term_first_group <- go_plot_data %>%
  group_by(Term) %>%
  summarise(
    first_group_idx = min(as.integer(Group)[!is.na(value)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(first_group_idx = ifelse(is.infinite(first_group_idx), length(groups) + 1, first_group_idx))

term_rank <- go_plot_data %>%
  group_by(Term) %>%
  summarise(best_p = min(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(best_p = ifelse(is.infinite(best_p), 1, best_p))

term_order <- term_first_group %>%
  left_join(term_rank, by = "Term") %>%
  arrange(first_group_idx, best_p) %>%          
  pull(Term)

go_plot_data$Term2 <- factor(go_plot_data$Term, levels = term_order)

GO_dotplot <- ggplot(go_plot_data, mapping = aes(x = Group, y = Term2, size = n, color = p)) +
  geom_point() +
  scale_colour_gradient(high="#990000", low="#FF9999", name = expression(-~log["10"]~italic(p))) + 
  scale_size(name = "Gene ratio", range = c(1, 3)) +
  guides(color = guide_colourbar(title.position = "top"),
         size = guide_legend(title.position = "top", ncol = 2)) +
  theme_bw() +
  theme(axis.title.x = element_text(hjust = 0.5),
        strip.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.direction = "horizontal")


## eQTL discovery rate and heritability of outlier genes ------------------

eqtl_table <- read.csv("results/eQTL_based_analysis/Whole_dataset/GTEx_v10/ICC_cis_eGene_data.txt", sep = '\t')
heritability_table <- read.csv("results/Heritability_analysis/Whole_dataset/ICC_Heritability_table.txt", sep = '\t')

gene_metrics <- gene_class %>%
  left_join(., eqtl_table[, c("gene_name", "ICC", "eGene", "eGene_binary", "slope", "n_eQTL")], by = c("gene"="gene_name")) %>%
  left_join(., heritability_table[, c("gene_name", "h2", "e2", "Heritable")], by = c("gene"="gene_name"))

gene_metrics <- gene_metrics %>%
  mutate(
    Category = case_when(
      is_outlier ~ "Outlier",
      TRUE ~ "Other"
    ),
    Category = factor(Category,
                      levels = c("Outlier",
                                 "Other"))
  )

summary_overlap <- gene_metrics %>%
  group_by(Category) %>%
  summarise(
    n_genes        = n(),
    ICC_median     = median(ICC, na.rm = TRUE),
    ICC_IQR        = IQR(ICC, na.rm = TRUE),
    
    eQTL_rate      = mean(eGene_binary, na.rm = TRUE),    # ← correct
    eQTL_rate_CI_l = binom.test(sum(eGene_binary, na.rm=TRUE),
                                n())$conf.int[1],
    eQTL_rate_CI_u = binom.test(sum(eGene_binary, na.rm=TRUE),
                                n())$conf.int[2],
    
    h2_median      = median(h2, na.rm = TRUE),
    h2_IQR         = IQR(h2, na.rm = TRUE),
    Heritable_rate = mean(Heritable == "YES", na.rm = TRUE)
  )

tab <- gene_metrics %>%
  count(Category, eGene_binary) %>%
  tidyr::pivot_wider(names_from = eGene_binary,
                     values_from = n,
                     values_fill = 0)

mat <- as.matrix(tab[, c("0", "1")])
rownames(mat) <- tab$Category

chisq_res <- chisq.test(mat, correct = FALSE)

pval <- chisq_res$p.value

# convert p-value to stars
pval_df <- data.frame(
  group1 = "Outlier",
  group2 = "Other",
  p = pval,
  p.signif = cut(
    pval,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  ),
  y.position = 0.95   # place stars above bars
)

my_colors <- c(
  "Other"              = "darkgray",
  "Outlier"       = "skyblue2"
)

## Bar plot for eQTL discovery rates as in Extended Data Fig. 9e -------------

gg_eqtl_bar <- ggplot(summary_overlap,
                      aes(x = Category, y = eQTL_rate, fill = Category, col = Category)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = eQTL_rate_CI_l, ymax = eQTL_rate_CI_u),
                width = 0.2, linewidth = 0.5) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  stat_pvalue_manual(
    pval_df,
    label = "p.signif",
    tip.length = 0.01,
    size = 5, 
    inherit.aes = FALSE
  ) +
  labs(
    x = "",
    y = "cis-eQTL discovery rate (% eGenes)",
    title = "eQTL discovery across gene categories"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

## Violin plot plot for heritability as in Extended Data Fig. 9f -------------

gg_heritability <- ggplot(gene_metrics,
                          aes(x = Category, y = h2, fill = Category)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", 
    hide.ns = TRUE
  ) + 
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  labs(title = "Heritability across gene categories", y = "Heritability", x = "")


