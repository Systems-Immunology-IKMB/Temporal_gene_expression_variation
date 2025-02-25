## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(DESeq2)
library(reshape2)
library(dplyr)
library(ICC)
library(rstatix)
library(ggpubr)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


## Loading data ------------------------------------------------------------

gene_attributes <- read.csv("info_files/filtered_genes_samples25_fpm5.txt", sep = '\t')
annotation <- gene_attributes[, c("Gene_ID", "gene_name")] %>% unique()
sex_chromosome_genes <- subset(gene_attributes, gene_attributes$Chr %in% c("chrX", "chrY"))$Gene_ID
gene_attributes <- gene_attributes %>%
  subset(!.$Gene_ID %in% sex_chromosome_genes)

col_clinical_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_clinical_data$Individual_ID <- factor(col_clinical_data$Individual_ID)
rownames(col_clinical_data) <- col_clinical_data$Sample_ID
col_clinical_data$Run <- factor(col_clinical_data$Run)

#import count data
count_data <- read.csv("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)
count_data <- count_data[as.character(gene_attributes$Gene_ID), as.character(col_clinical_data$Sample_ID)]





## Compute intraclass correlation coefficient (ICC) for each gene ------------------------

#normalize counts
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ 1)
dds_counts <- estimateSizeFactors(dds_counts)
quantvst <- vst(dds_counts)
quantvst <- assay(quantvst)

#compute icc
gene_icc_data <- data.frame()
for (i in 1:nrow(quantvst)) {
  test_data <- data.frame(counts = quantvst[i,], 
                          individual = col_clinical_data$Individual_ID, 
                          timepoint = col_clinical_data$Timepoint)
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

gene_icc_data <- left_join(gene_icc_data, gene_attributes, by = "Gene_ID")

write.table(gene_icc_data, "results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t', row.names = FALSE, quote = FALSE)





## Boxplots for inter- and intra-individual variance in gene expression as in Fig. 2a ---------

gene_icc_data <- read.csv("results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t')
gene_icc_data$Category <- ifelse(gene_icc_data$gene_type == "protein_coding", "Protein coding", "Non protein coding")
gene_icc_data$Group <- ifelse(gene_icc_data$ICC < 0.5, "Within", "Between")

plot_data <- melt(gene_icc_data[, c("gene_name", "Var_within", "Var_between", "Category")], id.vars = c("gene_name", "Category"))

wilcox_test <- plot_data %>%
  group_by(Category) %>%
  wilcox_test(value ~ variable, data = . , paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = 0.19)

boxplots <- ggplot(plot_data) + 
  geom_boxplot(mapping = aes(x=variable, y=value, fill = variable)) +
  facet_wrap(~ Category, scales = "free") + 
  scale_y_continuous(limits = quantile(plot_data$value, c(0.1, 0.9))) + 
  scale_x_discrete(labels = c("Intra-individual", "Inter-individual")) + 
  xlab("") + 
  ylab("Variance") + 
  scale_fill_manual(values = c('#CCCC00', '#009900'), name = "", labels = c("Var_within"="Intra-individual", "Var_between"="Inter-individual")) + 
  stat_pvalue_manual(data = wilcox_test, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) + 
  theme_bw() +
  theme(strip.background = element_blank(), 
        panel.grid = element_blank(),  
        legend.position = "top",
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())



## Dotplot for comparison of inter- and intra-individual variance as in Fig. 2b ---------

gene_icc_data_protein_coding <- subset(gene_icc_data, gene_icc_data$Category == "Protein coding")

scatter_plot <- ggplot(gene_icc_data_protein_coding, aes(x=Var_within, y=Var_between, col = Group)) + 
  geom_point(size = my_point_size/2) +
  geom_text_repel(data = gene_icc_data_protein_coding[gene_icc_data_protein_coding$Var_between > 2 | gene_icc_data_protein_coding$Var_within > 1, ],
                  mapping = aes(x=Var_within, y=Var_between, label=gene_name), 
                  col="black",
                  size = 1.7, 
                  min.segment.length = 0) +
  scale_color_manual(values = c('#009900', '#CCCC00'), name = "", labels = c("Between individuals", "Within individuals")) +
  ylab("Inter-individual variance") + 
  xlab("Intra-individual variance") +
  xlim(c(0,4)) + 
  ylim(c(0,4)) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        legend.title = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.7, 0.8))

dens1 <- ggplot(gene_icc_data_protein_coding, aes(x = Var_within, fill = Group)) + 
  geom_density(alpha = 0.4, lwd = my_line_width) + 
  scale_fill_manual(values = c('#009900', '#CCCC00')) +
  theme_void() + 
  theme(legend.position = "none")

dens2 <- ggplot(gene_icc_data_protein_coding, aes(x = Var_between, fill = Group)) + 
  geom_density(alpha = 0.4, lwd = my_line_width) + 
  scale_fill_manual(values = c('#009900', '#CCCC00')) + 
  coord_flip() +
  theme_void() + 
  theme(legend.position = "none")





## Density ICC plot as in Fig. 2c -----------------------

non_protein_coding_below_0.5 <- sum(gene_icc_data$Category == "Non protein coding" & gene_icc_data$ICC < 0.5) / sum(gene_icc_data$Category == "Non protein coding")
non_protein_coding_above_0.5 <- sum(gene_icc_data$Category == "Non protein coding" & gene_icc_data$ICC > 0.5) / sum(gene_icc_data$Category == "Non protein coding")
protein_coding_below_0.5 <- sum(gene_icc_data$Category == "Protein coding" & gene_icc_data$ICC < 0.5) / sum(gene_icc_data$Category == "Protein coding")
protein_coding_above_0.5 <- sum(gene_icc_data$Category == "Protein coding" & gene_icc_data$ICC > 0.5) / sum(gene_icc_data$Category == "Protein coding")


labels <- data.frame(
  Category = c("Non protein coding", "Non protein coding", "Protein coding", "Protein coding"),
  ICC_threshold = c("ICC < 0.5", "ICC > 0.5", "ICC < 0.5", "ICC > 0.5"),
  Percentage = c(non_protein_coding_below_0.5, non_protein_coding_above_0.5, protein_coding_below_0.5, protein_coding_above_0.5),
  x = c(0.05, 0.75, 0.05, 0.75),
  y = c(2, 2, 3, 3)
)

icc_density_plot <- ggplot(gene_icc_data, aes(x=ICC, col = Category)) + 
  geom_density(adjust = 1.2) + 
  geom_vline(xintercept = 0.5) +
  ylab("Density (genes)") + 
  scale_color_manual(values = c("#009999", "#FF6666"), name="") + 
  geom_text(data = labels, aes(x = x, y = y, label = scales::percent(Percentage), col = Category), 
            size = 1.7, hjust = 0, show.legend = FALSE) + 
  xlim(c(0, 1)) +
  theme_bw() + 
  theme(strip.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "top",
        panel.border = element_blank(),
        legend.title = element_blank())




## Ridgeplot for ORA results as in Fig. 2d ---------

ridgeplot_data <- read.csv("results/ICC_analysis/output/High_low_ICC_protein_coding_GO_data.txt", sep = '\t') %>%
  add_column("p" = -log10(.$Fisher.elim))

ridgeplot_data$Term <- str_wrap(ridgeplot_data$Term, width = 30)
ridgeplot_data$Term <- factor(ridgeplot_data$Term, levels = ridgeplot_data$Term %>% unique() %>% .[c(1, 6, 5, 2, 3, 4, 10, 7, 9, 8)] %>% rev())

ridgeplot <- ggplot(ridgeplot_data, aes(x = ICC, y = Term, fill = p)) +
  geom_density_ridges() +
  theme_ridges() + 
  scale_fill_gradient(high="#990000", low="#FF9999", name = expression(-~log["10"]~italic(p))) + 
  ylab("") + 
  guides(fill = guide_colourbar(title.position = "top")) +
  theme_bw() +
  theme(axis.title.x = element_text(hjust = 0.5),
        strip.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.direction = "horizontal")




## Boxplots for inter- and intra-individual variance comparison by expression level as in Extended Data Fig. 3a ------------------

gene_icc_data <- read.csv("results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t')
gene_icc_data$Category <- ifelse(gene_icc_data$gene_type == "protein_coding", "Protein coding", "Non protein coding")
gene_icc_data$Group <- ifelse(gene_icc_data$ICC < 0.5, "Within", "Between")
gene_icc_data_protein_coding <- subset(gene_icc_data, gene_icc_data$Category == "Protein coding")

gene_icc_data_protein_coding$Expression_group <- cut_number(gene_icc_data_protein_coding$Mean_expression, 3)

plot_data <- melt(gene_icc_data_protein_coding[, c("gene_name", "Expression_group","Var_within", "Var_between")], id.vars = c("gene_name", "Expression_group"))
expression_group_labels <- c("Low", "Mid", "High")
names(expression_group_labels) <- levels(plot_data$Expression_group)

wilcox_paired_fn <- function(df){
  wt <- wilcox.test(df$Var_within, df$Var_between, paired = TRUE)
  return(data.frame(group1 = "Var_within", group2 = "Var_between", statistic = wt$statistic, p = wt$p.value))
}

wilcox_test <- gene_icc_data_protein_coding %>%
  group_by(Expression_group) %>%
  do(wilcox_paired_fn(.)) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = 0.19)

var_by_expr_plot <- ggplot(plot_data) + 
  geom_boxplot(mapping = aes(x=variable, y=value, fill = variable)) +
  facet_wrap(~Expression_group, labeller = labeller(Expression_group = expression_group_labels)) +
  scale_y_continuous(limits = quantile(plot_data$value, c(0.1, 0.9))) +
  xlab("") +
  ylab("Variance") +
  scale_fill_manual(values = c('#CCCC00', '#009900'), name = "", labels = c("Var_within"="Intra-individual", "Var_between"="Inter-individual")) +
  stat_pvalue_manual(data = wilcox_test, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5, size = my_signif_size) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())







## Plot ICC for cell type proportions as in Extended Data Fig. 3b ----------------------

#compute ICC for cell type proportions
df_icc_cell_types <- col_clinical_data[, c("Individual_ID", "Sample_ID", "Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL.")] %>%
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

#prepare table for plot
df_plot <- df_icc %>%
  .[order(.$ICC, decreasing = FALSE), ] %>%
  add_column("within" = 1 - .$ICC, "between" = .$ICC) %>%
  melt(id.vars = c("Name"), measure.vars = c("within", "between"), value.name = "ICC") %>%
  add_column("Group" = paste0(.$Name))

df_plot$variable <- factor(df_plot$variable, levels = c("between", "within"))
df_plot$Name <- factor(df_plot$Name, levels = df_plot$Name %>% unique())
my_levels <- df_plot$Name %>% unique()

df_errorbars <- df_icc %>%
  add_column("CI_lower" = 1 - .$lower_ICC,
             "CI_upper" = 1 - .$upper_ICC)
df_errorbars$Name <- factor(df_errorbars$Name, levels = df_plot$Name %>% unique() %>% rev())

df_text <- df_plot %>%
  add_column("x" = .$ICC/2,
             "rounded_ICC" = .$ICC %>% round(2) %>% `*` (100))
df_text$x <- 0.05
df_text[df_text$variable == "between", ]$x <- 0.95

my_alpha <- 0.7

ICC_cell_types_plot <- ggplot(df_plot, aes(x = Name, y = ICC)) +
  geom_bar(position = "fill", stat = "identity", aes(fill = variable)) +
  geom_hline(yintercept = 0.5, color = "white") +
  geom_errorbar(data = df_errorbars, aes(ymin = CI_upper, ymax = CI_lower), width = 0.3, color = "black") +
  geom_text(data = df_text, aes(x = Name, y = x, label = rounded_ICC), color = "black") +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  labs(y = "Part of total variance (%)") +
  scale_fill_manual(breaks = c("within", "between"), values = c('#CCCC00', '#009900'), labels = c("Within individual", "Between individuals") ) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  coord_flip()







## Compute ICC for sub-cohorts (only healthy participants; only samples collected in the morning) -----------------------

#Compute ICC on samples collected in the morning
#only participants with at least two samples collected between 7 and 10am are included in this analysis
col_data_morning <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_data_morning <- subset(col_data_morning, col_data_morning$morning_sample == "YES")
n_morning_samples <- col_data_morning %>% group_by(Individual_ID) %>% dplyr::count() %>% data.frame()
col_data_morning <- subset(col_data_morning, col_data_morning$Individual_ID %in% n_morning_samples[n_morning_samples$n > 1, "Individual_ID"])

quantvst_morning <- quantvst[, as.character(col_data_morning$Sample_ID)]

gene_icc_data_morning <- data.frame()
for (i in 1:nrow(quantvst_morning)) {
  test_data <- data.frame(counts = quantvst_morning[i,], 
                          individual = col_data_morning$Individual_ID, 
                          timepoint = col_data_morning$Timepoint)
  total_var <- var(test_data$counts)
  res <- ICCest(test_data$individual, test_data$counts)
  gene_icc_data_morning <- rbind(gene_icc_data_morning, data.frame(Gene_ID = rownames(quantvst_morning)[i], 
                                                   ICC = res$ICC, 
                                                   CI_lower = res$LowerCI, 
                                                   CI_upper = res$UpperCI, 
                                                   Var_within = res$varw, 
                                                   Var_between = res$vara, 
                                                   Total_var = total_var, 
                                                   Mean_expression = mean(test_data$counts)))
}

gene_icc_data_morning <- left_join(gene_icc_data_morning, gene_attributes, by = "Gene_ID")
write.table(gene_icc_data_morning, "results/ICC_analysis/output/gene_icc_data_morning_samples.txt", sep = '\t', row.names = FALSE, quote = FALSE)



#Compute ICC on samples from healthy participants
#only participants without acute or ongoing healthy conditions are included in this analysis
col_data_healthy <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_data_healthy <- subset(col_data_healthy, col_data_healthy$healthy == "YES")

quantvst_healthy <- quantvst[, as.character(col_data_healthy$Sample_ID)]

gene_icc_data_healthy <- data.frame()
for (i in 1:nrow(quantvst_healthy)) {
  test_data <- data.frame(counts = quantvst_healthy[i,], 
                          individual = col_data_healthy$Individual_ID, 
                          timepoint = col_data_healthy$Timepoint)
  total_var <- var(test_data$counts)
  res <- ICCest(test_data$individual, test_data$counts)
  gene_icc_data_healthy <- rbind(gene_icc_data_healthy, data.frame(Gene_ID = rownames(quantvst_healthy)[i], 
                                                                   ICC = res$ICC, 
                                                                   CI_lower = res$LowerCI, 
                                                                   CI_upper = res$UpperCI, 
                                                                   Var_within = res$varw, 
                                                                   Var_between = res$vara, 
                                                                   Total_var = total_var, 
                                                                   Mean_expression = mean(test_data$counts)))
}

gene_icc_data_healthy <- left_join(gene_icc_data_healthy, gene_attributes, by = "Gene_ID")
write.table(gene_icc_data_healthy, "results/ICC_analysis/output/gene_icc_data_healthy.txt", sep = '\t', row.names = FALSE, quote = FALSE)






## Scatterplot for ICC comparisons between sub-cohorts as in Extended Data Fig. 3c -----------------------

#load data
gene_icc_data2 <- read.csv("results/ICC_analysis/output/gene_icc_data_healthy.txt", sep = '\t')
gene_icc_data2$Category <- ifelse(gene_icc_data2$gene_type == "protein_coding", "Protein coding", "Non protein coding")
gene_icc_data2$Group <- ifelse(gene_icc_data2$ICC < 0.5, "Within", "Between")

gene_icc_data3 <- read.csv("results/ICC_analysis/output/gene_icc_data_morning_samples.txt", sep = '\t')
gene_icc_data3$Category <- ifelse(gene_icc_data3$gene_type == "protein_coding", "Protein coding", "Non protein coding")
gene_icc_data3$Group <- ifelse(gene_icc_data3$ICC < 0.5, "Within", "Between")

#whole dataset vs healthy
gene_icc_data_combined12 <- full_join(gene_icc_data, gene_icc_data2, by = "Gene_ID", suffix = c(".whole", ".healthy"))

label_data1 <- gene_icc_data_combined12[abs(gene_icc_data_combined12$ICC.whole - gene_icc_data_combined12$ICC.healthy) > 0.15 &
                                          ((gene_icc_data_combined12$ICC.whole > 0.5 & gene_icc_data_combined12$ICC.healthy < 0.5) | 
                                             (gene_icc_data_combined12$ICC.whole < 0.5 & gene_icc_data_combined12$ICC.healthy > 0.5)), ]

p1 <- ggplot(gene_icc_data_combined12, aes(x=ICC.whole, y=ICC.healthy)) + 
  geom_point(col="skyblue2", alpha = 0.5) +
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) +
  geom_text_repel(data = label_data1, 
                  aes(x = ICC.whole, y = ICC.healthy, label = gene_name.whole), 
                  color = "black", max.overlaps = 20, size = 1.7) +
  stat_cor(method = "spearman", size = 1.7, cor.coef.name = "rho") +
  xlab("ICC (All samples)") + 
  ylab("ICC (Without medical conditions)") +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank())


#whole dataset vs morning samples
gene_icc_data_combined13 <- full_join(gene_icc_data, gene_icc_data3, by = "Gene_ID", suffix = c(".whole", ".morning"))

label_data2 <- gene_icc_data_combined13[abs(gene_icc_data_combined13$ICC.whole - gene_icc_data_combined13$ICC.morning) > 0.15 &
                                          ((gene_icc_data_combined13$ICC.whole > 0.5 & gene_icc_data_combined13$ICC.morning < 0.5) | 
                                             (gene_icc_data_combined13$ICC.whole < 0.5 & gene_icc_data_combined13$ICC.morning > 0.5)), ]

p2 <- ggplot(gene_icc_data_combined13, aes(x=ICC.whole, y=ICC.morning)) + 
  geom_point(col="skyblue2", alpha = 0.5) + 
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) + 
  geom_text_repel(data = label_data2, 
                  aes(x = ICC.whole, y = ICC.morning, label = gene_name.whole), 
                  color = "black", max.overlaps = 20, size = 1.7) + 
  stat_cor(method = "spearman", size = 1.7, cor.coef.name = "rho") + 
  xlab("ICC (All samples)") + 
  ylab("ICC (Morning samples)") + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank())






