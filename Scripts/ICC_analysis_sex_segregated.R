## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(DESeq2)
library(reshape2)
library(dplyr)
library(ICC)
library(limma)
library(ggpubr)
library(rstatix)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")



## Loading data ------------------------------------------------------------

gene_list <- read.csv("info_files/results_Neha/filtered_genes_samples25_fpm5.txt", sep = '\t')
gene_list <- subset(gene_list, !(gene_list$Chr %in% c("chrY")))
rownames(gene_list) <- gene_list$Gene_ID




## Sex segregated ICC analyses for men and women on age adjusted counts -------------------------

col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')

#load and adjust counts for age 
count_data <- read.csv("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)

count_data <- count_data[as.character(gene_list$Gene_ID), as.character(col_data$Sample_ID)]
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ 1)
dds_counts <- estimateSizeFactors(dds_counts)
quantvst <- vst(dds_counts)
quantvst <- assay(quantvst)

fit <- lmFit(quantvst, model.matrix(~ Age..years., col_data))
quantvst_residuals <- residuals(fit, quantvst)



### ICC analysis on men, morning samples -------------------------

#subset metadata
col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_data <- subset(col_data, col_data$Gender == "Man")
col_data <- subset(col_data, col_data$morning_sample == "YES")
n_morning_samples <- col_data %>% group_by(Individual_ID) %>% dplyr::count() %>% data.frame()
col_data <- subset(col_data, col_data$Individual_ID %in% n_morning_samples[n_morning_samples$n > 1, "Individual_ID"])

#run analysis
quantvst_residuals_temp <- quantvst_residuals[, col_data$Sample_ID]

gene_icc_data <- data.frame()
for (i in 1:nrow(quantvst_residuals_temp)) {
  
  test_data <- data.frame(counts = quantvst_residuals_temp[i,], 
                          individual = col_data$Individual_ID, 
                          timepoint = col_data$Timepoint)
  total_var <- var(test_data$counts)
  res <- ICCest(test_data$individual, test_data$counts)
  gene_icc_data <- rbind(gene_icc_data, data.frame(Gene_ID = rownames(quantvst_residuals_temp)[i], 
                                                   ICC = res$ICC, 
                                                   CI_lower = res$LowerCI, 
                                                   CI_upper = res$UpperCI, 
                                                   Var_within = res$varw, 
                                                   Var_between = res$vara, 
                                                   Total_var = total_var, 
                                                   Mean_expression = mean(test_data$counts)))
}

gene_icc_data <- left_join(gene_icc_data, gene_list, by = "Gene_ID")
write.table(gene_icc_data, "results/ICC_sex_differences/output/ICC_results_males_ageAdj.txt", sep = '\t', row.names = FALSE, quote = FALSE)





### ICC analysis on women, morning samples -------------------------

#subset metadata
col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_data <- subset(col_data, col_data$Gender == "Vrouw")
col_data <- subset(col_data, col_data$morning_sample == "YES")
n_morning_samples <- col_data %>% group_by(Individual_ID) %>% dplyr::count() %>% data.frame()
col_data <- subset(col_data, col_data$Individual_ID %in% n_morning_samples[n_morning_samples$n > 1, "Individual_ID"])

#run analysis
quantvst_residuals_temp <- quantvst_residuals[, col_data$Sample_ID]

gene_icc_data <- data.frame()
for (i in 1:nrow(quantvst_residuals_temp)) {
  
  test_data <- data.frame(counts = quantvst_residuals_temp[i,], 
                          individual = col_data$Individual_ID, 
                          timepoint = col_data$Timepoint)
  total_var <- var(test_data$counts)
  res <- ICCest(test_data$individual, test_data$counts)
  gene_icc_data <- rbind(gene_icc_data, data.frame(Gene_ID = rownames(quantvst_residuals_temp)[i], 
                                                   ICC = res$ICC, 
                                                   CI_lower = res$LowerCI, 
                                                   CI_upper = res$UpperCI, 
                                                   Var_within = res$varw, 
                                                   Var_between = res$vara, 
                                                   Total_var = total_var, 
                                                   Mean_expression = mean(test_data$counts)))
}

gene_icc_data <- left_join(gene_icc_data, gene_list, by = "Gene_ID")
write.table(gene_icc_data, "results/ICC_sex_differences/output/ICC_results_females_ageAdj.txt", sep = '\t', row.names = FALSE, quote = FALSE)





### ICC density plot as in Fig. 2e -----------------

gene_icc_data_males <- read.table("results/ICC_sex_differences/output/ICC_results_males_ageAdj.txt", sep = '\t', header = TRUE)
gene_icc_data_females <- read.table("results/ICC_sex_differences/output/ICC_results_females_ageAdj.txt", sep = '\t', header = TRUE)

merged_results_ageAdj <- gene_icc_data_males %>%
  rbind(gene_icc_data_females) %>%
  add_column("Category" = ifelse(.$gene_type == "protein_coding", "Protein coding", "Non protein coding"))

male_non_protein_coding_below_0.5 <- sum(merged_results_ageAdj$Category == "Non protein coding" & merged_results_ageAdj$ICC.males < 0.5) / sum(merged_results_ageAdj$Category == "Non protein coding")
male_non_protein_coding_above_0.5 <- sum(merged_results_ageAdj$Category == "Non protein coding" & merged_results_ageAdj$ICC.males > 0.5) / sum(merged_results_ageAdj$Category == "Non protein merged_results_ageAdj")
male_protein_coding_below_0.5 <- sum(merged_results_ageAdj$Category == "Protein coding" & merged_results_ageAdj$ICC.males < 0.5) / sum(merged_results_ageAdj$Category == "Protein coding")
male_protein_coding_above_0.5 <- sum(merged_results_ageAdj$Category == "Protein coding" & merged_results_ageAdj$ICC.males > 0.5) / sum(merged_results_ageAdj$Category == "Protein coding")

female_non_protein_coding_below_0.5 <- sum(merged_results_ageAdj$Category == "Non protein coding" & merged_results_ageAdj$ICC.females < 0.5) / sum(merged_results_ageAdj$Category == "Non protein coding")
female_non_protein_coding_above_0.5 <- sum(merged_results_ageAdj$Category == "Non protein coding" & merged_results_ageAdj$ICC.females > 0.5) / sum(merged_results_ageAdj$Category == "Non protein coding")
female_protein_coding_below_0.5 <- sum(merged_results_ageAdj$Category == "Protein coding" & merged_results_ageAdj$ICC.females < 0.5) / sum(merged_results_ageAdj$Category == "Protein coding")
female_protein_coding_above_0.5 <- sum(merged_results_ageAdj$Category == "Protein coding" & merged_results_ageAdj$ICC.females > 0.5) / sum(merged_results_ageAdj$Category == "Protein coding")

labels <- data.frame(
  Sex = c(rep("Male", 4), rep("Female", 4)), 
  Category = rep(c("Non protein coding", "Non protein coding", "Protein coding", "Protein coding"), 2),
  ICC_threshold = rep(c("ICC < 0.5", "ICC > 0.5", "ICC < 0.5", "ICC > 0.5"), 2),
  Percentage = c(male_non_protein_coding_below_0.5, male_non_protein_coding_above_0.5, male_protein_coding_below_0.5, male_protein_coding_above_0.5, 
                 female_non_protein_coding_below_0.5, female_non_protein_coding_above_0.5, female_protein_coding_below_0.5, female_protein_coding_above_0.5),
  x = rep(c(0, 0.75, 0, 0.75), 2),
  y = c(2, 2, 2.7, 2.7, 2.2, 2.2, 3, 3)
)

gene_icc_data_sex <- merged_results_ageAdj[, c(1, 2, 17, 32)] %>%
  melt(id.vars = c("Gene_ID", "Category"), value.name = "ICC") %>%
  add_column("Sex" = ifelse(.$variable == "ICC.males", "Male", "Female"))

density_plot_sex <- ggplot(gene_icc_data_sex, aes(x = ICC, col = Sex)) + 
  geom_density() + 
  ylab("Density (genes)") +
  facet_wrap(~ Category, scales = "free", ncol = 1, strip.position = "right") +
  geom_vline(xintercept = 0.5) +
  geom_text(data = labels, aes(x = x, y = y, label = scales::percent(Percentage)), 
            hjust = 0, size = 1.7, show.legend = FALSE) +
  scale_color_manual(values = c("Male"="#9933FF", "Female"="#FF9933"), name="") +
  xlim(c(0, 1)) +
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(), 
        panel.border = element_blank(),
        panel.grid = element_blank())




### ICC boxplots as in Fig. 2f -----------------

plot_data <- merged_results_ageAdj[, c("gene_name.males", "Category", "Var_within.males", "Var_between.males", "Var_within.females", "Var_between.females")] %>%
  melt(id.vars = c("gene_name.males", "Category")) %>%
  add_column("Sex" = ifelse(.$variable %in% c("Var_within.males", "Var_between.males"), "Male", "Female")) %>%
  mutate("variable" = ifelse(.$variable %in% c("Var_within.males", "Var_within.females"), "Var_within", "Var_between") %>%
           factor(levels = c("Var_within", "Var_between"))) 

wilcox_test <- plot_data %>%
  group_by(Category, variable) %>%
  wilcox_test(value ~ Sex, data = . , paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = 0.19)

boxplots_sex <- ggplot(plot_data) + 
  geom_boxplot(mapping = aes(x=Sex, y=value, fill = Sex)) +
  facet_grid(Category ~ variable, labeller = labeller(variable = c("Var_within"="Within\nindividuals", "Var_between"="Between\nindividuals"))) +
  xlab("") + 
  ylab("Variance") +
  scale_fill_manual(values = c("Male"="#9933FF", "Female"="#FF9933"), name="") +
  stat_pvalue_manual(data = wilcox_test, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
  scale_y_continuous(limits = quantile(plot_data$value, c(0.1, 0.9))) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(), 
        legend.key.size = my_legend_size,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")



### Scatter plot for sex-specific ICC comparison as in Extended Data Fig. 5a -------------------

label_data <- merged_results_ageAdj %>%
  .[abs(.$ICC.males - .$ICC.females) > 0.2 &
      ((.$ICC.males > 0.5 & .$ICC.females < 0.5) | 
         (.$ICC.males < 0.5 & .$ICC.females > 0.5)), ]
label_data1 <- rbind(label_data[order(label_data$ICC.males, decreasing = TRUE), ][1:15,], label_data[order(label_data$ICC.females, decreasing = TRUE), ][1:10,])

scatter_plot_sex <- ggplot(merged_results_ageAdj, aes(x=ICC.males, y=ICC.females)) + 
  geom_point(alpha = 0.5, col="skyblue2") + 
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) + 
  xlab("ICC (males)") + 
  ylab("ICC (females)") + 
  geom_text_repel(data = label_data1, 
                         aes(x = ICC.males, y = ICC.females, label = gene_name.males), 
                         color = "black", max.overlaps = 100, size = 1.7) + 
  stat_cor(method = "spearman", size = 1.7, cor.coef.name = "rho") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line())



### Functional enrichment dotplot as in Extended Data Fig. 5e ---------------------------

direction <- c("conserved_in_males", "conserved_in_females")
d <- 1
for(d in 1:length(direction)){
  
  if(direction[d] == "conserved_in_males") {
    degs <- merged_results_ageAdj %>%
      subset(.$ICC_diff > 0.2) %>%
      .[order(.$ICC_diff, decreasing = TRUE), ] %>%
      .$Gene_ID
  }else if(direction[d] == "conserved_in_females") {
    degs <- merged_results_ageAdj %>%
      subset(.$ICC_diff < -0.2) %>%
      .[order(.$ICC_diff, decreasing = FALSE), ] %>%
      .$Gene_ID
  }
  
  overallBaseMean <- as.matrix(base_mean[, "mean_expression", drop = F]) 
  topGOResults <- run_GO_ORA(degs, 
                             overallBaseMean, 
                             selected_database = "org.Hs.eg.db", 
                             back_num = 10,
                             annotation = annotation)
  
  go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
  
  file_filtered <- paste0("results/ICC_sex_differences/output/ICC_ORA_GO_", direction[d], ".csv")
  write.csv(go_results_filtered, file_filtered, row.names = FALSE)
}

#load data
direction <- direction <- c("conserved_in_males", "conserved_in_females")
GO_tables <- list()
d <- 1
for(d in 1:length(direction)){
  
  GO_tables <- append(GO_tables,
                      list(read_csv(paste0("results/ICC_sex_differences/output/ICC_ORA_GO_", direction[d], "ICC_diff_0.2.csv"))))
  
}

my_terms_1 <- GO_tables[[1]][c(4, 5, 6, 7, 8, 12, 14, 16, 20), ]$GO.ID
my_terms_2 <- GO_tables[[2]][c(3, 7, 10, 14, 15, 20, 21), ]$GO.ID
my_terms <- c(my_terms_1, my_terms_2)

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = c("conserved_in_males", "conserved_in_females"),
                     score_name = score_name,
                     reverse_bool = TRUE,
                     showTerms = my_terms,
                     multLines = TRUE,
                     numChar = 80)

GO_dotplot_sex_diff <- ggplot(df_plot, mapping = aes(x = Direction, y = Term3, size = Ratio, color = changed_scores)) +
  geom_point() +
  xlab("") + ylab("") +
  scale_colour_gradient(high="#990000", low="#FF9999", name = expression(-~log["10"]~italic(p))) + 
  scale_size(name = "Gene ratio", range = c(0, 3)) +
  guides(color = guide_colourbar(title.position = "top"),
         size = guide_legend(title.position = "top", ncol = 2)) +
  theme_bw() + 
  theme(axis.text.y = element_text(color = "black"), 
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1), 
        legend.text = element_text(),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        panel.grid = element_line(),
        plot.margin = unit(c(0, 0.3, 0, 0.3), "cm"))



## Sex segregated ICC analyses for women, pre- and post-menopausal women (not age adjusted) -------------------------

### ICC analysis on men, morning samples -------------------------

#subset metadata
col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_data <- subset(col_data, col_data$Gender == "Man")
col_data <- subset(col_data, col_data$morning_sample == "YES")
n_morning_samples <- col_data %>% group_by(Individual_ID) %>% dplyr::count() %>% data.frame()
col_data <- subset(col_data, col_data$Individual_ID %in% n_morning_samples[n_morning_samples$n > 1, "Individual_ID"])

#run analysis
count_data_temp <- count_data[as.character(gene_list$Gene_ID), as.character(col_data$Sample_ID)]
dds_counts <- DESeqDataSetFromMatrix(countData = count_data_temp, colData = col_data, design = ~ 1)
dds_counts <- estimateSizeFactors(dds_counts)
quantvst <- vst(dds_counts)
quantvst <- assay(quantvst)

gene_icc_data <- data.frame()
for (i in 1:nrow(quantvst)) {
  
  test_data <- data.frame(counts = quantvst[i,], 
                          individual = col_data$Individual_ID, 
                          timepoint = col_data$Timepoint)
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
write.table(gene_icc_data, "results/ICC_sex_differences/output/ICC_results_males.txt", sep = '\t', row.names = FALSE, quote = FALSE)




### ICC analysis on postmenopausal women, morning samples -------------------------

#subset metadata
col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_data <- subset(col_data, col_data$Gender == "Vrouw")
col_data <- subset(col_data, col_data$Individual_ID %in% sample_female_menopause)
col_data <- subset(col_data, col_data$morning_sample == "YES")
n_morning_samples <- col_data %>% group_by(Individual_ID) %>% dplyr::count() %>% data.frame()
col_data <- subset(col_data, col_data$Individual_ID %in% n_morning_samples[n_morning_samples$n > 1, "Individual_ID"])

#run analysis
count_data_temp <- count_data[as.character(gene_list$Gene_ID), as.character(col_data$Sample_ID)]
dds_counts <- DESeqDataSetFromMatrix(countData = count_data_temp, colData = col_data, design = ~ 1)
dds_counts <- estimateSizeFactors(dds_counts)
quantvst <- vst(dds_counts)
quantvst <- assay(quantvst)

gene_icc_data <- data.frame()
for (i in 1:nrow(quantvst)) {
  
  test_data <- data.frame(counts = quantvst[i,], 
                          individual = col_data$Individual_ID, 
                          timepoint = col_data$Timepoint)
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
write.table(gene_icc_data, "results/ICC_sex_differences/output/ICC_results_postmenopausal_females.txt", sep = '\t', row.names = FALSE, quote = FALSE)





### ICC analysis on premenopausal women, morning samples -------------------------

#subset metadata
col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
col_data <- subset(col_data, col_data$Gender == "Vrouw")
col_data <- subset(col_data, col_data$Individual_ID %in% sample_female_premenopause)
col_data <- subset(col_data, col_data$morning_sample == "YES")
n_morning_samples <- col_data %>% group_by(Individual_ID) %>% dplyr::count() %>% data.frame()
col_data <- subset(col_data, col_data$Individual_ID %in% n_morning_samples[n_morning_samples$n > 1, "Individual_ID"])

#run analysis
count_data_temp <- count_data[as.character(gene_list$Gene_ID), as.character(col_data$Sample_ID)]
dds_counts <- DESeqDataSetFromMatrix(countData = count_data_temp, colData = col_data, design = ~ 1)
dds_counts <- estimateSizeFactors(dds_counts)
quantvst <- vst(dds_counts)
quantvst <- assay(quantvst)

gene_icc_data <- data.frame()
for (i in 1:nrow(quantvst)) {
  
  test_data <- data.frame(counts = quantvst[i,], 
                          individual = col_data$Individual_ID, 
                          timepoint = col_data$Timepoint)
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
write.table(gene_icc_data, "results/ICC_sex_differences/output/ICC_results_premenopausal_females.txt", sep = '\t', row.names = FALSE, quote = FALSE)





### Density plot for pre-and post-menopause ICC comparison as in Extended Data Fig. 5b ------------------------------------------------------------

gene_icc_data_premenopause <- read.table("results/ICC_sex_differences/output/ICC_results_premenopausal_females.txt", header = TRUE, sep = "\t")
gene_icc_data_premenopause$group = "Female (Premenopausal)"
gene_icc_data_premenopause$Category <- ifelse(gene_icc_data_premenopause$gene_type == "protein_coding", "Protein coding", "Non protein coding")

gene_icc_data_postmenopause <- read.table("results/ICC_sex_differences/output/ICC_results_postmenopausal_females.txt", header = TRUE, sep = "\t")
gene_icc_data_postmenopause$group = "Female (Postmenopausal)"
gene_icc_data_postmenopause$Category <- ifelse(gene_icc_data_postmenopause$gene_type == "protein_coding", "Protein coding", "Non protein coding")

gene_icc_data_males <- read.table("results/ICC_sex_differences/output/ICC_results_males.txt", header = TRUE, sep = "\t")
gene_icc_data_males$group = "Male"
gene_icc_data_males$Category <- ifelse(gene_icc_data_males$gene_type == "protein_coding", "Protein coding", "Non protein coding")

merged_results <- gene_icc_data_males %>%
  rbind(gene_icc_data_premenopause) %>%
  rbind(gene_icc_data_postmenopause)

labels <- merged_results %>%
  add_column("ICC_threshold" = ifelse(.$ICC < 0.5, "ICC < 0.5", "ICC > 0.5")) %>%
  group_by(group, Category, ICC_threshold) %>%
  summarize("n" = n()) %>%
  group_by(group, Category) %>%
  mutate(Percentage =  100 * n/sum(n)) %>%
  add_column("x" = rep(c(0.05, 0.85, 0.05, 0.85), 3),
             "y" = rep(c(2.4, 2.6, 2.2, 2.4, 2.0, 2.2), each = 2))

density_plot_sex <- ggplot(merged_results, aes(x = ICC, col = group)) + 
  geom_density() +
  ylab("Density (genes)") +
  facet_wrap(~ Category, scales = "free", ncol = 1, strip.position = "right") +
  geom_vline(xintercept = 0.5) +
  geom_text(data = labels, aes(x = x, y = y, label = scales::percent(Percentage, scale = 1)), hjust = 0, size = 1.7, show.legend = FALSE) +
  scale_color_manual(values = c("Male" = "#9933FF", "Female (Premenopausal)" = "#FDB064", "Female (Postmenopausal)" = "#D67009"), name = "") +
  xlim(c(0, 1)) +
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(), 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank()) 




### Boxplot for pre-and post-menopause ICC comparison as in Extended Data Fig. 5c ------------------------------------------------------------

plot_data <- merged_results %>%
  .[, c("gene_name", "Category", "group", "Var_within", "Var_between")] %>%
  melt(id.vars = c("gene_name", "Category", "group"))

wilcox_test <- plot_data %>%
  group_by(Category, variable) %>%
  wilcox_test(value ~ group, data = . , paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = rep(c(0.19, 0.18, 0.17), 4))

boxplots_sex <- ggplot(plot_data) + 
  geom_boxplot(aes(x = group, y = value, fill = group)) +
  facet_grid(Category ~ variable, labeller = labeller(variable = c("Var_within"="Within\nindividuals", "Var_between"="Between\nindividuals"))) +
  xlab("") + ylab("Variance") +
  scale_fill_manual(values = c("Male" = "#9933FF", "Female (Premenopausal)" = "#FDB064", "Female (Postmenopausal)" = "#D67009"), name = "") +
  stat_pvalue_manual(data = wilcox_test, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) + 
  scale_y_continuous(limits = quantile(plot_data$value, c(0.1, 0.9))) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_blank())



