## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(DRIMSeq)
library(DESeq2)
library(psych)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(genefilter)
library(geneplotter)
library(topGO)
library(DBI)
library(plyr)
library(clusterProfiler)
library(msigdbr)
library(ICC)
library(rstatix)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")
source("results/enrichment_analysis_helper_functions.R")

##Loading data ------------------------------------------------------------

col_clinical_data <- read.table("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = "\t", header = TRUE)
col_clinical_data$Individual_ID <- factor(col_clinical_data$Individual_ID)
rownames(col_clinical_data) <- col_clinical_data$Sample_ID
col_clinical_data$Run <- factor(col_clinical_data$Run)

tx2gene <- read.table("info_files/tx2gene.txt")

gene_attributes <- read.table("info_files/filtered_genes_samples25_fpm5.txt", sep = "\t", header = TRUE)
annotation <- gene_attributes[, c("Gene_ID", "gene_name")] %>% unique()
protein_coding_genes <- subset(gene_attributes, gene_attributes$gene_type == "protein_coding")$Gene_ID
sex_chromosome_genes <- subset(gene_attributes, gene_attributes$Chr %in% c("chrX", "chrY"))$Gene_ID

transcript_attributes <- read.table("info_files/transcript_attributes.txt")
protein_coding_transcripts <- subset(transcript_attributes, transcript_attributes$transcript_type == "protein_coding")$transcript_id

#import transcript counts
count_data <- read.table("count_files/merged_transcript_counts_dtuScaledTPM.txt", header = TRUE, sep = '\t', row.names = 1)
col_clinical_data <- subset(col_clinical_data, col_clinical_data$Sample_ID %in% colnames(count_data))
count_data <- count_data[, as.character(col_clinical_data$Sample_ID)]
colnames(count_data) <- col_clinical_data$Sample_ID





## Filter and normalize transcript counts -------------------------------------------------

#keep only genes that passed the gene expression filter 
genes_to_keep <- gene_attributes$Gene_ID %>% unique()

transcripts_to_keep <- tx2gene %>%
  subset(.$gene_id %in% genes_to_keep) %>%
  .$t_name %>%
  unique() 

count_data <- count_data[transcripts_to_keep, ]
print(all(rownames(count_data) == transcripts_to_keep))

#filter out low expressed transcripts
n <- ncol(count_data)
n_25 <- n*0.25

count_data_drimseq <- rownames_to_column(count_data, var = "feature_id")
count_data_drimseq <- merge(tx2gene[, c("t_name", "gene_id")], count_data_drimseq, by.x = "t_name", by.y = "feature_id")

col_data_drimseq <- data.frame("sample_id" = col_clinical_data$Sample_ID, "group" = "all")
count_data_drimseq[, 3:ncol(count_data_drimseq)] <- count_data_drimseq[, col_data_drimseq$sample_id]
colnames(count_data_drimseq) <- c("feature_id", "gene_id", col_data_drimseq$sample_id)

d <- DRIMSeq::dmDSdata(counts = count_data_drimseq, samples = col_data_drimseq)

res  <-  dmFilter(d,
                  min_samps_feature_expr = n_25, min_feature_expr = 5,
                  min_samps_feature_prop = n_25, min_feature_prop = 0,
                  min_samps_gene_expr = n_25, min_gene_expr = 0)
filtered_counts <- counts(res)

#normalize filtered counts 
rownames(filtered_counts) <- filtered_counts$feature_id
filtered_counts <- filtered_counts[, c(-1, -2)]

normalized_counts <- filtered_counts %>%
  DESeqDataSetFromMatrix(colData = col_clinical_data, design = ~ 1) %>%
  estimateSizeFactors() %>%
  counts(normalized = TRUE) %>%
  as.data.frame()





## Compute contribution of alternative splicing and gene expression to the total transcript abundance variation -----------
#The method to compute the contribution of alternative splicing and gene expression to overall transcriptional variation was developed by Gonzalez-Porta et al. (http://www.genome.org/cgi/doi/10.1101/gr.121947.111)
#Lappalainen et al. extended it to link this model to within- and between-group variation (http://www.nature.com/doifinder/10.1038/nature12531)

#function to compute variation within and between groups
#X := a matrix with transcript expression values
#rows are the transcripts, columns are the samples
#group := a data.frame with one column samples and one column groups
get_var_within_between <- function(X, meta, sample_name, group_name) {
  
  #order X for sample names
  X <- X[, meta[[sample_name]]]
  colnames(X) <- meta[[sample_name]]
  
  grand_mean <- rowMeans(X)
  number_samples <- ncol(X)
  
  group <- factor(meta[[group_name]])
  group_split <- split(as.data.frame(t(X)), group)
  group_split <- sapply(group_split, t)
  group_split <- sapply(group_split, as.data.frame)
  group_length <- sapply(group_split, ncol)
  group_mean <- sapply(group_split, rowMeans)
  
  #Compute variation within groups
  var_within <- sapply(X = names(group_split), FUN = function(i)(sum((group_split[[i]] - group_mean[, i])^2)))
  var_within <- sum(var_within)  / (number_samples - 1) 
  
  #Compute variation between groups
  var_between <- as.data.frame(t(t((group_mean - grand_mean)^2) * group_length))
  var_between <- sum(var_between)  / (number_samples - 1) 
  
  return(c("var_within" = var_within, "var_between" = var_between))
}


#compute the statistic for every gene and all subsets of samples
morning_individuals <- col_clinical_data %>% 
  subset(.$morning_sample == "YES") %>%
  group_by(Individual_ID) %>%
  summarize("n" = n()) %>%
  subset(.$n > 1) %>%
  .$Individual_ID

morning_samples <- col_clinical_data %>% 
  subset(.$morning_sample == "YES") %>%
  subset(.$Individual_ID %in% morning_individuals) %>%
  .$Sample_ID

genes <- tx2gene[tx2gene$t_name %in% rownames(normalized_counts), ]$gene_id %>%
  unique()

my_groups <- c("all", "healthy", "morning")
for(i in 1:length(my_groups)){
  
  if(my_groups[i] == "all"){
    table_temp <- col_clinical_data
  }else if(my_groups[i] == "healthy"){
    table_temp <- col_clinical_data %>% subset(.$healthy == "YES")
  }else if(my_groups[i] == "morning"){
    table_temp <- col_clinical_data %>% subset(.$Sample_ID %in% morning_samples)
  }
  
  my_samples <- table_temp$Sample_ID
  print(paste0("Group: ", my_groups[i], ", number of samples: ", length(my_samples)))
  normalized_counts_temp <- normalized_counts[, my_samples]
  
  res_table <- data.frame()
  for(g in 1:length(genes)) {
    
    #transcript expression for all samples
    gene <- genes[g]
    tx_ids <- subset(tx2gene, tx2gene$gene_id == gene)
    X <- subset(normalized_counts_temp, rownames(normalized_counts_temp) %in% tx_ids$t_name)
    
    #singular value decomposition
    sing_val_dec <- svd(x = X, nu = 1, nv = 1)
    u_1 <- sing_val_dec$u
    v_1 <- sing_val_dec$v
    d_1 <- sing_val_dec$d[1]
    
    #transcript expression for samples projected on the least squares line
    X_line <- d_1 * (u_1 %*% t(v_1))
    colnames(X_line) <- colnames(X)
    rownames(X_line) <- rownames(X)
    var_line <- apply(X_line, MARGIN = 1, var) %>% sum()
    
    #get variation within and between individuals
    var_within_between_all <- get_var_within_between(X, meta = as.data.frame(table_temp[, c("Sample_ID", "Individual_ID")]), sample_name = "Sample_ID", group_name = "Individual_ID")
    var_within_between_line <- get_var_within_between(X_line, meta = as.data.frame(table_temp[, c("Sample_ID", "Individual_ID")]), sample_name = "Sample_ID", group_name = "Individual_ID")
    
    #compute the total variation in transcript expression (same as var_within_between_all[1] + var_within_between_all[2])
    sigma_mat <- cov(t(X))
    V_t <- tr(sigma_mat)
    
    #compute variation on the least squares line (same as var_within_between_line[1] + var_within_between_line[2])
    V_ls <- t(u_1) %*% sigma_mat %*% u_1
    
    res_table <- rbind(res_table,
                       data.frame("Gene_ID" = gene, 
                                  "number_of_transcripts" = nrow(X), 
                                  "total_var" = V_t, 
                                  "expr_var" = V_ls,
                                  "expr_var_within" = var_within_between_line["var_within"],
                                  "expr_var_between" = var_within_between_line["var_between"],
                                  "total_var_within" = var_within_between_all["var_within"],
                                  "total_var_between" = var_within_between_all["var_between"],
                                  "expr_var/total_var" =  V_ls / V_t,
                                  "expr_var_within/total_var_within" = var_within_between_line["var_within"] / var_within_between_all["var_within"],
                                  "expr_var_between/total_var_between" = var_within_between_line["var_between"] / var_within_between_all["var_between"],
                                  "total_var_between/total_var" = var_within_between_all["var_between"]/ V_t,
                                  "group" = my_groups[i]))
    
  }
  
  write.table(res_table, paste0("results/expr_var_over_total_var/output/expression_over_variation_statistics_", my_groups[i], ".txt"), sep = '\t', quote = FALSE, row.names = FALSE)
  
}



## Code for distribution plots of the contribution of alternative splicing and gene expression to the total variation as in Fig. 3 and Extended Data Fig. 6 ---------

res_expr_vs_splicing <- read.table("results/expr_var_over_total_var/output/expression_over_variation_statistics_all.txt", sep = '\t', header = TRUE)


#Plot density as in Figure 3a
df_expr_vs_spli <- res_expr_vs_splicing %>%
  add_column("protein" = {.$Gene_ID %in% protein_coding_genes} %>% 
               factor(levels = c(FALSE, TRUE), labels = c("Non protein coding", "Protein coding")))

expr_vs_splicing_plot <- ggplot(df_expr_vs_spli, aes(x = expr_var.total_var, color = protein)) +
  geom_density() + 
  labs(x = "Expression variation /\nTotal variation",
       y = "Density (genes)") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_color_manual(values = c("#00CCCC", "#FF9999")) +
  guides(color = guide_legend(nrow = 2)) + 
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())


#Plot density across different inter-individual variation levels as in Figure 3c
df_inter_var <- res_expr_vs_splicing 
breaks <-  summary(df_inter_var$total_var_between.total_var)
df_inter_var$Group <- cut(df_inter_var$total_var_between.total_var, 
                          breaks = c(breaks["Min."], breaks["1st Qu."], breaks["3rd Qu."], breaks["Max."]), 
                          right = TRUE, 
                          include.lowest = TRUE)
df_inter_var$Group_name <- factor(df_inter_var$Group, labels = c("low", "mid", "high"))

inter_ind_var_plot <- ggplot(df_inter_var, aes(x = expr_var.total_var, color = Group_name)) +
  geom_density() + 
  labs(x = "Expression variation /\nTotal variation",
       y = "Density (genes)",
       color = "Inter-individual variation /\nTotal variation") +
  scale_color_manual(values = c("low" = "#0072B2", "mid" = "#E69F00", "high" = "#CC79A7"), breaks = c("low", "mid", "high")) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())


#Plot histogram colored by mean gene expression as in Extended Data Fig. 6a
my_colors <- c("high" = "#CC79A7", "mid" = "#E69F00", "low" = "#0072B2")
mean_expression <- read.table("results/expr_var_over_total_var/output/mean_gene_expression.txt", header = TRUE)

df_expr_vs_spli <- res_expr_vs_splicing %>%
  merge(mean_expression[, c(1, 3)], by = "Gene_ID") %>%
  .[order(.$mean_expression, decreasing = FALSE), ] %>%
  add_column("expr_group" = cut_number(.$mean_expression, 3))

df_expr_vs_spli$expr_group <- factor(df_expr_vs_spli$expr_group, labels = c("low", "mid", "high") %>% rev())

avg_expression_plot <- ggplot(df_expr_vs_spli, aes(x = expr_var.total_var)) +
  geom_histogram(binwidth = 0.025, boundary = 0, color = "black", aes(fill = expr_group), alpha = 0.7) + 
  scale_fill_manual(values = my_colors, breaks = c("low", "mid", "high")) +
  labs(x = "Expression variation /\nTotal variation",
       y = "Number of genes",
       fill = "Mean expression") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(title.position = "top")) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())


#Plot density as in Extended Data Fig. 6b
my_colors <- c("expr_var_within.total_var_within" = "#CCCC00", "expr_var_between.total_var_between" = "#009900", "expr_var.total_var" = "#orange")

df_intra_inter <- res_expr_vs_splicing %>%
  .[, c(1, 10:11)] %>%
  melt(id.vars = "Gene_ID")

inter_intra_plot <- ggplot(df_intra_inter, aes(x = value, color = variable)) +
  geom_density() +  
  labs(y = "Genes (density)", color = "Expression variation/\nTotal variation") +  
  scale_color_manual(labels = c("expr_var_within.total_var_within" = "Within individuals", 
    "expr_var_between.total_var_between" = "Between individuals"), 
    values = my_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  guides(color = guide_legend(title.position = "top")) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())


#Plot scatterplots as in Extended Data Fig. 6d
table_all <- res_expr_vs_splicing

#healthy
table_healthy <- read.table("results/expr_var_over_total_var/output/expression_over_variation_statistics_healthy.txt", sep = '\t', header = TRUE)

df_plot <- merge(table_all, table_healthy,  by = c("Gene_ID"), suffixes = c("_all", "_healthy"))
df_plot <- df_plot %>%
  merge(annotation[, c("Gene_ID", "gene_name")], by  ="Gene_ID") %>%
  add_column("difference" = .$expr_var.total_var_all - .$expr_var.total_var_healthy) %>%
  add_column("label_genes" = NA)

genes_low <- df_plot %>%
  .[order(.$difference, decreasing = FALSE), ]
genes_high <- df_plot %>%
  .[order(.$difference, decreasing = TRUE), ]
idx <- which((df_plot$Gene_ID %in% head(genes_low, 10)$Gene_ID) | 
               (df_plot$Gene_ID %in% head(genes_high, 10)$Gene_ID))
df_plot[idx, ]$label_genes <- df_plot[idx, ]$gene_name

p_healthy <- ggplot(df_plot, aes(x = expr_var.total_var_all, y = expr_var.total_var_healthy)) +
  geom_point(color = "skyblue2", alpha = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_text_repel(aes(label = label_genes), min.segment.length = 0) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels =  c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(limits = c(0, 1), breaks =  c(0, 0.25, 0.5, 0.75, 1), labels =  c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  labs(x = "Expression variation/Total variation (All)", y = "Expression variation/Total variation (Healthy)") +
  stat_cor(method = "spearman") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  coord_equal()

#morning
table_morning <- read.table("results/expr_var_over_total_var/output/expression_over_variation_statistics_morning.txt", sep = '\t', header = TRUE)

df_plot <- merge(table_all, table_morning,  by = c("Gene_ID"), suffixes = c("_all", "_morning"))
df_plot <- df_plot %>%
  merge(annotation[, c("Gene_ID", "gene_name")], by  ="Gene_ID") %>%
  add_column("difference" = .$expr_var.total_var_all - .$expr_var.total_var_morning) %>%
  add_column("label_genes" = NA)

genes_low <- df_plot %>%
  .[order(.$difference, decreasing = FALSE), ]
genes_high <- df_plot %>%
  .[order(.$difference, decreasing = TRUE), ]
idx <- which((df_plot$Gene_ID %in% head(genes_low, 10)$Gene_ID) | 
               (df_plot$Gene_ID %in% head(genes_high, 10)$Gene_ID))
df_plot[idx, ]$label_genes <- df_plot[idx, ]$gene_name

p_morning <- ggplot(df_plot, aes(x = expr_var.total_var_all, y = expr_var.total_var_morning)) +
  geom_point(color = "skyblue2", alpha = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_text_repel(aes(label = label_genes), min.segment.length = 0) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels =  c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(limits = c(0, 1), breaks =  c(0, 0.25, 0.5, 0.75, 1), labels =  c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  labs(x = "Expression variation/Total variation (All)", y = "Expression variation/Total variation (Morning)") +
  stat_cor(method = "spearman") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  coord_equal()



## Code for GO dot plot as in Figure 3b ---------

#GO enrichment analysis
res_expr_vs_splicing <- read.table("results/expr_var_over_total_var/output/expression_over_variation_statistics_all.txt", sep = '\t', header = TRUE)
res_table_protein <- res_expr_vs_splicing %>%
  subset(.$Gene_ID %in% protein_coding_genes)
gene_universe <- res_table_protein$Gene_ID %>% unique()

mean_expression <- read.table("results/expr_var_over_total_var/output/mean_gene_expression.txt", header = TRUE)
rownames(mean_expression) <- mean_expression$Gene_ID
base_mean <- mean_expression %>%
  subset(rownames(.) %in% gene_universe)

direction <- c("expr_var.total_var_low", "expr_var.total_var_high")
for(d in 1:length(direction)){
  
  idx <- order(res_table_protein$expr_var.total_var, decreasing = FALSE)
  data <- res_table_protein[idx, ]
  
  if(direction[d] == "expr_var.total_var_low") {
    degs <- data[1:1000, ]$Gene_ID %>% unique() 
  }else if(direction[d] == "expr_var.total_var_high") {
    degs <- data[(nrow(data)-999):nrow(data), ]$Gene_ID %>% unique() 
  }
  
  overallBaseMean <- as.matrix(base_mean[, "mean_expression", drop = F])
  
  topGOResults <- run_GO_ORA(degs, 
                             overallBaseMean, 
                             selected_database = "org.Hs.eg.db", 
                             back_num = 10,
                             annotation = annotation)
  
  go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
  
  file_filtered <- paste0("results/expr_var_over_total_var/output/ORA_GO_", direction[d], ".csv")
  write.csv(go_results_filtered, file_filtered, row.names = FALSE)
}

#plot GO results
GO_tables <- list()
d <- 1
for(d in 1:length(direction)){
  
  GO_tables <- append(GO_tables,
                      list(read_csv(paste0("results/expr_var_over_total_var/output/ORA_GO_", direction[d], ".csv"))))
  
}

my_terms_1 <- GO_tables[[1]][c(1, 2, 3, 4, 5, 6, 8, 11), ]$GO.ID
my_terms_2 <- GO_tables[[2]][c(1, 2, 3, 4, 6, 11, 15, 17), ]$GO.ID
my_terms <- c(my_terms_1, my_terms_2)

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = c("expr_var.total_var_low", "expr_var.total_var_high"),
                     score_name = score_name,
                     reverse_bool = TRUE,
                     showTerms = my_terms,
                     multLines = TRUE,
                     numChar = 80)

GO_dotplot <- ggplot(df_plot, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  xlab("Expression variation /\n Total variation") +
  scale_x_discrete(labels = c("low", "high")) +
  scale_colour_gradient(high="#990000", low="#FF9999", name = expression(-~log["10"]~italic(p))) + 
  scale_size(name = "Gene ratio", range = c(0, 3)) +
  theme_bw() +
  theme(axis.title.y = element_blank())





## Code for barplot of ICC of splice ratios and transcript expression as in Extended Data Fig. 6c ---------

morning_individuals <- col_clinical_data %>% 
  subset(.$morning_sample == "YES") %>%
  group_by(Individual_ID) %>%
  summarize("n" = n()) %>%
  subset(.$n > 1) %>%
  .$Individual_ID

morning_samples <- col_clinical_data %>% 
  subset(.$morning_sample == "YES") %>%
  subset(.$Individual_ID %in% morning_individuals) %>%
  .$Sample_ID

genes <- tx2gene[tx2gene$t_name %in% rownames(normalized_counts), ]$gene_id %>%
  unique()

my_groups <- c("all", "healthy", "morning")

#compute ICC of splice ratios
for(i in 1:length(my_groups)){
  
  if(my_groups[i] == "all"){
    table_temp <- col_clinical_data
  }else if(my_groups[i] == "healthy"){
    table_temp <- col_clinical_data %>% subset(.$healthy == "YES")
  }else if(my_groups[i] == "morning"){
    table_temp <- col_clinical_data %>% subset(.$Sample_ID %in% morning_samples)
  }
  
  my_samples <- table_temp$Sample_ID
  print(paste0("Group: ", my_groups[i], ", number of samples: ", length(my_samples)))
  normalized_counts_temp <- normalized_counts[, my_samples]
  
  res_table <- data.frame()
  g <- 1
  for (g in 1:length(genes)) {
    
    gene <- genes[g]
    genename <- tx2gene[tx2gene$gene_id == gene, ]$gene_name[1]
    tx_ids <- subset(tx2gene, tx2gene$gene_id == gene)

    splice_ratios <- normalized_counts_temp %>%
      subset(rownames(.) %in% tx_ids$t_name) %>%
      as.matrix() %>%
      prop.table(margin = 2) %>%
      t() 
    
    df_meta <- splice_ratios %>%
      as.data.frame() %>%
      add_column("Sample_ID" = rownames(.), .before = 1) %>%
      merge(table_temp[, c("Sample_ID", "Individual_ID")], by = "Sample_ID") %>%
      na.omit()
    
    for(t in 1:ncol(splice_ratios)){
      
      transcript <- colnames(splice_ratios)[t]
      transcript_data <- data.frame(df_meta[, c("Sample_ID", "Individual_ID", transcript)])
      colnames(transcript_data) <- c("Sample_ID", "Individual_ID", "Splice_ratio")
      
      res <- ICCest(transcript_data$Individual_ID, transcript_data$Splice_ratio)
      
      res_table <- rbind(res_table, 
                         data.frame("Transcript_ID" = transcript, 
                                    "Gene_ID" = gene, 
                                    "Gene_name" = genename,
                                    "ICC" = res$ICC,
                                    "CI_lower" = res$LowerCI, 
                                    "CI_upper" = res$UpperCI, 
                                    "Var_within" = res$varw, 
                                    "Var_between" = res$vara,
                                    "expressed_in_samples" = nrow(df_meta),
                                    "number_of_isoforms" = ncol(splice_ratios),
                                    "group" = my_groups[i]))
    }
  }
  write.table(res_table, paste0("results/expr_var_over_total_var/output/ICC_transcript_ratio_", my_groups[i], ".txt"), sep = '\t', quote = FALSE, row.names = FALSE)
}

#compute ICC of transcript expression
for(i in 1:length(my_groups)){
  
  if(my_groups[i] == "all"){
    table_temp <- col_clinical_data
  }else if(my_groups[i] == "healthy"){
    table_temp <- col_clinical_data %>% subset(.$healthy == "YES")
  }else if(my_groups[i] == "morning"){
    table_temp <- col_clinical_data %>% subset(.$Sample_ID %in% morning_samples)
  }
  
  my_samples <- table_temp$Sample_ID
  print(paste0("Group: ", my_groups[i], ", number of samples: ", length(my_samples)))
  normalized_counts_temp <- normalized_counts[, my_samples]
  
  res_table <- data.frame()
  g <- 1
  for (g in 1:length(genes)) {
    
    gene <- genes[g]
    genename <- tx2gene[tx2gene$gene_id == gene, ]$gene_name[1]
    tx_ids <- subset(tx2gene, tx2gene$gene_id == gene)
    
    transcript_expression <- normalized_counts_temp %>%
      subset(rownames(.) %in% tx_ids$t_name) %>%
      t() 
    
    df_meta <- transcript_expression %>%
      as.data.frame() %>%
      add_column("Sample_ID" = rownames(.), .before = 1) %>%
      merge(table_temp[, c("Sample_ID", "Individual_ID")], by = "Sample_ID") %>%
      na.omit()
    
    for(t in 1:ncol(transcript_expression)){
      
      transcript <- colnames(transcript_expression)[t]
      transcript_data <- data.frame(df_meta[, c("Sample_ID", "Individual_ID", transcript)])
      colnames(transcript_data) <- c("Sample_ID", "Individual_ID", "Expression")
      
      expressed_in_samples <- transcript_data %>%
        .[.$Expression != 0, ] %>%
        nrow()
      
      res <- ICCest(transcript_data$Individual_ID, transcript_data$Expression)
      
      res_table <- rbind(res_table, 
                         data.frame("Transcript_ID" = transcript, 
                                    "Gene_ID" = gene, 
                                    "Gene_name" = genename,
                                    "ICC" = res$ICC,
                                    "CI_lower" = res$LowerCI, 
                                    "CI_upper" = res$UpperCI, 
                                    "Var_within" = res$varw, 
                                    "Var_between" = res$vara,
                                    "expressed_in_samples" = expressed_in_samples,
                                    "number_of_isoforms" = ncol(transcript_expression),
                                    "group" = my_groups[i]))
    }
  }
  write.table(res_table, paste0("results/expr_var_over_total_var/output/ICC_transcript_expression_", my_groups[i], ".txt"), sep = '\t', quote = FALSE, row.names = FALSE)
}


#plot ICC barplots
my_colors <- c("Var_within" = "#CCCC00", "Var_between" = "#009900")

res_splice_ratio <- read.table("results/expr_var_over_total_var/output/ICC_transcript_ratio_all.txt", sep = '\t', header = TRUE)
res_transcript_expr <- read.table("results/expr_var_over_total_var/output/ICC_transcript_expression_all.txt", sep = '\t', header = TRUE)

df_combined <- res_splice_ratio %>%
  add_column("Category" = "splice_ratio") %>%
  rbind(res_transcript_expr %>% add_column("Category" = "transcript_expr")) %>%
  subset(!(.$Gene_ID %in% sex_chromosome_genes)) %>%
  add_column("coding" = ifelse(.$Transcript_ID %in% protein_coding_transcripts, "Protein coding", "Non protein coding")) %>%
  .[, c("Transcript_ID", "Var_within", "Var_between", "Category", "coding")] %>%
  melt(id.vars = c("Transcript_ID", "Category", "coding"))

df_combined_subset <- df_combined %>%
  split(~ .$Category) %>%
  lapply(function(x){res <- x %>%
    subset(.$value > quantile(.$value, 0.1) & .$value < quantile(.$value, 0.9))
  return(res)}) %>%
  do.call(rbind, .)

df_test <- df_combined %>%
  group_by(Category, coding) %>%
  wilcox_test(value ~ variable, data = . , paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = c(0.028, 15000) %>% rep(each = 2))

var_names <- c("splice_ratio" = "Splice ratio",
               "transcript_expr" = "Transcript expression",
               "Protein coding" = "Protein coding",
               "Non protein coding" = "Non protein coding")

icc_plot <- ggplot(df_combined_subset) + 
  geom_boxplot(aes(x = variable, y = value, fill = variable)) +
  stat_pvalue_manual(data = df_test, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
  facet_grid(Category ~ coding, scales = "free_y", labeller = as_labeller(var_names)) +
  labs(y = "Variance") +
  scale_fill_manual(values = my_colors, labels = c("Var_within"="Within\nindividuals", "Var_between"="Between\nindividuals")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line())
