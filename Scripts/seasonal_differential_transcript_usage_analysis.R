## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(DRIMSeq)
library(DESeq2)
library(lme4)
library(Maaslin2)
library(Gviz)
library(ape)
library(GenomicFeatures)
library(VennDiagram)
library(ggpubr)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")
source("results/enrichment_analysis_helper_functions.R")


## Loading data ------------------------------------------------------------

col_clinical_data <- read.table("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = "\t", header = TRUE)
col_clinical_data$Individual_ID <- factor(col_clinical_data$Individual_ID)
rownames(col_clinical_data) <- col_clinical_data$Sample_ID
col_clinical_data$Run <- factor(col_clinical_data$Run)

#add factor to align timepoints of sampling
col_clinical_data$Time <- ifelse(col_clinical_data$Month %in% c("January", "February", "March"), 1, 
                                 ifelse(col_clinical_data$Month %in% c("April", "May", "October", "November"), 2, 3))

#Get RUVSeq factors of unwanted variation
unwanted_factors <- read.table("results/RUVSeq/output/unwanted_factors_RUVg_2_3.txt", header = TRUE)
col_clinical_data <- col_clinical_data %>%
  merge(unwanted_factors, by = "Sample_ID")

col_clinical_data$W_2 <- cut(col_clinical_data$W_2, 5)
col_clinical_data$W_3 <- cut(col_clinical_data$W_3, 5)

# import gene and transcript annotation
tx2gene <- read.table("info_files/tx2gene.txt")
gene_attributes <- read.table("info_files/filtered_genes_samples25_fpm5.txt", sep = "\t", header = TRUE)
annotation <- gene_attributes[, c("Gene_ID", "gene_name")] %>% unique()

protein_coding <- gene_attributes %>%
  subset(.$gene_type == "protein_coding") %>%
  .$Gene_ID

protein_coding_transcripts <- tx2gene %>%
  subset(.$gene_id %in% protein_coding) %>%
  .$t_name

#import transcript counts
count_data <- read.table("count_files/merged_transcript_counts_dtuScaledTPM.txt", header = TRUE, sep = '\t', row.names = 1)
col_clinical_data <- subset(col_clinical_data, col_clinical_data$Sample_ID %in% colnames(count_data))
count_data <- count_data[, as.character(col_clinical_data$Sample_ID)]
colnames(count_data) <- col_clinical_data$Sample_ID

rownames(col_clinical_data) <- col_clinical_data$Sample_ID



## Prepare count data -------------------------------------------------

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


#normalize filtered counts and exclude non protein coding genes
rownames(filtered_counts) <- filtered_counts$feature_id
filtered_counts <- filtered_counts[, c(-1, -2)]

normalized_counts <- filtered_counts %>%
  DESeqDataSetFromMatrix(colData = col_clinical_data, design = ~ Individual_ID + W_2 + W_3 + Time) %>%
  estimateSizeFactors() %>%
  counts(normalized = TRUE) %>%
  as.data.frame()

normalized_counts <- normalized_counts[protein_coding_transcripts, ]


#compute splice ratios 
genes <- tx2gene[tx2gene$t_name %in% rownames(normalized_counts), ]$gene_id %>%
  unique()
splice_ratios <- data.frame()
g <- 1
for (g in 1:length(genes)) {
  
  gene <- genes[g]
  genename <- tx2gene[tx2gene$gene_id == gene, ]$gene_name[1]
  tx_ids <- subset(tx2gene, tx2gene$gene_id == gene)

  splice_ratios_temp <- normalized_counts %>%
    subset(rownames(.) %in% tx_ids$t_name) %>%
    as.matrix() %>%
    prop.table(margin = 2) 
  
  splice_ratios <- rbind(splice_ratios, splice_ratios_temp)
}




## Identify genes with season-dependent splicing patterns ----------------------------------------------------------
#Identify genes with DTU between annual seasons by testing the splice ratios of each transcript for an effect of season with a LMM.
#This is done for every subgroup of samples

morning_individuals <- col_clinical_data %>% 
  subset(.$morning_sample == "YES") %>%
  group_by(Individual_ID) %>%
  dplyr::summarize("n" = n()) %>%
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
  splice_ratios_temp <- splice_ratios[, my_samples]

  fit_data <- Maaslin2(input_data = splice_ratios_temp, 
                       input_metadata = table_temp, 
                       output = paste0("results/DTU/output/LMM_Maaslin2_DTU_analysis_", my_groups[i], "/"),
                       transform = "NONE",
                       fixed_effects = c("W_2", "W_3", "Time"),
                       random_effects = c("Individual_ID"),
                       normalization = 'NONE',
                       standardize = FALSE, 
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)
}





## Code for plots as in Figure 4 and Extended Data Fig. 2d ---------------------------

#get significant genes (DTU)
DTU_all <- read_tsv(paste0("results/DTU/output/LMM_Maaslin2_DTU_analysis_all/all_results.tsv")) 

DTU <- DTU_all %>%
  merge(tx2gene, by.x = "feature", by.y = "t_name") %>%
  dplyr::rename("Gene_ID" = "gene_id", "Gene_name" = "gene_name") %>%
  subset(.$metadata == "Time") 

sig_DTU <- DTU %>%
  subset(.$qval < 0.05) 
DTU_genes <- sig_DTU$Gene_ID %>% unique()

#get significantly genes (DE)
DE <- read.table("results/DGE_seasons/output/all/all_results.txt", header = TRUE) %>%
  dplyr::rename("Gene_ID" = "feature", "Gene_name" = "gene_name") 

sig_DE <- DE %>%
  subset(.$qval < 0.05) %>%
  subset(.$metadata == "Time")
DE_genes <- sig_DE$Gene_ID %>% unique()




### DE and DTU dotplot as in Figure 4d ----------------

df_DE_DTU <- DTU %>%
  .[order(.$qval, decreasing = FALSE), ] %>%
  .[match(.$Gene_ID %>% unique(), .$Gene_ID), ] %>%
  .[, c("Gene_ID", "Gene_name", "qval")] %>%
  unique() %>%
  merge(DE, by = c("Gene_ID", "Gene_name"), suffixes = c("_DTU", "_DE")) %>%
  add_column("log_qval_DTU" = -log(.$qval_DTU, base = 10),
             "log_qval_DE" = -log(.$qval_DE, base = 10),
             "group" = ifelse(.$Gene_ID %in% DE_genes & .$Gene_ID %in% DTU_genes, "both", 
                              ifelse(.$Gene_ID %in% DE_genes, "DE", 
                                     ifelse(.$Gene_ID %in% DTU_genes, "DTU", "no"))),
             "label_genes" = NA) 

#label top genes upregulated in summer
df_summer_up <- df_DE_DTU %>%
  subset(.$coef > 0) %>%
  .[order(.$qval_DE, decreasing = TRUE), ]

DTU_label <- df_summer_up %>%
  subset(.$group == "DTU") %>%
  .[order(.$qval_DTU, decreasing = FALSE), ]
DE_label <- df_summer_up %>%
  subset(.$group == "DE") %>%
  .[order(.$qval_DE, decreasing = FALSE), ]
both_label_DE <- df_summer_up %>%
  subset(.$group == "both") %>%
  .[order(.$qval_DE, decreasing = FALSE), ]
both_label_DTU <- df_summer_up %>%
  subset(.$group == "both") %>%
  .[order(.$qval_DTU, decreasing = FALSE), ]
idx <- which((df_summer_up$Gene_ID %in% head(DTU_label, 2)$Gene_ID) | 
               (df_summer_up$Gene_ID %in% head(DE_label, 2)$Gene_ID) |
               (df_summer_up$Gene_ID %in% head(both_label_DE, 2)$Gene_ID) |
               (df_summer_up$Gene_ID %in% head(both_label_DTU, 2)$Gene_ID))

df_summer_up[idx, ]$label_genes <- df_summer_up[idx, ]$Gene_name

#label top genes upregulated in winter
df_winter_up <- df_DE_DTU %>%
  subset(.$coef < 0) %>%
  .[order(.$qval_DE, decreasing = TRUE), ]

DTU_label <- df_winter_up %>%
  subset(.$group == "DTU") %>%
  .[order(.$qval_DTU, decreasing = FALSE), ]
DE_label <- df_winter_up %>%
  subset(.$group == "DE") %>%
  .[order(.$qval_DE, decreasing = FALSE), ]
both_label_DE <- df_winter_up %>%
  subset(.$group == "both") %>%
  .[order(.$qval_DE, decreasing = FALSE), ]
both_label_DTU <- df_winter_up %>%
  subset(.$group == "both") %>%
  .[order(.$qval_DTU, decreasing = FALSE), ]
idx <- which((df_winter_up$Gene_ID %in% head(DTU_label, 2)$Gene_ID) | 
               (df_winter_up$Gene_ID %in% head(DE_label, 2)$Gene_ID) |
               (df_winter_up$Gene_ID %in% head(both_label_DE, 2)$Gene_ID) |
               (df_winter_up$Gene_ID %in% head(both_label_DTU, 2)$Gene_ID))

df_winter_up[idx, ]$label_genes <- df_winter_up[idx, ]$Gene_name

#final data frame
df_volcano_DTU <- rbind(data.frame(df_summer_up, "direction" = "up"),
                        data.frame(df_winter_up, "direction" = "down")) %>%
  add_column("log_qval_gene_plot" = ifelse(.$direction == "up", .$log_qval_DTU, - .$log_qval_DTU))

df_volcano_DTU$group <- factor(df_volcano_DTU$group, levels = c("DE", "both", "DTU"), labels = c("DEG", "DEG & DTU", "DTU"))

#labels for plot
labels <- data.frame(
  Direction = c("High in winter", "High in summer"),
  Label = c(paste0("High in winter"), paste0("High in summer")),
  x = c(-6, 3),
  y = c(27, 27)
)

volcano_plot_DTU <- ggplot(df_volcano_DTU, aes(x = log_qval_gene_plot, y = log_qval_DE)) +
  geom_point(aes(color = group), pch = 20) +
  geom_text(data = labels, mapping = aes(x = x, y = y, label = Label)) + 
  geom_vline(xintercept = -log(0.05, base = 10), lty = 3) +
  geom_vline(xintercept = log(0.05, base = 10), lty = 3) +
  geom_hline(yintercept = -log(0.05, base = 10),lty = 3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_text_repel(aes(color = group, label = label_genes), min.segment.length = 0, show.legend = FALSE) +
  scale_x_continuous(breaks = c(-7.5, -5.0, -2.5, 0, 2.5, 5, 7.5, 10), labels = c("7.5", "5.0", "2.5", "0.0", "2.5", "5.0", "7.5", "10.0")) +
  scale_color_manual(values = c("#f5bf00", "#FF7F00", "#A65628"), breaks = c("DEG", "DEG & DTU", "DTU")) +
  labs(x = "Transcript usage (-Log10 Q-value)", 
       y = "Gene expression\n(-Log10 Q-value)") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())



### Functional enrichment dotplot as in Figure 4e ----------------

gene_universe <- intersect(DE$Gene_ID, DTU$Gene_ID) %>% unique()
sig_DE_sub <- DE_genes %>%
  .[. %in% gene_universe]
sig_DTU_sub <- DTU_genes %>%
  .[. %in% gene_universe]

mean_expression <- read.table("results/expr_var_over_total_var/output/mean_gene_expression.txt", header = TRUE)
rownames(mean_expression) <- mean_expression$Gene_ID
base_mean <- mean_expression %>%
  subset(rownames(.) %in% gene_universe)

direction <- c("DTU_without_DE", "DTU_with_DE")
for(d in 1:length(direction)){
  
  if(direction[d] == "DTU_without_DE") {
    degs <- setdiff(sig_DTU_sub, sig_DE_sub)
  }else if(direction[d] == "DTU_with_DE") {
    degs <- intersect(sig_DTU_sub, sig_DE_sub)
  }
  
  overallBaseMean <- as.matrix(base_mean[, "mean_expression", drop = F]) 
  topGOResults <- run_GO_ORA(degs, 
                             overallBaseMean, 
                             selected_database = "org.Hs.eg.db", 
                             back_num = 20,
                             annotation = annotation)
  
  go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
  
  file_filtered <- paste0("results/DTU/output/ORA_GO_", direction[d], ".csv")
  write.csv(go_results_filtered, file_filtered, row.names = FALSE)
}


direction <- c("DTU_without_DE", "DTU_with_DE")
GO_tables <- list()
for(d in 1:length(direction)){
  
  GO_tables <- append(GO_tables,
                      list(read_csv(paste0("results/DTU/output/ORA_GO_", direction[d], ".csv"))))
}

terms_only_DTU <- GO_tables[[1]] %>% .$GO.ID %>% .[c(1, 2, 3, 4, 7, 8, 11, 12)]
terms_DTU_and_DE <- GO_tables[[2]] %>% .$GO.ID %>% .[c(3, 5, 6, 8, 9, 18, 20)]
showTerms <- c(terms_only_DTU, terms_DTU_and_DE)

df_plot <- get_GO_df(GO_tables = GO_tables,
                     score_name = "Fisher.elim",
                     onto = "ALL",
                     direction = c("only DTU", "DTU and DEG"),
                     showTerms = showTerms,
                     multLines = TRUE,
                     numChar = 80)

GO_dotplot <- ggplot(df_plot, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  scale_colour_gradient(high="#990000", low="#FF9999", name = expression(-~log["10"]~italic(p))) + 
  scale_size(name = "Gene ratio", range = c(0, 3)) +
  theme_bw() +
  theme(axis.title.y = element_blank())



### Expression and splice ratio bar plot for TREM1 as in Figure 4f ----------------

my_colors <- c("Winter" = '#3399FF', "Spring/Autumn" = '#FFB266', "Summer" = '#FF6666')
my_gene_id <- subset(tx2gene, tx2gene$gene_name == "TREM1")[1, "gene_id"]
tx_ids <- subset(tx2gene, tx2gene$gene_id == my_gene_id)


#### Gene expression bar plot -----
#import, filter and normalize gene-level counts
count_data <- read.table("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)
col_clinical_data <- subset(col_clinical_data, col_clinical_data$Sample_ID %in% colnames(count_data))
count_data <- count_data[, as.character(col_clinical_data$Sample_ID)]
colnames(count_data) <- col_clinical_data$Sample_ID

count_data <- count_data[as.character(gene_attributes$Gene_ID %>% unique()), ]
rownames(count_data) <- gene_attributes$Gene_ID %>% unique()

vst_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ 1) %>%
  estimateSizeFactors() %>%
  vst() %>%
  assay()

#get expression of TREM1
data_trem1_expr <- vst_counts[my_gene_id, ] %>%
  data.frame("Sample_ID" = names(.), "expression" = .) %>%
  merge(col_clinical_data[, c("Sample_ID", "Time", "Individual_ID")], by = "Sample_ID")

data_trem1_expr$Time <- factor(data_trem1_expr$Time, levels = 1:3, labels = c("Winter", "Spring/Autumn", "Summer"))
data_trem1_expr$Individual_ID <- factor(data_trem1_expr$Individual_ID, levels = unique(data_trem1_expr$Individual_ID))
data_trem1_expr$group <- as.numeric(data_trem1_expr$Individual_ID)
data_trem1_expr$variable <- gene

#barplot of TREM1 expression at gene level
trem1_expression_plot <- ggplot(data_trem1_expr, aes(x = variable, y = expression, fill = Time, color = Time)) +
  geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(y = "Variance stabilized gene counts", title = gene) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "plain"),
        title = element_text(face = "italic"),
        panel.grid = element_blank(),
        legend.position = "none")


#### Splice ratio bar plot -----
#get splice ratios for TREM1
data_trem1_splice <- normalized_counts %>%
  subset(rownames(normalized_counts) %in% tx_ids$t_name) %>%
  t() %>%
  `/` (rowSums(.)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  melt(id.vars = "sample", value.name = "splice_ratio") %>%
  merge(col_clinical_data[, c("Sample_ID", "Time", "Individual_ID")], by.x = "sample", by.y = "Sample_ID")

data_trem1_splice$Time <- factor(data_trem1_splice$Time, levels = 1:3, labels = c("Winter", "Spring/Autumn", "Summer"))
data_trem1_splice$Individual_ID <- factor(data_trem1_splice$Individual_ID, levels = unique(data_trem1_splice$Individual_ID))
data_trem1_splice$group <- as.numeric(data_trem1_splice$Individual_ID)
data_trem1_splice$variable <- data_trem1_splice %>%
  group_by(variable) %>%
  dplyr::summarize("mean" = mean(splice_ratio)) %>%
  .[order(.$mean, decreasing = TRUE), ] %>%
  {factor(data_trem1_splice$variable, levels = .$variable)}

#barplot of TREM1 splice ratios
trem1_splice_plot <- ggplot(data_trem1_splice, aes(x = variable, y = splice_ratio, fill = Time, color = Time)) +
  geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(y = "Splice ratio") + 
  ylim(c(0, 1)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1.0),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


####Genomic region plot for TREM1 ------
txdb <- makeTxDbFromGFF("info_files/gencode.v25.annotation.gff3", format = "gff3")
geneTrack <- GeneRegionTrack(txdb, chromosome = "chr6", start = 41267926, end = 41286682)

ranges(geneTrack)$gene <- ranges(geneTrack)$gene %>% substr(1, 15)
ranges(geneTrack)$transcript <- ranges(geneTrack)$transcript %>% substr(1, 15)
ranges(geneTrack) <- ranges(geneTrack) %>%
  subset(.$transcript %in% c("ENST00000244709", "ENST00000334475"))

plotTracks(geneTrack, 
           chromosome = "chr6", 
           showId = TRUE, 
           transcriptAnnotation = "transcript",
           reverseStrand = TRUE,
           extend.right = 0.2,
           fontsize.group = 20,
           fontface.group = "plain",
           fontcolor.group = "black")



### Venn diagram as in Extended Data Fig. 7d ----------------

DTU_DE_colors <- c("#f5bf00", "#FF7F00", "#A65628")
my_cex <- 0.5

p_venn <- venn.diagram(x = list(DE_genes, DTU_genes), 
                       filename = NULL,
                       output = TRUE,
                       disable.logging = TRUE,
                       height = 7,
                       width = 7, 
                       resolution = 500,
                       units = "in",
                       margin = 0, 
                       compression = "lzw",
                       inverted = TRUE, 
                       cat.pos = c(145, 215),
                       col= c("black", "black"),
                       fill = DTU_DE_colors[c(1, 3)],
                       cex = c(my_cex, my_cex, my_cex),
                       fontfamily = "sans",
                       category.names = c("DEG", "DTU"),
                       cat.cex = my_cex,
                       cat.fontface = "plain",
                       cat.default.pos = "outer",
                       cat.dist = c(0.2,0.2), 
                       cat.fontfamily = "sans",
                       cat.col = c("black", "black"),
                       alpha = c(0.7, 0.7),
                       label.col = "black",
                       ext.text = FALSE, 
                       fontfamiliy = "sans")



### Scatter plots for splice ratios as in Extended Data Fig. 2e, middle right and right ----------------

#import results
DTU_list <- list()
my_groups <- c("all", "healthy", "morning")
my_groups_names <- c("All", "Healthy", "Morning")

for(g in 1:length(my_groups)){
  
  res_table <- read_tsv(paste0("results/DTU/output/LMM_Maaslin2_DTU_analysis_", my_groups[g], "/all_results.tsv")) %>%
    merge(tx2gene, by.x = "feature", by.y = "t_name") %>%
    dplyr::rename("Gene_ID" = "gene_id", "Gene_name" = "gene_name") %>%
    subset(.$metadata == "Time") 
  
  DTU_list[[g]] <- res_table
  names(DTU_list)[g] <- my_groups_names[g]
}

table_DTU <- merge(DTU_list[[2]], DTU_list[[3]], by = c("Gene_ID", "Gene_name", "feature"), suffixes = c("_healthy", "_morning")) %>%
  merge(DTU_list[[1]], by = c("Gene_ID", "Gene_name", "feature"))

df_DTU_scatter <- table_DTU %>%
  add_column("diff_all_healthy" = .$coef - .$coef_healthy,
             "diff_all_morning" = .$coef - .$coef_morning)

#healthy
p_healthy <- ggplot(df_DTU_scatter, aes(x = coef, y = coef_healthy)) +
  geom_point(color = "skyblue2", alpha = 0.5) +
  geom_vline(xintercept = 0.0) +
  geom_hline(yintercept = 0.0) +
  scale_x_continuous(limits = c(-0.055, 0.055), breaks = c(-0.05, 0.00, 0.05), labels =  c("-0.05", "0.00", "0.05")) +
  scale_y_continuous(limits = c(-0.09, 0.09), breaks = c(-0.08, 0.00, 0.08), labels =  c("-0.08", "0.00", "0.08")) +
  xlab("Splice ratio change (All samples)") + ylab("Splice ratio change (Without medical conditions)") +
  stat_cor(method = "spearman") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.x = element_blank())


#morning
p_morning <- ggplot(df_DTU_scatter, aes(x = coef, y = coef_morning)) +
  geom_point(color = "skyblue2", alpha = 0.5) +
  geom_vline(xintercept = 0.0) +
  geom_hline(yintercept = 0.0) +
  scale_x_continuous(limits = c(-0.055, 0.055), breaks = c(-0.05, 0.00, 0.05), labels =  c("-0.05", "0.00", "0.05")) +
  scale_y_continuous(limits = c(-0.09, 0.09), breaks = c(-0.08, 0.00, 0.08), labels =  c("-0.08", "0.00", "0.08")) +
  xlab("Splice ratio change (All samples)") + ylab("Splice ratio change (Morning samples)") +
  stat_cor(method = "spearman") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.title.x = element_blank())
