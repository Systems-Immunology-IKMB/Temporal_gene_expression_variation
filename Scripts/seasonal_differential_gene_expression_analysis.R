## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(DESeq2)
library(magrittr)
library(compiler)
library(stringr)
library(tidyverse)
library(Maaslin2)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


gene_list <- read.table("info_files/results_Neha/filtered_genes_samples25_fpm5.txt", sep = "\t", header = TRUE)


## Loading data ------------------------------------------------------------

gene_list <- read.csv("info_files/filtered_genes_samples25_fpm5.txt", sep = '\t')
gene_list <- subset(gene_list, gene_list$gene_type == "protein_coding")

count_data <- read.csv("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)
col_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')
unwanted_factors <- read.csv("results/RUVSeq/output/unwanted_factors_RUVg_2_3.txt", sep = '\t')

col_data <- inner_join(col_data, unwanted_factors, by = "Sample_ID")

selected_col_clinical_data <- col_data
selected_col_clinical_data$W_2 <- cut(selected_col_clinical_data$W_2, 5)
selected_col_clinical_data$W_3 <- cut(selected_col_clinical_data$W_3, 5)
selected_col_clinical_data$Time <- ifelse(selected_col_clinical_data$Month %in% c("January", "February", "March"), 1, 
                                          ifelse(selected_col_clinical_data$Month %in% c("April", "May", "October", "November"), 2, 3))

count_data <- count_data[as.character(gene_list$Gene_ID), as.character(selected_col_clinical_data$Sample_ID)]



## Identify genes with season-dependent expression patterns ----------------------------------------------------------
#Identify genes with differential gene expression between annual seasons by testing the expression of each for an effect of season with a LMM.
#This is done for every subgroup of samples

morning_individuals <- selected_col_clinical_data %>% 
  subset(.$morning_sample == "YES") %>%
  group_by(Individual_ID) %>%
  dplyr::summarize("n" = n()) %>%
  subset(.$n > 1) %>%
  .$Individual_ID

morning_samples <- selected_col_clinical_data %>% 
  subset(.$morning_sample == "YES") %>%
  subset(.$Individual_ID %in% morning_individuals) %>%
  .$Sample_ID

my_groups <- c("all", "healthy", "morning")

for(i in 1:length(my_groups)){
  
  if(my_groups[i] == "all"){
    table_temp <- selected_col_clinical_data
  }else if(my_groups[i] == "healthy"){
    table_temp <- selected_col_clinical_data %>% subset(.$healthy == "YES")
  }else if(my_groups[i] == "morning"){
    table_temp <- selected_col_clinical_data %>% subset(.$Sample_ID %in% morning_samples)
  }
  
  count_data_temp <- count_data[, as.character(table_temp$Sample_ID)]
  dds_counts <- DESeqDataSetFromMatrix(countData = count_data_temp, colData = table_temp, design = ~ Individual_ID + W_2 + W_3 + Time)
  dds_counts <- estimateSizeFactors(dds_counts)
  vst_counts <- assay(vst(dds_counts))
  rownames(table_temp) <- table_temp$Sample_ID
  
  fit_data <- Maaslin2(vst_counts, 
                       selected_col_clinical_data, 
                       output = paste0("results/DGE_seasons/output/", my_groups[i], "/"), 
                       transform = "NONE",
                       fixed_effects = c("W_2", "W_3", "Time"),
                       random_effects = c("Individual_ID"),
                       normalization = 'NONE',
                       standardize = FALSE, 
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)
}




## Volcanoplot as in Fig. 4b ---------------------------

lmm_result <- read.csv("results/DGE_seasons/output/all/all_results.txt", sep = '\t')
sig_lmm_result <- subset(lmm_result, lmm_result$qval < 0.05)
sig_lmm_result <- sig_lmm_result[order(sig_lmm_result$qval, decreasing = FALSE), ]

high_in_winter <- subset(sig_lmm_result, sig_lmm_result$coef < 0)
high_in_summer <- subset(sig_lmm_result, sig_lmm_result$coef > 0)

sig_lmm_result$Direction <- ifelse(sig_lmm_result$coef > 0, "High in summer", "High in winter")
sig_lmm_result$log_padj <- -1 * log10(sig_lmm_result$qval)

labels <- data.frame(
  Direction = c("High in winter", "High in summer"),
  Label = c(paste0("High in winter (", nrow(high_in_winter), ")"), paste0("High in summer (", nrow(high_in_summer), ")")),
  x = c(-0.15, 0.05),
  y = c(500, 500)
)

p <- ggplot(data = sig_lmm_result, aes(x=coef, y=log_padj, col=Direction)) + geom_point(pch=20) 
p <- p + geom_text_repel(data = sig_lmm_result[1:30,], aes(x=coef, y=log_padj, label=gene_name), min.segment.length = 0)
p <- p + geom_vline(xintercept = 0, lty=3)
p <- p + scale_color_manual(values = c("High in winter"='#3399FF',"High in summer"='#FF6666'))
p <- p + xlab("Expression (Beta coefficient)") + ylab("- Log10 Q-value")
p <- p + theme_bw()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               panel.border = element_blank(),
               axis.line = element_line())

dens1 <- ggplot(sig_lmm_result, aes(x = coef, fill = Direction)) + 
  geom_histogram(alpha = 0.4, col="black", bins=50) + 
  scale_fill_manual(values = c("High in winter"='#3399FF',"High in summer"='#FF6666')) +
  geom_text(data = labels, mapping = aes(x=x, y=y, label=Label, col = Direction), hjust = 0) + 
  scale_color_manual(values = c("High in winter"='#3399FF',"High in summer"='#FF6666')) +
  theme_void() + 
  theme(legend.position = "none",
        plot.tag = element_text())

volcano_plot <- dens1 + plot_spacer() + p +
  plot_layout(ncol = 2, nrow = 2, widths = c(5, 0), heights = c(1.5, 3.5))



## Average GO term heatmap as in Fig. 4c ---------------------------

heatmap_matrix <- read.csv("results/DGE_seasons/output/seasonal_degs_GO_average_heatmap_matrix.txt", sep = '\t', row.names = NULL)
heatmap_rownames <- heatmap_matrix$row.names
heatmap_matrix$row.names <- NULL
heatmap_matrix <- as.matrix(heatmap_matrix)
rownames(heatmap_matrix) <- heatmap_rownames

go_genes <- read.csv("results/DGE_seasons/output/seasonal_degs_GO_genes.txt", sep = '\t')
go_genes$Term <- str_wrap(go_genes$Term, width = 20)
rownames(sig_lmm_result) <- sig_lmm_result$feature
go_genes_beta <- sig_lmm_result[as.character(go_genes$Gene_ID), "coef"]
matDataHeat <- scale(t(heatmap_matrix))

row_labels <- heatmap_rownames[abs(go_genes_beta) > 0.03]

row_ha = rowAnnotation(df = data.frame(Timepoints = rownames(matDataHeat)),
                       col = list(Timepoints = c('1' = '#3399FF', '2' = '#FFB266', '3' = '#FF6666')), 
                       annotation_legend_param = list("Timepoints" = 
                                                        list(nrow = 1,
                                                             title = "Season",
                                                             at = c("1", "2", "3"),
                                                             labels = c("Winter", "Spring/Autumn", "Summer"),
                                                             title_position = "leftcenter")),
                       gp = gpar(col = "black"),
                       border = TRUE,
                       simple_anno_size = unit(0.2, "cm"),
                       show_annotation_name = FALSE,
                       name = "")

col_fun = colorRamp2(seq(from = -0.08, to = 0.08, by = 0.02), rev(brewer.pal(9, "Greens")))
col_ha = HeatmapAnnotation(Beta = go_genes_beta, 
                           col=list(Beta = col_fun),
                           annotation_legend_param = list("Beta" = 
                                                            list(direction = "horizontal",
                                                                 title_position = "lefttop")),
                           
                           annotation_name_rot = 0,
                           annotation_name_side = "left")
col_ha2 <- HeatmapAnnotation(gene_mark = anno_mark(at = which(abs(go_genes_beta) > 0.03), 
                                                   labels = row_labels, 
                                                   padding = unit(0.2, "mm"),
                                                   which = "column", 
                                                   side = "bottom", 
                                                   link_width = unit(0.35, "cm")))

DEG_heatmap <- Heatmap(matDataHeat, 
                       col = rev(brewer.pal(11, "RdYlBu")), 
                       top_annotation = col_ha, 
                       left_annotation = row_ha, 
                       cluster_rows = FALSE, 
                       show_row_names = FALSE,
                       column_split = go_genes$Term, 
                       column_title_rot = 0, 
                       column_gap = unit(0.7, "mm"), 
                       border = TRUE, 
                       heatmap_legend_param = list(direction = "horizontal",
                                                   title_position = "lefttop"),
                       column_dend_height = unit(0.5, "cm"),
                       cluster_column_slices = FALSE, 
                       show_column_names = TRUE,
                       name = "Scaled normalized expression")



## GO dotplot as in Extended Data Fig. 7a ---------------------------

go_plot_data <- read.csv("results/DGE_seasons/output/Winter_Summer_GO_plot_data.txt", sep = '\t')

go_plot_data$Term2 <- str_wrap(go_plot_data$Term, width = 57)
go_plot_data$Term3 <- reorder(go_plot_data$Term2, c(11:20, 1:10, 31:40, 21:30) %>% rev())
go_plot_data$Group <- factor(go_plot_data$Group, levels = c("Winter", "Summer"))

go_plot_Leuven <- ggplot(go_plot_data, mapping = aes(x = Group, y = Term3, size = n, color = p)) + 
  geom_point() + 
  xlab("") + ylab("") +
  scale_colour_gradient(high="#990000", low="#FF9999", name = expression(-~log["10"]~italic(p))) +
  scale_size(name = "Gene\nratio", range = c(0, 3)) +
  guides(color = guide_colorbar(order = 1), 
         size = guide_legend(order = 2)) +
  theme_bw() +
  theme(axis.text.y = element_text(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"), 
        legend.position = "right",
        legend.direction = "vertical")




## GO term heatmap as in Extended Data Fig. 7b ---------------------------

#load cell type specificity scores
df_cell_type_specificity <- read.table("redults/cell_type_specificity/output/tau_lowRes_celltypes.txt", sep = "\t", header = TRUE) %>%
  dplyr::rename("gene_name" = "Gene_ID") %>%
  merge(annotation, by = "gene_name")
my_breaks <- seq(from = 0, to = 1, by = 0.25)
df_cell_type_specificity$cell_type_specificity_binned <- df_cell_type_specificity$cell_type_specificity %>%
  cut(breaks = my_breaks, include.lowest = TRUE, right = TRUE)

#load heatmap data
heatmap_matrix <- read.csv("results/DGE_seasons/output/seasonal_degs_GO_heatmap_matrix.txt", sep = '\t', row.names = NULL)
heatmap_rownames <- heatmap_matrix$row.names
heatmap_matrix$row.names <- NULL
heatmap_matrix <- as.matrix(heatmap_matrix)
rownames(heatmap_matrix) <- heatmap_rownames

lmm_result <- read.csv("results/DGE_seasons/output/all/all_results.txt", sep = '\t')
sig_lmm_result <- subset(lmm_result, lmm_result$qval < 0.05)
sig_lmm_result <- sig_lmm_result[order(sig_lmm_result$qval, decreasing = FALSE), ]

go_genes <- read.csv("results/DGE_seasons/output/seasonal_degs_GO_genes.txt", sep = '\t')
go_genes$Term <- str_wrap(go_genes$Term, width = 17)
rownames(sig_lmm_result) <- sig_lmm_result$feature
go_genes_beta <- sig_lmm_result[as.character(go_genes$Gene_ID), "coef"]
matDataHeat <- t(scale(t(heatmap_matrix)))
rownames(matDataHeat) <- go_genes$Gene_name

icc_res <- read.csv("results/ICC/output/gene_icc_data_whole_dataset.txt", sep = '\t')
go_genes_icc <- left_join(go_genes, icc_res, by="Gene_ID")
col_data_temp<- selected_col_clinical_data
col_data_temp$Annual_season <- ifelse(col_data_temp$Time == '1', "Winter", 
                                 ifelse(col_data_temp$Time == '2', "Spring/Autumn", "Summer"))
col_data_temp$Annual_season <- factor(col_data_temp$Annual_season, levels = c("Winter", "Spring/Autumn", "Summer"))

bar_color <- ifelse(go_genes_icc$ICC > 0.5, '#009900', '#CCCC00')
col_fun = colorRamp2(seq(from = -0.08, to = 0.08, by = 0.02) ,rev(brewer.pal(9, "Greens")))
my_ht_colors <- structure(c("#FEE0D2", "#FCBBA1", "#FB6A4A", "#A50F15"), names = c("[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]"))

row_ha1 <- rowAnnotation("ICC" = anno_barplot(go_genes_icc$ICC, 
                                              which="row", 
                                              gp = gpar(fill=bar_color), 
                                              width = unit(0.5, "cm")), 
                         annotation_name_rot = 0,
                         annotation_name_side = "top")

row_anno <- data.frame("Cell type\nspecificity" = df_cell_type_specificity %>%
                         .[match(rownames(matDataHeat), .$gene_name), ] %>%
                         .$cell_type_specificity_binned,
                       "Beta" = go_genes_beta,
                       check.names = FALSE)

row_ha2 <- rowAnnotation(df = row_anno,
                         col = list("Cell type\nspecificity" = my_ht_colors,
                                    "Beta" = col_fun),
                         annotation_legend_param = list("Beta" = 
                                                          list(direction = "vertical",
                                                               title_position = "topleft"),
                                                        "Cell type\nspecificity" = 
                                                          list(direction = "vertical",
                                                               title_position = "topleft")),
                         annotation_name_rot = 90,
                         simple_anno_size = unit(0.351, "cm"),
                         annotation_name_side = "top",
                         show_legend = FALSE)

col_ha = HeatmapAnnotation(df = data.frame(Seasons = col_data_temp$Annual_season),
                           col = list(Seasons = c('Winter' = '#3399FF', 'Spring/Autumn' = '#FFB266', 'Summer' = '#FF6666')), 
                           border = FALSE,
                           simple_anno_size = unit(0.2, "cm"),
                           show_annotation_name = FALSE,
                           name = "",
                           show_legend = FALSE)

DEG_heatmap <- Heatmap(matDataHeat, 
                       col = colorRamp2(seq(from = -2.5, to = 2.5, by = 0.5) ,rev(brewer.pal(11, "RdYlBu"))), 
                       top_annotation = col_ha, 
                       left_annotation = row_ha2, 
                       right_annotation = row_ha1, 
                       column_split = col_data_temp$Annual_season, 
                       cluster_column_slices = FALSE, 
                       show_column_names = FALSE,
                       row_split = go_genes$Term, 
                       row_title_rot = 90, 
                       row_gap = unit(0.7, "mm"), 
                       cluster_row_slices = FALSE,
                       border = TRUE, 
                       border_gp = gpar(col = "black"),
                       heatmap_legend_param = list(direction = "vertical",
                                                   title_position = "topcenter"),
                       show_heatmap_legend = FALSE, 
                       column_dend_height = unit(0.5, "cm"),
                       row_dend_width = unit(0.5, "cm"),
                       name = "Scaled\nnormalized\nexpression")





## Functional enrichment dotplot as in Extended Data Fig. 7c ---------------------------

mean_expression <- read.table("info_files/mean_gene_expression.txt", header = TRUE)
rownames(mean_expression) <- mean_expression$Gene_ID
base_mean <- mean_expression %>%
  subset(rownames(.) %in% gene_universe)

#load results Leuven cohort
df_seasonal_genes_LS <- read.table("results/DGE_seasons/output/all/all_results.txt", header = TRUE, sep = "\t") %>%
  .[, -10]
df_seasonal_genes_LS$direction <- case_when(df_seasonal_genes_LS$qval < 0.05 & df_seasonal_genes_LS$coef < 0 ~ "high in winter",
                                            df_seasonal_genes_LS$qval < 0.05 & df_seasonal_genes_LS$coef > 0 ~ "high in summer",
                                            .default = "not_significant")

#load results Rhineland cohort (validation cohort)
df_seasonal_genes_RS <- read.table("results/Rhineland_cohort/output/all/all_results.txt") %>%
  .[, -10]
df_seasonal_genes_RS$direction <- case_when(df_seasonal_genes_RS$qval < 0.05 & df_seasonal_genes_RS$coef < 0 ~ "high in winter",
                                            df_seasonal_genes_RS$qval < 0.05 & df_seasonal_genes_RS$coef > 0 ~ "high in summer",
                                            .default = "not_significant")

df_merged <- df_seasonal_genes_LS %>%
  merge(df_seasonal_genes_RS, by = "feature", suffixes = c("_LS", "_RS")) %>%
  merge(annotation, by.x = "feature", by.y ="Gene_ID")

direction <- c("high_in_summer_both", "high_in_winter_both", "shared_RS_LS")
for(d in 1:length(direction)){
  
  if(direction[d] == "high_in_summer_both"){
    degs <- df_merged %>%
      subset(.$direction_RS == "high in summer") %>%
      subset(.$direction_LS == "high in summer") %>%
      .$feature %>%
      unique()
  }else if(direction[d] == "high_in_winter_both"){
    degs <- df_merged %>%
      subset(.$direction_RS == "high in winter") %>%
      subset(.$direction_LS == "high in winter") %>%
      .$feature %>%
      unique()
  }else if(direction[d] == "shared_RS_LS"){
    degs <- df_merged %>%
      subset(.$direction_RS == "high in winter" | .$direction_RS == "high in summer") %>%
      subset(.$direction_LS == "high in winter" | .$direction_LS == "high in summer") %>%
      .$feature %>%
      unique()
  }
  
  overallBaseMean <- as.matrix(base_mean[, "mean_expression", drop = F]) 
  topGOResults <- run_GO_ORA(degs, 
                             overallBaseMean, 
                             selected_database = "org.Hs.eg.db", 
                             back_num = 10,
                             annotation = annotation)
  
  go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
  
  file_filtered <- paste0("results/DGE_seasons/output/validation_seasonal_ORA_GO_", direction[d], ".csv")
  write.csv(go_results_filtered, file_filtered, row.names = FALSE)
}

direction <- c("high_in_winter_both", "high_in_summer_both", "shared_RS_LS")
GO_tables <- list()
d <- 1
for(d in 1:length(direction)){
  
  GO_tables <- append(GO_tables,
                      list(read_csv(paste0("results/DGE_seasons/output/validation_seasonal_ORA_GO_", direction[d], ".csv"))))
  
}

my_terms_1 <- GO_tables[[1]][c(1, 2, 3, 4, 8, 10), ]$GO.ID
my_terms_2 <- GO_tables[[2]][c(1:5), ]$GO.ID
my_terms_3 <- GO_tables[[3]][c(1:5, 10), ]$GO.ID
my_terms <- c(my_terms_1, my_terms_2, my_terms_3)

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = c("high_in_winter_both", "high_in_summer_both", "shared_RS_LS"),
                     score_name = score_name,
                     reverse_bool = TRUE,
                     showTerms = my_terms,
                     multLines = TRUE,
                     numChar = 60)

go_plot_data$Direction <- factor(go_plot_data$Direction, 
                                 levels = c("high_in_winter_both", "high_in_summer_both", "shared_RS_LS"),
                                 labels = c("Shared in winter", "Shared in summer", "All shared genes"))

go_plot_validation <- ggplot(go_plot_data, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  scale_colour_gradient(high = "#990000", low = "#FF9999", name = expression(-~log["10"]~italic(p))) +
  scale_size(name = "Gene\nratio", range = c(0, 3)) +
  theme_bw() +
  guides(color = guide_colorbar(order = 1), 
         size = guide_legend(order = 2)) +
  theme(axis.text.y = element_text(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        legend.position = "right",
        legend.direction = "vertical")




## Scatterplots as in Extended Data Fig. 7e, left and middle left ---------------------------

all_cohorts_combined_DEG_result <- read.csv("results/DGE_seasons/output/All_cohorts_combined_seasonal_results.txt", sep = '\t')

p1 <- ggplot(all_cohorts_combined_DEG_result, aes(x=coef.All, y=coef.Healthy)) + geom_point(alpha = 0.5, col="skyblue2")
p1 <- p1 + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p1 <- p1 + xlab("Expression change (All samples)") + ylab("Expression change (Without medical conditions)")
p1 <- p1 + stat_cor(method = "spearman", cor.coef.name = "rho")
p1 <- p1 + theme_pubr() + theme(axis.text = element_text(color = "black"), 
                                axis.title = element_text(), 
                                legend.text = element_text(), 
                                legend.title = element_text(), 
                                panel.grid = element_blank(), 
                                axis.line = element_line(),
                                axis.title.x = element_blank())

p2 <- ggplot(all_cohorts_combined_DEG_result, aes(x=coef.All, y=coef.Morning)) + geom_point(alpha = 0.5, col="skyblue2")
p2 <- p2 + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p2 <- p2 + xlab("Expression change (All samples)") + ylab("Expression change (Morning samples)")
p2 <- p2 + stat_cor(method = "spearman", cor.coef.name = "rho")
p2 <- p2 + theme_pubr() + theme(axis.text = element_text(color = "black"), 
                                axis.title = element_text(), 
                                legend.text = element_text(), 
                                legend.title = element_text(), 
                                panel.grid = element_blank(), 
                                axis.line = element_line(),
                                axis.title.x = element_blank())




