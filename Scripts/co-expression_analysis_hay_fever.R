## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(rstatix)
library(ggpubr)
library(DESeq2)
library(RColorBrewer)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")


##Loading data ------------------------------------------------------------

#import metadata
col_clinical_data <- read.table("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = "\t", header = TRUE)
col_clinical_data$Individual_ID <- factor(col_clinical_data$Individual_ID)
rownames(col_clinical_data) <- col_clinical_data$Sample_ID
col_clinical_data$Run <- factor(col_clinical_data$Run)

gene_attributes <- read.table("info_files/filtered_genes_samples25_fpm5.txt", sep = "\t", header = TRUE)
annotation <- gene_attributes[, c("Gene_ID", "gene_name")] %>% unique()
protein_coding_genes <- subset(gene_attributes, gene_attributes$gene_type == "protein_coding")$Gene_ID
sex_chromosome_genes <- subset(gene_attributes, gene_attributes$Chr %in% c("chrX", "chrY"))$Gene_ID

#add factor for annual season
col_clinical_data$Annual_season <- ifelse(col_clinical_data$Month %in% c("January", "February", "March"), "Winter", 
                                          ifelse(col_clinical_data$Month %in% c("July", "August", "September"), "Summer",
                                                 ifelse(col_clinical_data$Month %in% c("April", "May"), "Spring", "Autumn")))
col_clinical_data$Annual_season <- factor(col_clinical_data$Annual_season, levels = c("Winter", "Spring", "Summer", "Autumn"))

#Get RUVSeq factors of unwanted variation
unwanted_factors <- read.table("results/RUVSeq/output/unwanted_factors_RUVg_2_3.txt", header = TRUE)
col_clinical_data <- col_clinical_data %>%
  merge(unwanted_factors, by = "Sample_ID")

col_clinical_data$W_2 <- cut(col_clinical_data$W_2, 5)
col_clinical_data$W_3 <- cut(col_clinical_data$W_3, 5)


#import count data
count_data <- read.table("count_files/merged_gene_counts_all_samples.txt", sep = '\t')

colnames(count_data) <- substr(colnames(count_data), 1, 6)
col_clinical_data <- subset(col_clinical_data, col_clinical_data$Sample_ID %in% colnames(count_data))
count_data <- count_data[, as.character(col_clinical_data$Sample_ID)]
colnames(count_data) <- col_clinical_data$Sample_ID

count_data <- count_data[as.character(gene_attributes$Gene_ID %>% unique()), ]
rownames(count_data) <- gene_attributes$Gene_ID %>% unique()


#import gene co-expression results
res_module <- read.table("results/Coexpression/output/gene_module.txt", header = TRUE) %>%
  merge(annotation, by.x = "Gene", by.y = "Gene_ID")
df_eigengene <- read.table("results/Coexpression/output/module_eigengenes.txt", header = TRUE)
df_module_list <- read.table("results/Coexpression/output/module_list.txt", header = TRUE)
df_membership <- read.table("results/Coexpression/output/module_membership_score.txt", header = TRUE)
df_icc_modules <- read.table("results/Coexpression/output/ICC_module_eigengenes.txt", sep = "\t", header = TRUE)





## Select samples with hay fever and matched controls --------------------------

#Get samples with hay fever
df_hay_fever <- col_clinical_data %>%
  subset(.$hay_fever == "YES") %>%
  subset(.$asthma == "NO") %>%
  .[order(.$Individual_ID), ] %>%
  add_column("group" = "hay_fever") 

allergy_individuals <- df_hay_fever %>%
  .$Individual_ID %>%
  unique()

length(allergy_individuals) #51

df_hay_fever[, c("Individual_ID", "Season")] %>%
  unique() %>%
  {table(.$Season)}
# A  B 
#45  6

df_hay_fever[, c("Individual_ID", "Annual_season")] %>%
  unique() %>%
  {table(.$Annual_season)}
#Winter Spring Summer Autumn 
#    49     45     49      6 

#select controls matched for age, sex and annual season
set.seed(123)

control_samples <- c()
control_individuals <- c()
i <- 1
for(i in 1:length(allergy_individuals)){
  
  samples_temp <- df_hay_fever %>%
    subset(.$Individual_ID == allergy_individuals[i])
  
  my_individuals <- col_clinical_data %>%
    subset(! .$Individual_ID %in% allergy_individuals) %>%
    subset(.$hay_fever == "NO") %>%
    subset(.$asthma == "NO") %>%
    subset(.$allergies_season == "NO") %>%
    subset(.$Gender == samples_temp[1, "Gender"]) %>%
    subset(.$Age..years. > as.numeric(samples_temp[1, "Age..years."]) - 3) %>%
    subset(.$Age..years. < as.numeric(samples_temp[1, "Age..years."]) + 3) %>%
    subset(.$Season == samples_temp[1, "Season"]) %>%
    subset(! .$Individual_ID %in% control_individuals) %>%
    .[, c("Individual_ID", "Annual_season")] %>%
    unique() %>%
    group_by(Individual_ID) %>%
    summarize("n" = n()) %>%
    subset(.$n == 3) %>%
    .$Individual_ID %>%
    unique()
  
  if(length(my_individuals) == 0){
    my_individuals <- col_clinical_data %>%
      subset(! .$Individual_ID %in% allergy_individuals) %>%
      subset(.$hay_fever == "NO") %>%
      subset(.$asthma == "NO") %>%
      subset(.$allergies_season == "NO") %>%
      subset(.$Gender == samples_temp[1, "Gender"]) %>%
      subset(.$Age..years. > as.numeric(samples_temp[1, "Age..years."]) - 10) %>%
      subset(.$Age..years. < as.numeric(samples_temp[1, "Age..years."]) + 10) %>%
      subset(.$Season == samples_temp[1, "Season"]) %>%
      subset(! .$Individual_ID %in% control_individuals) %>%
      .[, c("Individual_ID", "Annual_season")] %>%
      unique() %>%
      group_by(Individual_ID) %>%
      summarize("n" = n()) %>%
      subset(.$n ==3) %>%
      .$Individual_ID %>%
      unique()
  }
  
  selected_individual <- my_individuals[sample(1:length(my_individuals), 1, replace = FALSE)]
  control_individuals <- c(control_individuals, as.character(selected_individual))
  
  selected_samples <- col_clinical_data %>%
    subset(.$Individual_ID == selected_individual) %>%
    subset(.$Annual_season %in% samples_temp$Annual_season) %>%
    .[match( samples_temp$Annual_season, .$Annual_season), ] %>%
    .$Sample_ID
  
  control_samples <- c(control_samples, selected_samples)
}

samples_allergy_control <- col_clinical_data %>% 
  subset(.$Sample_ID %in% control_samples) %>%
  .[match( control_samples, .$Sample_ID), ] %>%
  add_column("group" = "control")

all(df_hay_fever$Annual_season == samples_allergy_control$Annual_season) #TRUE
samples_allergy_control$Sample_ID %>% unique() %>% length() #149

col_data_allergy <- rbind(df_hay_fever, samples_allergy_control) 
write.table(col_data_allergy, "results/allergies/output/samples_hay_fever.txt", row.names = FALSE, quote = FALSE, sep = "\t")




## Code to plot hay fever analysis results as in Figure 5 and Extended Data Fig. 8 ---------

my_colors <- c("control" = "grey", "hay_fever" = "#7031B9", "NO" = "grey", "YES" = "#7031B9", "Winter" = '#3399FF', "Spring" = '#FFB266', "Summer" = '#FF6666')

### Code for module eigengene boxplots as in Figure 5c and 5d ------------------------------

#import individuals
df_hay_fever <- read.table("results/allergies/output/samples_hay_fever.txt", header = TRUE, sep = "\t")
df_hay_fever$group <- factor(df_hay_fever$group, levels = c("control", "hay_fever"))

my_theme_hay_fever <- theme(strip.placement = "outside",
                            strip.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.title = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            legend.title = element_blank(),
                            legend.position = "bottom")

#compare across seasons
df_plot_temp <- df_eigengene %>%
  add_column("Description" = rownames(.)) %>%
  merge(df_hay_fever, by = "Description") %>%
  dplyr::rename("value" = "MEplum4") %>%
  subset(.$Season == "A")

df_plot_temp$Annual_season <- factor(df_plot_temp$Annual_season, levels = c("Winter", "Spring", "Summer")) 

df_test <- df_plot_temp %>%
  subset(.$Season == "A") %>%
  group_by(Annual_season) %>%
  wilcox_test(value ~ group, data = . , paired = FALSE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = df_plot_temp$value %>% {max(., na.rm = TRUE) - (max(., na.rm = TRUE) - min(., na.rm = TRUE))*0.1})

hay_fever_eigengene_seasons <- ggplot(df_plot_temp) +
  geom_boxplot(width = 0.8, alpha = 0.5, aes(x = group, y = value, color = group, fill = group)) +
  scale_color_manual(values = my_colors, labels = c("Control", "Hay fever")) +
  scale_fill_manual(values = my_colors, labels = c("Control", "Hay fever")) +
  stat_pvalue_manual(data = df_test, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
  labs(x = "Hay fever\n(self-reported)", y = "Eosinophil module eigengene") +
  facet_wrap(~Annual_season, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  my_theme_hay_fever


#plot individual trajectories
var_names <- c("control" = "Control", 
               "hay_fever" = "Hay fever")

df_plot <- df_eigengene %>%
  add_column("Description" = rownames(.)) %>%
  merge(df_hay_fever, by = "Description") %>%
  dplyr::rename("value" = "MEplum4") %>%
  subset(.$Season == "A") %>%
  add_column("group_plot" = ifelse(.$group == "control", "control", .$Annual_season))

df_plot$my_order <- factor(df_plot$Annual_season, levels = c("Winter", "Spring", "Summer")) %>%
  as.numeric()

df_gene <- df_plot 
df_gene$Annual_season <- factor(df_gene$Annual_season, levels = c("Winter", "Spring", "Summer"))

df_test_wide <- df_gene[, c("Individual_ID", "value", "my_order")] %>%
  pivot_wider(id_cols = "Individual_ID", names_from = "my_order") %>%
  merge(df_plot[, c("Individual_ID", "group")] %>% unique(), by = "Individual_ID")

df_test <- df_test_wide %>%
  melt(id.vars = c("Individual_ID", "group"), value.name = "value", variable.name = "my_order") %>%
  group_by(group) %>%
  wilcox_test(value ~ my_order, data = . , paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = df_test_wide[, 2:4] %>%
               {c(max(., na.rm = TRUE) - (max(., na.rm = TRUE) - min(., na.rm = TRUE))*0.05,
                  max(., na.rm = TRUE) - (max(., na.rm = TRUE) - min(., na.rm = TRUE))*0.1,
                  max(., na.rm = TRUE) - (max(., na.rm = TRUE) - min(., na.rm = TRUE))*0.15)} %>% rep(2)) %>%
  mutate(group1 = group1 %>% as.numeric(),
         group2 = group2 %>% as.numeric())

hay_fever_eigengene_trajectories <- ggplot(df_gene, aes(x = my_order)) +
  geom_line(data = df_gene, aes(y = value, group = Individual_ID), color = "lightgray") +
  facet_wrap(~ group, strip.position = "bottom", labeller = as_labeller(var_names)) +
  stat_pvalue_manual(data = df_test, y.position = "max", label = "p.adj_signif", tip.length = 0.01, bracket.shorten = 0.5) +
  geom_boxplot(aes(x = my_order, y = value, group = Annual_season, color = Annual_season, fill = Annual_season), 
               width = 0.5, 
               alpha = 0.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(labels = c("Winter", "Spring", "Summer"), breaks = 1:3) +
  labs(y = "Eosinophil module eigengene") +
  theme_bw() +
  my_theme_hay_fever +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


table(df_plot$group, df_plot$Annual_season)
#          Winter Spring Summer
#control       44     45     43
#hay_fever     44     45     43





### Code for heatmap of module genes as in Extended Data Fig. 8d ------------------------------

normalized_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ 1) %>%
  estimateSizeFactors() %>%
  counts(normalized = TRUE) %>%
  as.data.frame()

my_module <- res_module %>%
  subset(.$ModuleColor %in% c("plum4")) %>%
  .[order(.$ModuleLabel), ] 

my_genes <- my_module %>%
  .$Gene

df_hay_fever <- read.table("results/allergies/output/samples_hay_fever.txt", header = TRUE, sep = "\t") %>%
  subset(.$Season == "A")

df_plot <- normalized_counts[my_genes, ] %>%
  t() %>%
  as.data.frame() %>%
  add_column("Sample_ID" = rownames(.)) %>%
  melt(varnames = c("Sample_ID")) %>%
  merge(df_hay_fever, by = "Sample_ID") 

df_mean <- df_plot %>%
  group_by(variable, Annual_season, group) %>%
  summarize("mean" = mean(value)) %>%
  add_column("my_order" = paste0(.$group, "_", .$Annual_season)) %>%
  .[, c(-2, -3)]

matDataHeat <- df_mean %>%
  pivot_wider(id_cols = "variable", names_from = "my_order", values_from = "mean") %>%
  as.data.frame() %>%
  .[, c("variable", paste0(rep(c("control_", "hay_fever_"), each = 3), rep(c("Winter", "Spring", "Summer"), 2)) %>% rev())]
rownames(matDataHeat) <- matDataHeat$variable
matDataHeat <- matDataHeat[, -1]

matDataHeat <- scale(t(matDataHeat))

col_anno <- data.frame("Status" = c(rep(c("Control", "Hay fever"), each = 3) %>% rev()),
                       "Season" = rep(c("Winter", "Spring", "Summer"), 2) %>% rev(), 
                       check.names = FALSE)

color_hay_fever <- c("Control" = "grey", "Hay fever" = "#7031B9")
color_season <- c("Winter" = '#3399FF', "Spring" = '#FFB266', "Summer" = '#FF6666')

col_ha <- rowAnnotation(df = col_anno,
                        col = list("Status" = color_hay_fever,
                                   "Season" = color_season),
                        annotation_legend_param = list("Status" = 
                                                         list(title = "Status"),
                                                       "Season" = 
                                                         list(labels = c("Winter", "Spring", "Summer"),
                                                              at = c("Winter", "Spring", "Summer"),
                                                              legend_gp = gpar(fill = c("#3399FF","#FFB266", "#FF6666")))),
                        border = TRUE,
                        simple_anno_size = unit(0.25, "cm"),
                        annotation_name_gp = gpar(size = 7),
                        show_annotation_name = FALSE,
                        name = "")

col_split <- c(rep("A", 3), rep("B", 3)) %>%
  factor()

my_names <- matDataHeat %>% 
  colnames() %>%
  data.frame("Gene_ID"= . ) %>%
  merge(annotation, by = "Gene_ID", sort = FALSE) %>%
  .[match(matDataHeat %>% colnames(), .$Gene_ID), ]

row_ha <- rowAnnotation(gene_names = anno_text(my_names$gene_name, 
                                               gp = gpar(fontsize = 7)))

row_ha <- HeatmapAnnotation(gene_names = anno_text(my_names$gene_name, 
                                                   gp = gpar(fontsize = 7),
                                                   location = unit(1, 'npc'),
                                                   just = "right"))

p <- Heatmap(matDataHeat, 
             name = "Scaled\nnormalized\nexpression",
             col = colorRamp2(seq(from = -1.7, to = 1.7, by = 0.34), rev(brewer.pal(11, "RdYlBu"))),
             cluster_columns = TRUE,
             cluster_rows = FALSE,
             bottom_annotation = row_ha, 
             left_annotation = col_ha,
             row_title = NULL,
             row_names_gp = gpar(fontsize = 7),
             column_dend_side = "top",
             row_gap = unit(1, "mm"), 
             border = TRUE, 
             show_row_names = FALSE,
             show_column_names = FALSE,
             row_split = col_split,
             row_title_gp = gpar(fontsize = 10),
             row_title_rot = 0,
             cluster_row_slices = FALSE,
             column_title = NULL,
             show_column_dend = TRUE,
             cluster_column_slices = FALSE)
p_heatmap <- draw(p, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")


