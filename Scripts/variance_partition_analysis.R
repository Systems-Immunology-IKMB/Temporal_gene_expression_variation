## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(DESeq2)
library(variancePartition)
library(reshape2)
library(patchwork)

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

#add factor for annual season
col_clinical_data$Annual_season <- ifelse(col_clinical_data$Month %in% c("January", "February", "March"), "Winter", 
                                          ifelse(col_clinical_data$Month %in% c("July", "August", "September"), "Summer",
                                                 ifelse(col_clinical_data$Month %in% c("April", "May"), "Spring", "Autumn")))
col_clinical_data$Annual_season <- factor(col_clinical_data$Annual_season, levels = c("Winter", "Spring", "Summer", "Autumn"))

#import count data
count_data <- read.csv("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)
count_data <- count_data[as.character(gene_attributes$Gene_ID), as.character(col_clinical_data$Sample_ID)]



## Compute variance partition ------------------------

dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ 1)
dds_counts <- estimateSizeFactors(dds_counts)
quantvst <- vst(dds_counts)
quantvst <- assay(quantvst)

form <- ~ (1|Individual_ID) + (1|Run) + (1|Annual_season) + (1|Season) + Time_of_sampling + (1|Infection_type) + (1|Vaccine_last_month) + (1|allergies_season) + (1|allergies_food) + Age..years. + (1|Gender) + BMI..kg.m2. + Hemoglobin..g.dL. + Neutrophils...µL. + Eosinophils...µL. + Lymphocytes...µL. + Monocytes...µL. + Thrombocytes..x1000.µL. + hsCRP..mg.L. + Hemoglobine.A1c...... + Creatinine..mg.dL. + Uric.acid..mg.dL. + Albumin..g.L. + Triglycerides..mg.dL.

scaled_clinical_data <- col_clinical_data
scaled_clinical_data[, c("Age..years.", "BMI..kg.m2.", "Hemoglobin..g.dL.", "Neutrophils...µL.", "Eosinophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L.", "Hemoglobine.A1c......", "Creatinine..mg.dL.", "Uric.acid..mg.dL.", "Albumin..g.L.", "Gamma.globulins..g.L.", "Triglycerides..mg.dL.", "Cholesterol.LDL..measured...mg.dL.",    "Bilirubin.total..mg.dL.", "Gamma.GT..U.L.", "CK..U.L.", "Time_of_sampling")] <- scale(scaled_clinical_data[, c("Age..years.", "BMI..kg.m2.", "Hemoglobin..g.dL.", "Neutrophils...µL.", "Eosinophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L.", "Hemoglobine.A1c......", "Creatinine..mg.dL.", "Uric.acid..mg.dL.", "Albumin..g.L.", "Gamma.globulins..g.L.", "Triglycerides..mg.dL.", "Cholesterol.LDL..measured...mg.dL.", "Bilirubin.total..mg.dL.", "Gamma.GT..U.L.", "CK..U.L.", "Time_of_sampling")])

varPart <- fitExtractVarPartModel(quantvst, form, scaled_clinical_data)
varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
write.table(varPart_sorted, "results/Variance_partition/output/variance_partition.txt", sep = '\t', quote = FALSE)

varMean <- colMeans(varPart)
varMean <- as.data.frame(varMean)
varMean$Parameter <- rownames(varMean)
varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]
write.table(varMean, "results/Variance_partition/output/variance_partition_mean.txt", sep = '\t', quote = FALSE, row.names = FALSE)




## Violinplot of variance partition as in Fig. 1e ---------

#load labeller and set colors
column_labels <- read.csv("results/Column_labels.csv", row.names = 1)
column_labels_vector <- c(column_labels$Label, "Residuals")
names(column_labels_vector) <- c(rownames(column_labels), "Residuals")

c25 <- c("Individual_ID"="skyblue2", 
         "Run"="darkturquoise", 
         "Season" = "green4",
         "Gender" = "#6A3D9A", 
         "Annual_season" = "#FF7F00", 
         "Residuals" = "gray70", 
         "Age..years." = "gold1", 
         "BMI..kg.m2." = "dodgerblue2", 
         "Hemoglobin..g.dL." = "#FB9A99", 
         "Neutrophils...µL." = "palegreen2", 
         "Eosinophils...µL." = "#CAB2D6", 
         "Lymphocytes...µL." = "#FDBF6F", 
         "hsCRP..mg.L." = "darkorange4",  
         "Monocytes...µL." = "khaki2", 
         "Thrombocytes..x1000.µL." = "maroon", 
         "Creatinine..mg.dL." = "orchid1", 
         "Hemoglobine.A1c......" = "deeppink1", 
         "Uric.acid..mg.dL." = "blue1", 
         "Albumin..g.L." = "steelblue4", 
         "Gamma.globulins..g.L." = "#E31A1C", 
         "Triglycerides..mg.dL." = "green1", 
         "Cholesterol.LDL..measured...mg.dL." = "yellow4", 
         "Bilirubin.total..mg.dL." = "yellow3", 
         "allergies_season" = "slategray", 
         "allergies_food" = "brown",
         "Time_of_sampling" = "darkgreen", 
         "Infection_type" = "skyblue", 
         "Vaccine_last_month" = "darkorange3"
)

varPart <- read.table( "results/Variance_partition/output/variance_partition.txt", sep = "\t") 
varPart <- varPart %>%
  add_column("ID" = "ID") %>%
  melt(id.vars = "ID")
varPart$variable <- factor(varPart$variable)
varPart$value <- varPart$value * 100

varpart_plot <- ggplot(varPart, aes(x = variable,  y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 1, scale = "width") +
  geom_boxplot(width = 0.07, fill = "grey") +
  labs(y = "Variance explained (%)") +
  ylim(c(0, 100)) +
  scale_x_discrete(labels = as_labeller(column_labels_vector)) +
  scale_fill_manual(values = c25) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid = element_blank())





## Barplot of mean variance explained as in Extended Data Fig. 2b ---------

varPart <- read.table("results/Variance_partition/output/variance_partition.txt", sep = "\t") 
varMean <- colMeans(varPart)
varMean2 <- varMean[order(-varMean)]
varMean2 <- as.data.frame(varMean)
varMean2$Parameter <- rownames(varMean2)

varMean_plot <- ggplot(varMean2[varMean2$Parameter != "Residuals", ], aes(x=reorder(Parameter, -varMean), y=varMean, fill = Parameter)) + 
  geom_bar(stat="identity", color="black") + 
  scale_x_discrete(labels=column_labels[rownames(varMean2)[!rownames(varMean2) == "Residuals"], "Label"]) +
  xlab("Parameter") + ylab("Mean variation explained") + scale_fill_manual(values = c25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())







## Barplot of variance partition for top 10 cell type-specific genes as in Fig. 1f -------------------------

#load data
varPart <- read.table("results/Variance_partition/output/variance_partition.txt", sep = "\t") 
varPart <- varPart %>%
  add_column("Gene_ID" = rownames(.)) %>%
  merge(annotation, by = "Gene_ID")

df_cell_type_specificity <- read.table("results/cell_type_specificity/output/tau_lowRes_celltypes.txt", sep = "\t", header = TRUE) %>%
  dplyr::rename("gene_name" = "Gene_ID") %>%
  merge(annotation, by = "gene_name")

#convert to binned values to have a discrete scale
my_breaks <- seq(from = 0, to = 1, by = 0.25)
df_cell_type_specificity$cell_type_specificity_binned <- df_cell_type_specificity$cell_type_specificity %>%
  cut(breaks = my_breaks, include.lowest = TRUE, right = TRUE)

#get top 10 genes with highest variance explained by the 4 different cell types
my_variables <- c("Neutrophils...µL.", "Monocytes...µL.", "Lymphocytes...µL.", "Eosinophils...µL.")
df_celltypes_top10 <- data.frame()
for(v in 1:length(my_variables)){
  
  df_variable <- varPart %>%
    .[order(.[[my_variables[v]]], decreasing = TRUE), ] %>%
    .[1:10, ] %>%
    add_column("group" = my_variables[v])
  
  df_celltypes_top10 <- rbind(df_celltypes_top10, df_variable)
}

rownames(df_celltypes_top10) <- df_celltypes_top10$gene_name
df_celltypes_top10 <- df_celltypes_top10[, c(-1, -27)]
df_celltypes_top10 <- df_celltypes_top10[, c(1:3, 5, 8, 13, 4, 6:7, 9:12, 14:26)]

#order colors
my_colors <- c25[match(colnames(df_celltypes_top10)[1:25], names(c25))]

#set parameters
my_ht_colors <- structure(c("#FEE0D2", "#FCBBA1", "#FB6A4A", "#A50F15"), names = c("[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]"))
my_width_barplot <- 3.5 #in cm
my_with_tau_score <- 0.3 #in cm
my_width_text <- 3/5 

#plot legend for cell type specificity
lgd_tau_score <- Legend(at = c("[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]"), 
                        title = "Cell type specificity",
                        legend_gp = gpar(fill = c("#FEE0D2", "#FCBBA1", "#FB6A4A", "#A50F15")),
                        border = "black",
                        labels_gp = gpar(col = "black"),
                        title_gp = gpar(col = "black"))

#plot barplot with cell type specificity annotation for top 10 neutrophil genes (plots for monocytes, lymphocytes and eosinophils are analogous)
varPart_neutrophils <- df_celltypes_top10 %>%
  subset(.$group == "Neutrophils...µL.") %>%
  .[, -26]

tau_neutrophils <- df_cell_type_specificity %>%
  .[match(rownames(varPart_neutrophils), .$gene_name), ] %>%
  {data.frame("tau_neutrophils" = .$cell_type_specificity_binned %>% as.character(), row.names = .$gene_name)}

text_neutrophils <- rowAnnotation(gene_names = anno_text(rownames(varPart_neutrophils), 
                                                         location = unit(1, 'npc'),
                                                         width = unit(my_width_text, "cm"),
                                                         just = "right"))

bar_neutrophils <- rowAnnotation("Neutrophils" = anno_barplot(varPart_neutrophils, 
                                                              gp = gpar(fill = my_colors),
                                                              border = TRUE,
                                                              ylim = c(0, 1),
                                                              width = unit(my_width_barplot, "cm"),
                                                              bar_width = 1,
                                                              axis_param = list(
                                                                side = "bottom",
                                                                at = c(0, 0.25, 0.5, 0.75, 1), 
                                                                labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
                                                                labels_rot = 0)),
                                 annotation_name_rot = 0,
                                 annotation_name_side = "top",
                                 gp = gpar(col = "black"))

ht_neutrophils <- Heatmap(tau_neutrophils, 
                          name = "Cell type specificity", 
                          col = my_ht_colors,
                          show_row_names = FALSE,
                          show_column_names = FALSE,
                          show_row_dend = FALSE, 
                          cluster_rows = FALSE, 
                          show_heatmap_legend = FALSE,
                          border = TRUE,
                          border_gp = gpar(col = "black"),
                          heatmap_width = unit(my_with_tau_score, "cm"))

draw(text_neutrophils + ht_neutrophils + bar_neutrophils)




## Barplot of variance partition for top 20 individual-specific genes as in Fig. 1g -------------------------

#load data
varPart <- read.table("results/Variance_partition/output/variance_partition.txt", sep = "\t") 

column_labels_vector <- c(column_labels$Label, "Residuals")
names(column_labels_vector) <- c(rownames(column_labels), "Residuals")

#get order for legend
varMean <- colMeans(varPart)
order_legend <- varMean[order(varMean, decreasing = TRUE)] %>% names()

#get top 20 individual-specific genes
df_individual_top20 <- varPart %>%
  add_column("Gene_ID" = rownames(.)) %>%
  merge(annotation, by = "Gene_ID") %>%
  .[order(.$Individual_ID, decreasing = TRUE), ] %>%
  .[1:20, ] %>%
  melt(id.vars = c("Gene_ID", "gene_name"))

df_individual_top20$gene_name <- factor(df_individual_top20$gene_name[1:20], levels = df_individual_top20$gene_name[1:20] %>% rev())
df_individual_top20$variable <- factor(df_individual_top20$variable, levels = levels(df_individual_top20$variable) %>% .[c(2, 1, 3:25)] %>% rev())

individual_top20_plot <- ggplot(df_individual_top20, aes(x = value, y = gene_name)) +
  geom_bar(position = "stack", stat = "identity", aes(fill = variable), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c25, labels = as_labeller(column_labels_vector), breaks = order_legend) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.25, 0.1, 0, 0), "cm"),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank())






