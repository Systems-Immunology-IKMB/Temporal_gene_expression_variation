## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(dplyr)
library(DESeq2)
library(mixOmics)
library(variancePartition)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")


## Loading data ------------------------------------------------------------

gene_list <- read.csv("info_files/filtered_genes_samples25_fpm5.txt", sep = '\t')

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
count_data <- count_data[as.character(gene_list$Gene_ID), as.character(col_clinical_data$Sample_ID)]



## Run PCA - multilevel approach ----------------------

vst_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ 1) %>%
  estimateSizeFactors() %>%
  vst() %>%
  assay() %>%
  t() %>%
  as.data.frame()


#within individuals
design <- data.frame(sample = col_clinical_data$Individual_ID)
vst_counts_within <- withinVariation(vst_counts, design)

pca <- vst_counts_within %>%
  prcomp()

all(col_clinical_data$Sample_ID == rownames(pca$x))
pc_df_within <- cbind(col_clinical_data, pca$x)

pc_variance_within <- summary(pca) %>%
  .$importance %>%
  .["Proportion of Variance", ] %>% 
  `*`(100) %>%
  round(2) %>%
  sprintf(fmt = '%.2f')


#between individuals
vst_counts_between <- vst_counts - vst_counts_within

#subset to only one sample per individual (values are identical)
set.seed(123)
col_clinical_data_between <- col_clinical_data %>%
  group_by(Individual_ID) %>%
  sample_n(1) %>%
  ungroup()

vst_counts_between <- vst_counts_between[col_clinical_data_between$Sample_ID, ]

pca <- vst_counts_between %>%
  prcomp()

all(col_clinical_data_between$Sample_ID == rownames(pca$x))
pc_df_between <- cbind(col_clinical_data_between, pca$x)

pc_variance_between <- summary(pca) %>%
  .$importance %>%
  .["Proportion of Variance", ] %>% 
  `*`(100) %>%
  round(2) %>%
  sprintf(fmt = '%.2f')




## Plot correlation between variables and between-individual PCs as in Fig. 1a ---------

#compute canonical correlation between variables and the first 6 BPCs
form <- ~ Individual_ID + Neutrophils...µL. + Lymphocytes...µL. + Monocytes...µL. + Eosinophils...µL. + Annual_season + Time_of_sampling + Vaccine_last_month + allergies_season + Age..years. + Gender + BMI..kg.m2. + Run + Season + PC1 + PC2 + PC3 + PC4 + PC5 + PC6

Corr_mat_between <- canCorPairs(form, pc_df_between) %>%
  .[paste0(rep("PC", 6), 1:6), c("Individual_ID", "Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Annual_season", "Time_of_sampling", "Vaccine_last_month", "allergies_season", "Age..years.", "Gender", "BMI..kg.m2.", "Run", "Season")]
colnames(Corr_mat_between) <- c("Individual", "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Annual season", "Time of sampling", "Vaccination", "Seasonal allergies", "Age", "Sex", "BMI", "Sequencing Run", "Cohort")

Corr_mat_between <- Corr_mat_between %>%
  as.data.frame() %>%
  add_column("Components" = rownames(.)) %>%
  melt(id.vars = "Components") %>%
  add_column("value_rounded" = round(.$value, digits = 2)) %>%
  subset(!.$variable %in% c("Vaccination", "Seasonal allergies"))

Corr_mat_between$Components <- factor(Corr_mat_between$Components, levels = paste0(rep("PC", 6), 1:6) %>% rev())
Corr_mat_between$variable <- factor(Corr_mat_between$variable, levels = c("Individual", "Sex", "Age", "BMI", "Sequencing Run", "Cohort", "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Annual season", "Time of sampling"))

#add grid lines
Corr_mat_between <- Corr_mat_between %>%
  mutate(line_positions_y = as.numeric(factor(Components, levels = unique(Components))), 
         line_positions_y= line_positions_y + .5,  
         line_positions_y = ifelse(line_positions_y == max(line_positions_y), NA, line_positions_y)) %>%
  mutate(line_positions_x = as.numeric(factor(variable, levels = unique(variable))), 
         line_positions_x= line_positions_x + .5,  
         line_positions_x = ifelse(line_positions_x == max(line_positions_x), NA, line_positions_x)) 

Corr_plot_between <- ggplot(Corr_mat_between, aes(x = variable, y = Components)) +
  geom_point(aes(color = value, fill = value, size = value_rounded), shape = 21) +
  geom_text(aes(label = value_rounded), color = "black") +
  geom_vline(aes(xintercept = line_positions_x)) + 
  geom_hline(aes(yintercept = line_positions_y)) + 
  scale_y_discrete(position = "left", 
                   labels = paste0("B", names(pc_variance_between[1:6]), "\n(", pc_variance_between[1:6], rep("%)", 6)), 
                   breaks = paste0(rep("PC", 6), 1:6), 
                   expand = c(0, 0.5)) + 
  scale_x_discrete(position = "bottom", expand = c(0, 0.5)) +
  scale_size(range = c(1, 6), guide = "none", limits = c(0, 1)) +
  scale_fill_gradient(low = "#F7FBFF", high = "#053061", limits = c(0, 1)) +
  scale_color_gradient(low = "#F7FBFF", high = "#053061", limits = c(0, 1)) +
  labs(fill = "Canonical\ncorrelation", 
       color = "Canonical\ncorrelation", 
       size = "Canonical\ncorrelation", 
       y = "Between-individual PCA", 
       x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        panel.grid.major = element_blank()) +
  coord_fixed()





## Plot correlation between variables and within-individual PCs as in Fig. 1b ---------

#compute canonical correlation between variables and the first W6 PCs
form <- ~ Individual_ID + Neutrophils...µL. + Lymphocytes...µL. + Monocytes...µL. + Eosinophils...µL. + Annual_season + Time_of_sampling + Vaccine_last_month + allergies_season + Age..years. + Gender + BMI..kg.m2. + Run + Season + PC1 + PC2 + PC3 + PC4 + PC5 + PC6

Corr_mat_within <- canCorPairs(form, pc_df_within) %>%
  .[paste0(rep("PC", 6), 1:6), c("Individual_ID", "Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Annual_season", "Time_of_sampling", "Vaccine_last_month", "allergies_season", "Age..years.", "Gender", "BMI..kg.m2.", "Run", "Season")]
colnames(Corr_mat_within) <- c("Individual", "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Annual season", "Time of sampling", "Vaccination", "Seasonal allergies", "Age", "Sex", "BMI", "Sequencing Run", "Cohort")

Corr_mat_within <- Corr_mat_within %>%
  as.data.frame() %>%
  add_column("Components" = rownames(.)) %>%
  melt(id.vars = "Components") %>%
  add_column("value_rounded" = round(.$value, digits = 2)) %>%
  subset(!.$variable %in% c("Vaccination", "Seasonal allergies"))

Corr_mat_within$Components <- factor(Corr_mat_within$Components, levels = paste0(rep("PC", 6), 1:6) %>% rev())
Corr_mat_within$variable <- factor(Corr_mat_within$variable, levels = c("Individual", "Sex", "Age", "BMI", "Sequencing Run", "Cohort", "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Annual season", "Time of sampling")) 

#add grid lines
Corr_mat_within <- Corr_mat_within %>%
  mutate(line_positions_y = as.numeric(factor(Components, levels = unique(Components))), 
         line_positions_y= line_positions_y + .5,  
         line_positions_y = ifelse(line_positions_y == max(line_positions_y), NA, line_positions_y)) %>%
  mutate(line_positions_x = as.numeric(factor(variable, levels = unique(variable))), 
         line_positions_x= line_positions_x + .5,  
         line_positions_x = ifelse(line_positions_x == max(line_positions_x), NA, line_positions_x)) 

Corr_plot_within <- ggplot(Corr_mat_within, aes(x = variable, y = Components)) +
  geom_point(aes(color = value, fill = value, size = value_rounded), shape = 21) +
  geom_text(aes(label = value_rounded), color = "black") +
  geom_vline(aes(xintercept = line_positions_x)) + 
  geom_hline(aes(yintercept = line_positions_y)) + 
  scale_y_discrete(position = "left", 
                   labels = paste0("W", names(pc_variance_within[1:6]), "\n(", pc_variance_within[1:6], rep("%)", 6)), 
                   breaks = paste0(rep("PC", 6), 1:6), 
                   expand = c(0, 0.5)) + 
  scale_x_discrete(position = "bottom", 
                   expand = c(0, 0.5)) +
  scale_size(range = c(1, 6), guide = "none", limits = c(0, 1)) +
  scale_fill_gradient(low = "#F7FBFF", high = "#053061", limits = c(0, 1)) +
  scale_color_gradient(low = "#F7FBFF", high = "#053061", limits = c(0, 1)) +
  labs(fill = "Canonical\ncorrelation", 
       color = "Canonical\ncorrelation", 
       size = "Canonical\ncorrelation", 
       y = "Within-individual PCA", 
       x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        panel.grid.major = element_blank()) +
  coord_fixed()




## Plot BPCs colored by sex, sequencing run, age, and diurnal time of sampling as in Fig. 1c and Extended Data Fig. 2a ---------

#Fig. 1c
ggplot(data = pc_df_between, mapping = aes(x = PC1, y = PC4, color=Gender)) + 
  geom_point(alpha = 0.6) + 
  labs(x = paste("BPC1 (", pc_variance_between[1], "% variance)", sep = ''),
       y = paste("BPC4 (", pc_variance_between[4], "% variance)", sep = ''),
       color = "Sex") + 
  scale_color_manual(values = c("Man" = "#9933FF", "Vrouw" = "#FF9933"), 
                     labels = c("Man" = "Male", "Vrouw" = "Female")) +
  scale_x_continuous(position = "top") +
  stat_ellipse() +
  theme_bw() + 
  theme(legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        panel.grid = element_blank())

#Extended Data Fig. 2a
ggplot(data = pc_df_between, mapping = aes(x = PC1, y = PC6, color=as.character(Run))) + 
  geom_point(alpha = 0.6) + 
  labs(x = paste("BPC1 (", pc_variance_between[1], "% variance)", sep = ''),
       y = paste("BPC6 (", pc_variance_between[6], "% variance)", sep = ''),
       color = "Sequencing\nrun") + 
  scale_color_manual(values = c("151" = "#BF616A","152" = "#EBCB8B","125" = "#A3BE8C","8" = "#B48EAD"),
                     labels = c("151" = "Run 1", "152" = "Run 2", "125" = "Run 3", "8" = "Run 4"),
                     breaks = c("151", "152", "125", "8")) +
  scale_x_continuous(position = "bottom") +
  stat_ellipse() +
  theme_bw() + 
  theme(legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        panel.grid = element_blank())


ggplot(data = pc_df_between, mapping = aes(x = PC3, y = PC4, color=Age..years.)) + 
  geom_point(alpha = 0.6) + 
  labs(x = paste("BPC3 (", pc_variance_between[3], "% variance)", sep = ''),
       y = paste("BPC4 (", pc_variance_between[4], "% variance)", sep = ''),
       color = "Age") + 
  scale_color_gradient(low = "#99CCFF", high = "#0000CC") +
  scale_x_continuous(position = "bottom") +
  theme_bw() + 
  theme(legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        panel.grid = element_blank())


ggplot(data = pc_df_between, mapping = aes(x = PC4, y = PC5, color=Time_of_sampling)) + 
  geom_point(alpha = 0.6) + 
  labs(x = paste("BPC4 (", pc_variance_between[4], "% variance)", sep = ''),
       y = paste("BPC5 (", pc_variance_between[5], "% variance)", sep = ''),
       color = "Time of\nsampling") + 
  scale_color_gradient(low = "#FDD0A2", high = "#7F2704") +
  scale_x_continuous(position = "bottom") +
  theme_bw() + 
  theme(legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        panel.grid = element_blank())




## Plot WPCs colored by annual season as in Fig. 1d and Extended Data Fig. 2a ---------

#Fig. 1d
ggplot(data = pc_df_within, mapping = aes(x = PC1, y = PC2, color = Annual_season, fill = Annual_season)) + 
  geom_point(alpha = 0.6) + 
  labs(x = paste("WPC1 (", pc_variance_within[1], "% variance)", sep = ''),
       y = paste("WPC2 (", pc_variance_within[2], "% variance)", sep = ''),
       color = "Annual\nseason",
       fill = "Annual\nseason") + 
  scale_color_manual(values = c("Winter" = '#3399FF', "Autumn" = '#FFB266', "Summer" = '#FF6666', "Spring" = '#FF99CC'), breaks = c("Winter", "Spring", "Summer", "Autumn")) +
  scale_fill_manual(values = c("Winter" = '#3399FF', "Autumn" = '#FFB266', "Summer" = '#FF6666', "Spring" = '#FF99CC'), breaks = c("Winter", "Spring", "Summer", "Autumn")) +
  stat_ellipse() +
  scale_x_continuous(position = "top") +
  theme_bw() + 
  theme(legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        panel.grid = element_blank())

#Extended Data Fig. 2a
ggplot(data = pc_df_within, mapping = aes(x = PC4, y = PC6, color = Annual_season, fill = Annual_season)) + 
  geom_point(alpha = 0.6) + 
  labs(x = paste("WPC4 (", pc_variance_within[4], "% variance)", sep = ''),
       y = paste("WPC6 (", pc_variance_within[6], "% variance)", sep = ''),
       color = "Annual\nseason",
       fill = "Annual\nseason") + 
  scale_color_manual(values = c("Winter" = '#3399FF', "Autumn" = '#FFB266', "Summer" = '#FF6666', "Spring" = '#FF99CC'), breaks = c("Winter", "Spring", "Summer", "Autumn")) +
  scale_fill_manual(values = c("Winter" = '#3399FF', "Autumn" = '#FFB266', "Summer" = '#FF6666', "Spring" = '#FF99CC'), breaks = c("Winter", "Spring", "Summer", "Autumn")) +
  stat_ellipse() +
  scale_x_continuous(position = "bottom") +
  theme_bw() + 
  theme(legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"), 
        panel.grid = element_blank())

