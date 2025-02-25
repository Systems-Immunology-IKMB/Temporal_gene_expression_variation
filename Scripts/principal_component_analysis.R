## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(stringr)
library(dplyr)
library(DESeq2)
library(variancePartition)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


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




## Run PCA ----------------------

pca <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ 1) %>%
  estimateSizeFactors() %>%
  vst() %>%
  assay() %>%
  t() %>%
  prcomp()

all(col_clinical_data$Sample_ID == rownames(pca$x))
pc_df <- cbind(col_clinical_data, pca$x)

pc_variance <- summary(pca) %>%
  .$importance %>%
  .["Proportion of Variance", ] %>% 
  `*`(100) %>%
  round(2) %>%
  sprintf(fmt = '%.2f')




## Plot correlation between variables and PCs as in Fig. 1a ---------

#compute canonical correlation between variables and the first 6 PCs
form <- ~ Individual_ID + Sequencing_Run + Cohort + Annual_season + Time_of_sampling + Age + Sex + BMI + PC1 + PC2 + PC3 + PC4 + PC5 + PC6

Corr_mat <- canCorPairs(form, pc_df) %>%
  .[paste0(rep("PC", 6), 1:6), c("Individual_ID", "Sequencing_Run", "Cohort", "Annual_season", "Time_of_sampling", "Age", "Sex", "BMI")]
colnames(Corr_mat) <- c("Individual", "Sequencing run", "Cohort", "Annual season", "Time of sampling", "Age", "Sex", "BMI")

Corr_mat <- Corr_mat %>%
  as.data.frame() %>%
  add_column("Components" = rownames(.)) %>%
  melt(id.vars = "Components") %>%
  add_column("value_rounded" = round(.$value, digits = 2))

Corr_mat$Components <- factor(Corr_mat$Components, levels = paste0(rep("PC", 6), 1:6) %>% rev())
Corr_mat$variable <- factor(Corr_mat$variable, levels = c("Individual", "Sequencing run", "Cohort", "Annual season", "Time of sampling", "Age", "Sex", "BMI"))

#add grid lines
Corr_mat <- Corr_mat %>%
  mutate(line_positions_y = as.numeric(factor(Components, levels = unique(Components))), 
         line_positions_y= line_positions_y + .5,  
         line_positions_y = ifelse(line_positions_y == max(line_positions_y), NA, line_positions_y)) %>%
  mutate(line_positions_x = as.numeric(factor(variable, levels = unique(variable))), 
         line_positions_x= line_positions_x + .5,  
         line_positions_x = ifelse(line_positions_x == max(line_positions_x), NA, line_positions_x)) 

Corr_plot <- ggplot(Corr_mat, aes(x = variable, y = Components)) +
  geom_point(aes(color = value, fill = value, size = value), shape = 21) +
  geom_text(aes(label = value_rounded), color = "black") +
  geom_vline(aes(xintercept = line_positions_x)) + 
  geom_hline(aes(yintercept = line_positions_y)) + 
  scale_y_discrete(position = "left", 
                   breaks = paste0(rep("PC", 6), 1:6)) + 
  scale_x_discrete(position = "bottom") +
  scale_size(range = c(1, 5), guide = "none") +
  scale_fill_gradient(low = "#F7FBFF", high = "#053061", limits = c(0, 1)) +
  scale_color_gradient(low = "#F7FBFF", high = "#053061", limits = c(0, 1)) +
  labs(fill = "Canonical\ncorrelation", color = "Canonical\ncorrelation", size = "Canonical\ncorrelation", x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank()) +
  coord_fixed()


#barplot for variance explained by first 6 PCs
df_pc_variance <- summary(pc) %>%
  .$importance %>%
  .["Proportion of Variance", ] %>% 
  `*`(100) 
df_pc_variance <- data.frame("Components" = paste0(rep("PC", 6), 1:6), "Variance" = df_pc_variance[1:6]) 
df_pc_variance$Components <- factor(df_pc_variance$Components, levels = paste0(rep("PC", 6), 1:6) %>% rev())

variance_barplot <- ggplot(df_pc_variance, aes(x = Variance, y = Components)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7, fill = "#053061") + 
  scale_y_discrete(position = "left", 
                   breaks = paste0(rep("PC", 6), 1:6)) +
  labs(x = "Variance (%)") +
  scale_x_continuous(breaks = c(0, 15, 30), limits = c(0, 30)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        legend.title = element_blank(),
        axis.ticks.y = element_blank())








## Plot PCs colored by cohort, age,sequencing run and sex as in Fig. 1b ---------

P1 <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC2, col=Cohort)) + 
  geom_point(alpha = 0.6) + xlab("") + 
  ylab(paste("PC2 (", pc_variance[2], "% variance)", sep = '')) 
P1 <- P1 + scale_color_manual(values = c("#7D4B19", "#4B324b"))
P1 <- P1 + stat_ellipse()
P1 <- P1 + theme_bw() + 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        panel.grid = element_blank())

P2 <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC2, color=Sequencing_Run)) + 
  geom_point(alpha = 0.6) + xlab(paste("PC1 (", pc_variance[1], "% variance)", sep = '')) + 
  ylab(paste("PC2 (", pc_variance[2], "% variance)", sep = '')) 
P2 <- P2 + scale_color_manual(values = c("#BF616A","#EBCB8B","#A3BE8C","#B48EAD"))
P2 <- P2 + stat_ellipse()
P2 <- P2 + labs(color = "Sequencing\nRun") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid = element_blank())


P3 <- ggplot(data = pc_df, mapping = aes(x=PC3, y=PC4, color=Age)) + 
  geom_point(alpha = 0.6) + xlab("") + 
  ylab(paste("PC4 (", pc_variance[4], "% variance)", sep = '')) 
P3 <- P3 + scale_color_gradient(low = "#99CCFF", high = "#0000CC")
P3 <- P3 + theme_bw() + 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        panel.grid = element_blank())

P4 <- ggplot(data = pc_df, mapping = aes(x=PC3, y=PC4, color=Sex)) + 
  geom_point(alpha = 0.6) + xlab(paste("PC3 (", pc_variance[3], "% variance)", sep = '')) + 
  ylab(paste("PC4 (", pc_variance[4], "% variance)", sep = '')) 
P4 <- P4 + scale_color_manual(values = c("Male"="#9933FF", "Female"="#FF9933"))
P4 <- P4 + stat_ellipse()
P4 <- P4 + theme_bw() + 
  theme(legend.position = "right",
        panel.grid = element_blank())




## Plot PCs colored by annual season, time of sampling and BMI as in Extended Data Fig. 2a ---------

P1 <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC2, color=Annual_season, fill=Annual_season)) + 
  geom_point(alpha = 0.6) + xlab(paste("PC1 (", pc_variance[1], "% variance)", sep = '')) + 
  ylab(paste("PC2 (", pc_variance[2], "% variance)", sep = '')) 
P1 <- P1 + stat_ellipse()
P1 <- P1 + scale_color_manual(values = c("Winter" = '#3399FF', "Autumn" = '#FFB266', "Summer" = '#FF6666', "Spring" = '#FF99CC'), breaks = c("Winter", "Spring", "Summer", "Autumn"))
P1 <- P1 + scale_fill_manual(values = c("Winter" = '#3399FF', "Autumn" = '#FFB266', "Summer" = '#FF6666', "Spring" = '#FF99CC'), breaks = c("Winter", "Spring", "Summer", "Autumn"))
P1 <- P1 + labs(color = "Annual\nseason", fill = "Annual\nseason")
P1 <- P1 + theme_bw() + 
  theme(legend.position = "right",
        plot.margin = unit(c(0, 0.2, 0.2, 0), "cm"), 
        panel.grid = element_blank())

P2 <- ggplot(data = pc_df, mapping = aes(x=PC3, y=PC4, color=Annual_season, fill=Annual_season)) + 
  geom_point(alpha = 0.6) + xlab(paste("PC3 (", pc_variance[3], "% variance)", sep = '')) + 
  ylab(paste("PC4 (", pc_variance[4], "% variance)", sep = '')) 
P2 <- P2 + stat_ellipse()
P2 <- P2 + scale_color_manual(values = c("Winter" = '#3399FF', "Autumn" = '#FFB266', "Summer" = '#FF6666', "Spring" = '#FF99CC'))
P2 <- P2 + theme_bw() + 
  theme(legend.position = "right",
        plot.margin = unit(c(0, 0, 0.2, 0), "cm"), 
        panel.grid = element_blank())

P3 <- ggplot(data = pc_df[!is.na(pc_df$Time_of_sampling), ], mapping = aes(x=PC3, y=PC4, color=Time_of_sampling)) + 
  geom_point(alpha = 0.6) + xlab(paste("PC3 (", pc_variance[3], "% variance)", sep = '')) + 
  ylab(paste("PC4 (", pc_variance[4], "% variance)", sep = '')) 
P3 <- P3 + scale_color_gradient(low = "#FDD0A2", high = "#7F2704")
P3 <- P3 + labs(color = "Time of\nsampling", fill = "Time of\nsampling")
P3 <- P3 + theme_bw() + 
  theme(legend.position = "right",
        plot.margin = unit(c(0, 0.2, 0, 0), "cm"), 
        panel.grid = element_blank())


P4 <- ggplot(data = pc_df, mapping = aes(x=PC3, y=PC4, color=BMI)) + 
  geom_point(alpha = 0.6) + xlab(paste("PC3 (", pc_variance[3], "% variance)", sep = '')) + 
  ylab(paste("PC4 (", pc_variance[4], "% variance)", sep = '')) 
P4 <- P4 + scale_color_gradient(low = "#C7E9C0", high = "#00441B")
P4 <- P4 + theme_bw() + 
  theme(legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid = element_blank())

