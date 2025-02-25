## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggmosaic)
library(rstatix)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


## Loading data ------------------------------------------------------------

heritability_scores <- read.csv("info_files/MultiMuTHER_Heritability.txt", sep = '\t')
icc_results <- read.csv("results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t')

icc_heritability <- inner_join(icc_results, heritability_scores, by = c("Gene_ID" = "EnsemblID"))

icc_heritability$Category <- ifelse(icc_heritability$gene_type == "protein_coding", "Protein coding", "Non protein coding")
icc_heritability$Group <- ifelse(icc_heritability$ICC < 0.5, "Within", "Between")
icc_heritability$Heritable <- ifelse(icc_heritability$h2_LowerCI > 0, "YES", "NO")




## Heritability mosaic plot as in Fig. 2i ----------

mosaic_plot_heritability <- ggplot(icc_heritability) + 
  geom_mosaic(aes(x=product(Group, Heritable), fill=Group)) +
  geom_text(data = layer_data(last_plot(), 1) %>% filter(.wt > 0),
            aes(x = (xmin + xmax) / 2,
                y = (ymin + ymax) / 2,
                label = .wt), 
            size = 4) + 
  scale_fill_manual(values = c("Within"='#CCCC00', "Between"='#009900'), name="") + 
  xlab("") + ylab("") + 
  scale_x_productlist(labels = c("Not Heritable", "Heritable")) + 
  scale_y_productlist(labels = c("High inter-\nindividual variance", "High intra-\nindividual variance")) +
  theme_void() +
  theme(axis.text = element_text(),
        axis.title = element_text(),
        axis.text.y = element_text(angle = 90), 
        legend.text = element_text(),
        legend.position = "none",
        aspect.ratio = 1,
        legend.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"))

chisq.test(icc_heritability$Group, icc_heritability$Heritable)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  icc_heritability$Group and icc_heritability$Heritable
#X-squared = 1164.4, df = 1, p-value < 2.2e-16





## Heritability violin plot as in Extended Data Fig. 4d ----------

Heritability_wt <- icc_heritability %>%
  wilcox_test(h2 ~ Group, data = . , paired = FALSE) %>%
  add_significance(p.col = "p", 
                   output.col = "p_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = 0.95)

violin_plot <- ggplot(icc_heritability) + 
  geom_violin(alpha = 0.5, mapping = aes(x=Group, y=h2, fill = Group)) + 
  geom_boxplot(width = 0.1, mapping = aes(x=Group, y=h2, fill = Group)) + 
  xlab("") + 
  ylab("Heritability") + 
  scale_fill_manual(values = c("Within"='#CCCC00', "Between"='#009900'), name="", 
                    labels = c("Between"="High inter-individual\nvariance genes", "Within"="High intra-individual\nvariance genes")) + 
  stat_pvalue_manual(data = Heritability_wt, y.position = "max", label = "p_signif", tip.length = 0, bracket.shorten = 0.5) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())



