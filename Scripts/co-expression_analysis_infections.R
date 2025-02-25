## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(rstatix)
library(DESeq2)
library(ggpubr)
library(DBI)
library(topGO)
library(plyr)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")
source("results/Enrichment_analysis_helper_functions.R")


##Loading data ------------------------------------------------------------

#import metadata
col_clinical_data <- read.table("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = "\t", header = TRUE)
col_clinical_data$Individual_ID <- factor(col_clinical_data$Individual_ID)
rownames(col_clinical_data) <- col_clinical_data$Sample_ID
col_clinical_data$Run <- factor(col_clinical_data$Run)

gene_attributes <- read.table("info_files/filtered_genes_samples25_fpm5.txt", sep = "\t", header = TRUE)
annotation <- gene_attributes[, c("Gene_ID", "gene_name")] %>% unique()
rownames(annotation) <- annotation$Gene_ID
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



## Select samples with acute infections and matched controls --------------------------

### Select samples with viral infections and matched controls ---------
infected_individuals <- col_clinical_data %>%
  subset(.$Infection_now == "YES") %>%
  .$Individual_ID

samples_viral_infection <- col_clinical_data %>%
  subset(.$Infection_type == "Viral") %>%
  subset(! .$Which_infection %in% c("HERPES_GENITALIS", "HERPES")) %>%
  add_column("group" = "infected")

viral_individuals <- samples_viral_infection %>%
  .$Individual_ID %>%
  unique()

length(viral_individuals) #14

#get controls matched for age, sex and annual season
set.seed(1234)

control_samples <- c()
i <- 3
for(i in 1:length(samples_viral_infection$Sample_ID)){
  
  my_samples <- col_clinical_data %>%
    subset(! .$Individual_ID %in% infected_individuals) %>%
    subset(.$Infection_now == "NO") %>%
    subset(.$Cold_last_month == "NO") %>%
    subset(.$healthy == "YES") %>%
    subset(.$Gender == samples_viral_infection[i, "Gender"]) %>%
    subset(.$Age..years. == samples_viral_infection[i, "Age..years."]) %>%
    subset(.$Annual_season == samples_viral_infection[i, "Annual_season"]) %>%
    subset(! .$Sample_ID %in% control_samples)
  
  if(nrow(my_samples) == 0){
    my_samples <- col_clinical_data %>%
      subset(! .$Individual_ID %in% infected_individuals) %>%
      subset(.$Infection_now == "NO") %>%
      subset(.$Cold_last_month == "NO") %>%
      subset(.$healthy == "YES") %>%
      subset(.$Gender == samples_viral_infection[i, "Gender"]) %>%
      subset(.$Age..years. > as.numeric(samples_viral_infection[i, "Age..years."]) - 5) %>%
      subset(.$Age..years. < as.numeric(samples_viral_infection[i, "Age..years."]) + 5) %>%
      subset(.$Annual_season == samples_viral_infection[i, "Annual_season"]) %>%
      subset(! .$Sample_ID %in% control_samples)
  }
  
  selected_sample <- my_samples[sample(1:nrow(my_samples), 1, replace = FALSE), "Sample_ID"]
  control_samples <- c(control_samples, selected_sample)
}

control_samples %>% unique() %>% length() #14

samples_viral_control <- col_clinical_data %>% 
  subset(.$Sample_ID %in% control_samples) %>%
  .[match( control_samples, .$Sample_ID), ] %>%
  add_column("group" = "control")

col_data_viral <- rbind(samples_viral_infection, samples_viral_control) 
write.table(col_data_viral, "results/infections/output/samples_viral_infections_healthy.txt", row.names = FALSE, quote = FALSE, sep = "\t")






### Get samples with bacterial infections and matched controls -----------------
infected_individuals <- col_clinical_data %>%
  subset(.$Infection_now == "YES") %>%
  .$Individual_ID

samples_bacterial_infection <- col_clinical_data %>%
  subset(.$Infection_type == "Bacterial")  %>%
  add_column("group" = "infected")

bacterial_individuals <- samples_bacterial_infection %>%
  .$Individual_ID %>%
  unique()

length(bacterial_individuals) #15

#get controls matched for age, sex and annual season
set.seed(1234)

control_samples <- c()
i <- 3
for(i in 1:length(samples_bacterial_infection$Sample_ID)){
  
  my_samples <- col_clinical_data %>%
    subset(! .$Individual_ID %in% infected_individuals) %>%
    subset(.$Infection_now == "NO") %>%
    subset(.$Cold_last_month == "NO") %>%
    subset(.$healthy == "YES") %>%
    subset(.$Gender == samples_bacterial_infection[i, "Gender"]) %>%
    subset(.$Age..years. == samples_bacterial_infection[i, "Age..years."]) %>%
    subset(.$Annual_season == samples_bacterial_infection[i, "Annual_season"]) %>%
    subset(! .$Sample_ID %in% control_samples)
  
  if(nrow(my_samples) == 0){
    my_samples <- col_clinical_data %>%
      subset(! .$Individual_ID %in% infected_individuals) %>%
      subset(.$Infection_now == "NO") %>%
      subset(.$Cold_last_month == "NO") %>%
      subset(.$healthy == "YES") %>%
      subset(.$Gender == samples_bacterial_infection[i, "Gender"]) %>%
      subset(.$Age..years. > as.numeric(samples_bacterial_infection[i, "Age..years."]) - 5) %>%
      subset(.$Age..years. < as.numeric(samples_bacterial_infection[i, "Age..years."]) + 5) %>%
      subset(.$Annual_season == samples_bacterial_infection[i, "Annual_season"]) %>%
      subset(! .$Sample_ID %in% control_samples)
  }
  
  selected_sample <- my_samples[sample(1:nrow(my_samples), 1, replace = FALSE), "Sample_ID"]
  control_samples <- c(control_samples, selected_sample)
}

control_samples %>% unique() %>% length() #18

samples_bacterial_control <- col_clinical_data %>% 
  subset(.$Sample_ID %in% control_samples) %>%
  .[match( control_samples, .$Sample_ID), ] %>%
  add_column("group" = "control")

col_data_bacterial <- rbind(samples_bacterial_infection, samples_bacterial_control) 
write.table(col_data_bacterial, "results/infections/output/samples_bacterial_infections_healthy.txt", row.names = FALSE, quote = FALSE, sep = "\t")





## Get remaining longitudinal samples from participants with infections ---------------------------------------------------

infection <- c("Viral", "Bacterial")
df_infection <- data.frame()
j <- 1
t <- 1
for(t in 1:length(infection)){
  
  if(infection[t] == "Viral"){
    samples_infected_control <- read.table("results/infections/output/samples_viral_infections_healthy.txt", sep = "\t", header = TRUE)
  }else if(infection[t] == "Bacterial"){
    samples_infected_control <- read.table("results/infections/output/samples_bacterial_infections_healthy.txt", sep = "\t", header = TRUE)
  }
  
  col_data_infection <- col_clinical_data %>%
    subset(.$Individual_ID %in% (samples_infected_control %>% subset(.$Infection_type == infection[t]) %>% .$Individual_ID) ) %>%
    subset(! is.na(.$Infection_now))
  
  col_data_infection$Individual_ID <- col_data_infection$Individual_ID %>% droplevels()
  
  #get before, infected or after samples
  my_individuals <- col_data_infection$Individual_ID %>% unique()
  df_group <- data.frame()
  
  i <- 1
  for(i in 1:length(my_individuals)){
    
    my_individual <- my_individuals[i]
    col_data_temp <- col_data_infection %>%
      subset(.$Individual_ID == my_individual)
    
    my_samples <- col_data_temp$Sample_ID
    
    for(j in 1:length(my_samples)){
      
      my_timepoint <- col_data_temp[col_data_temp$Sample_ID == my_samples[j], ]$Timepoint
      infected_timepoint <- col_data_temp[col_data_temp$Infection_type == infection[t], ]$Timepoint %>% max()
      
      col_data_temp_after <- col_data_temp %>%
        subset(.$Timepoint < my_timepoint) %>%
        subset(.$Infection_type == infection[t])
      
      col_data_temp_after6 <- col_data_temp %>%
        subset(my_timepoint > infected_timepoint + 3)
      
      col_data_temp_before <- col_data_temp %>%
        subset(.$Timepoint > my_timepoint) %>%
        subset(.$Infection_type == infection[t])
      
      if(col_data_temp[col_data_temp$Sample_ID == my_samples[j], ]$Infection_type == infection[t]){
        df_group <- rbind(df_group, data.frame("Sample_ID" = my_samples[j],
                                               "group" = "infected"))
      }else if(col_data_temp[col_data_temp$Sample_ID == my_samples[j], ]$Infection_now == "NO" &
               nrow(col_data_temp_before) > 0){
        df_group <- rbind(df_group, data.frame("Sample_ID" = my_samples[j],
                                               "group" = "before"))
      }else if(col_data_temp[col_data_temp$Sample_ID == my_samples[j], ]$Infection_now == "NO" &
               nrow(col_data_temp_after6) > 0){
        df_group <- rbind(df_group, data.frame("Sample_ID" = my_samples[j],
                                               "group" = "after6"))
      }else if(col_data_temp[col_data_temp$Sample_ID == my_samples[j], ]$Infection_now == "NO" &
               nrow(col_data_temp_after) > 0){
        df_group <- rbind(df_group, data.frame("Sample_ID" = my_samples[j],
                                               "group" = "after3"))
      }
    }
  }
  
  col_data_infection <- col_data_infection %>%
    merge(df_group, by = "Sample_ID")
  
  col_data_infection <- col_data_infection %>%
    rbind(samples_infected_control[, 1:143] %>% subset(.$Infection_type != infection[t]) %>% add_column("group" = "control")) %>%
    add_column("Infection_group" = infection[t])

  col_data_infection$group <- factor(col_data_infection$group, 
                                     levels = c("control", "before", "infected", "after3", "after6"))
  
  df_infection <- rbind(df_infection, col_data_infection)
}

write.table(df_infection, "results/infections/output/metadata_infections_healthy.txt", sep = "\t", quote = FALSE, row.names = FALSE)

table(df_infection$group, df_infection$Infection_group)
#         Bacterial Viral
#control         18    14
#before          13     5
#infected        18    14
#after3           8    11
#after6           5     8






## Code to plot module analysis results as in Figure 5 and Extended Data Fig. 8 ---------

### Plot eigengene heatmap as in Figure 5e ------------------------------

df_module_list <- df_module_list %>%
  add_column("variable" = paste0("ME", .$Module_color))

df_infection <- read.table("results/infections/output/metadata_infections_healthy.txt", sep = "\t", header = TRUE) %>%
  subset(.$group != "before") 

infection_modules <- c("palevioletred1", "skyblue4", "tan4", "lavenderblush1", "brown2", "firebrick4", "lightsteelblue1")

df_plot <- df_eigengene %>%
  add_column("Description" = rownames(.)) %>%
  melt(varnames = c("Description")) %>%
  merge(df_infection, by = "Description") %>%
  subset(.$variable %in% paste0("ME", infection_modules)) %>%
  merge(df_module_list, by = "variable")

df_mean <- df_plot %>%
  group_by(Module_name, group, Infection_group) %>%
  summarize("mean" = mean(value)) %>%
  add_column("my_order" = paste0(.$Infection_group, "_", .$group)) %>%
  .[, c(-2, -3)]

matDataHeat <- df_mean %>%
  pivot_wider(id_cols = "Module_name", names_from = "my_order", values_from = "mean") %>%
  as.data.frame() %>%
  .[, c("Module_name", paste0(rep(c("Bacterial_", "Viral_"), each = 4), rep(c("control", "infected", "after3", "after6"), 2)))]
rownames(matDataHeat) <- matDataHeat$Module_name
matDataHeat <- matDataHeat[, -1]

col_anno <- data.frame("Infection type" = c(rep(c("Bacterial", "Viral"), each = 4)),
                       "Infection status" = c("control", "infected_bacterial", "after3", "after6", "control", "infected_viral", "after3", "after6"), 
                       check.names = FALSE)

color_infection <- c("control" = "grey", "infected_bacterial" = "#CC79A7", "infected_viral" = "#E69F00", "after3" = "#56B4E9", "after6" = "#0072B2")
color_type <- c("Bacterial" = "#CC79A7", "Viral" = "#E69F00")

col_ha <- HeatmapAnnotation(df = col_anno,
                            col = list("Infection type" = color_type,
                                       "Infection status" = color_infection),
                            annotation_legend_param = list("Infection status" = 
                                                             list(labels = c("Control", "Bacterial infection", "Viral infection", "3 months after", "6 months after"),
                                                                  at = c("control", "infected_bacterial", "infected_viral", "after3", "after6"),
                                                                  legend_gp = gpar(fill = c("grey", "#CC79A7", "#E69F00", "#56B4E9", "#0072B2")),
                                                                  ncol = 2,
                                                                  border_gp = gpar(col = "black"))),
                            border = TRUE,
                            show_annotation_name = FALSE,
                            show_legend = c(FALSE, TRUE))

col_split <- c(rep("A", 4), rep("B", 4)) %>%
  factor()

my_icc <- df_icc_modules %>%
  add_column("my_color" = ifelse(.$ICC > 0.5, "between", "within")) %>%
  .[match(matDataHeat %>% rownames(), .$Module_name), ] 

color_icc <- c("within" = "#CCCC00", "between" = "#009900")

row_ha <- rowAnnotation("ICC" = anno_simple(my_icc$my_color, 
                                            col = color_icc,
                                            border = TRUE),
                        module_names = anno_text(matDataHeat %>% rownames()),
                        annotation_name_rot = 0,
                        annotation_name_side = "top",
                        border =  gpar(col = "black"),
                        gp = gpar(col = "black"))

ht_infection_eigengenes <- Heatmap(matDataHeat, 
                                   name = "Eigengene",
                                   col = colorRamp2(seq(from = -0.03, to = 0.03, by = 0.006), rev(brewer.pal(11, "RdYlBu"))),
                                   cluster_columns = FALSE,
                                   cluster_rows = TRUE,
                                   top_annotation = col_ha, 
                                   right_annotation = row_ha,
                                   row_gap = unit(0.5, "mm"), 
                                   border = TRUE, 
                                   border_gp =  gpar(col = "black"),
                                   show_row_names = FALSE,
                                   show_column_names = FALSE,
                                   column_split = col_split,
                                   row_title_rot = 0,
                                   heatmap_legend_param = list(direction = "horizontal",
                                                               title_position = "lefttop"),
                                   row_dend_width = unit(0.25, "cm"),
                                   height = unit(1, "in"),   
                                   cluster_row_slices = FALSE,
                                   column_title = NULL,
                                   show_column_dend = TRUE,
                                   cluster_column_slices = FALSE)
draw(ht_infection_eigengenes, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")



### Plot heatmap of module genes as in Figure 5g ------------------------------

res_icc <- read.table("results/ICC/output/gene_icc_data_whole_dataset.txt", header = TRUE, sep = "\t")
df_infection <- read.table("results/infections/output/metadata_infections_healthy.txt", sep = "\t", header = TRUE) %>%
  subset(.$group != "before")

#Get genes with highest membership score for each module
my_modules <- c("skyblue4", "lightsteelblue1", "tan4")
selected_genes <- c()
for(i in 1:length(my_modules)){
  
  my_genes <- res_module %>%
    subset(.$ModuleColor == my_modules[i]) %>%
    .$Gene
  
  selected_genes_temp <- df_membership %>%
    add_column("Gene_ID" = rownames(.)) %>%
    .[my_genes, c("Gene_ID", paste0("kME", my_modules[i]))] %>%
    .[order(.[[paste0("kME", my_modules[i])]], decreasing = TRUE), ] %>%
    subset(.[[paste0("kME", my_modules[i])]] >= 0.8) %>%
    .$Gene_ID
  
  selected_genes <- c(selected_genes, selected_genes_temp)
}

my_module <- res_module %>%
  subset(.$ModuleColor %in% c("skyblue4", "lightsteelblue1", "tan4")) %>%
  mutate(ModuleColor = factor(ModuleColor, levels = c("skyblue4", "lightsteelblue1", "tan4"))) %>%
  .[order(.$ModuleColor), ]

my_genes <- my_module %>%
  .$Gene

df_plot <- normalized_counts[my_genes, ] %>%
  t() %>%
  as.data.frame() %>%
  add_column("Sample_ID" = rownames(.)) %>%
  melt(varnames = c("Sample_ID")) %>%
  merge(df_infection, by = "Sample_ID") 

df_mean <- df_plot %>%
  group_by(variable, group, Infection_group) %>%
  summarize("mean" = mean(value)) %>%
  add_column("my_order" = paste0(.$Infection_group, "_", .$group)) %>%
  .[, c(-2, -3)]

matDataHeat <- df_mean %>%
  pivot_wider(id_cols = "variable", names_from = "my_order", values_from = "mean") %>%
  as.data.frame() %>%
  .[, c("variable", paste0(rep(c("Bacterial_", "Viral_"), each = 4), rep(c("control", "infected", "after3", "after6"), 2)))]
rownames(matDataHeat) <- matDataHeat$variable
matDataHeat <- matDataHeat[my_genes, -1]
matDataHeat <- scale(t(matDataHeat))
matDataHeat <- matDataHeat[rownames(matDataHeat) %>% rev(), ]

row_anno <- data.frame("Infection type" = c(rep(c("Bacterial", "Viral"), each = 4)) %>% rev(),
                       "Infection status" = c("control", "infected_bacterial", "after3", "after6", "control", "infected_viral", "after3", "after6") %>% rev(), 
                       check.names = FALSE)

color_infection <- c("control" = "grey", "infected_bacterial" = "#CC79A7", "infected_viral" = "#E69F00", "after3" = "#56B4E9", "after6" = "#0072B2")
color_type <- c("Bacterial" = "#CC79A7", "Viral" = "#E69F00")

row_ha <- HeatmapAnnotation(df = row_anno,
                            col = list("Infection type" = color_type,
                                       "Infection status" = color_infection),
                            annotation_legend_param = list("Infection status" = 
                                                             list(labels = c("Control", "Bacterial infection", "Viral infection", "After 3", "After 6"),
                                                                  at = c("control", "infected_bacterial", "infected_viral", "after3", "after6"),
                                                                  legend_gp = gpar(fill = c("grey", "#CC79A7", "#E69F00", "#56B4E9", "#0072B2")))),
                            border = TRUE,
                            simple_anno_size = unit(0.2, "cm"),
                            annotation_name_gp = gpar(size = 7),
                            show_legend = FALSE, 
                            show_annotation_name = FALSE,
                            name = "",
                            which = "row")

row_split <- c(rep("A", 4), rep("B", 4)) %>%
  factor()

my_module <- res_module %>%
  subset(.$ModuleColor %in% c("skyblue4", "lightsteelblue1", "tan4")) %>%
  mutate(ModuleColor = factor(ModuleColor, levels = c("skyblue4", "lightsteelblue1", "tan4"))) %>%
  .[order(.$ModuleColor), ] %>%
  merge(df_module_list, by.x = "ModuleColor", by.y = "Module_color", sort = FALSE)

col_split <- my_module$Module_name %>%
  factor()

my_names <- matDataHeat %>% 
  colnames() %>%
  data.frame("Gene_ID"= . ) %>%
  merge(annotation, by = "Gene_ID", sort = FALSE) %>%
  .[match(matDataHeat %>% colnames(), .$Gene_ID), ] 

my_idx <- which(my_names$Gene_ID %in% selected_genes)
my_names_idx <- my_names[my_idx, "gene_name"]

my_icc <- matDataHeat %>% 
  colnames() %>%
  data.frame("Gene_ID"= . ) %>%
  merge(res_icc, by = "Gene_ID", sort = FALSE) %>%
  add_column("my_color" = ifelse(.$ICC > 0.5, "between", "within")) %>%
  .[match(matDataHeat %>% colnames(), .$Gene_ID), ] 

color_icc <- c("within" = "#CCCC00", "between" = "#009900")

col_ha <- HeatmapAnnotation("ICC" = anno_simple(my_icc$my_color, 
                                                col = color_icc,
                                                simple_anno_size = unit(0.2, "cm"),
                                                border = TRUE),
                            gene_marks = anno_mark(at = my_idx,
                                                   labels = my_names_idx,
                                                   padding = unit(0.2, "mm"),
                                                   which = "column", 
                                                   link_width = unit(0.35, "cm"),
                                                   side = "bottom"),
                            annotation_name_rot = 0,
                            annotation_name_side = "left",
                            gp = gpar(col = "black"))

ht_infection_genes <- Heatmap(matDataHeat, 
                              name = "Scaled normalized expression",
                              col = colorRamp2(seq(from = -2, to = 2, by = 0.4), rev(brewer.pal(11, "RdYlBu"))),
                              cluster_columns = TRUE,
                              cluster_rows = FALSE,
                              bottom_annotation = col_ha, 
                              left_annotation = row_ha,
                              column_dend_side = "top",
                              row_gap = unit(0.7, "mm"), 
                              border = TRUE, 
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              column_split = col_split,
                              row_split = row_split,
                              column_title_rot = 0,
                              column_title_side = "top",
                              cluster_row_slices = FALSE,
                              row_title = NULL,
                              column_gap = unit(0.7, "mm"),
                              heatmap_legend_param = list(direction = "horizontal",
                                                          title_position = "lefttop"),
                              column_dend_height = unit(0.5, "cm"),
                              show_column_dend = TRUE,
                              show_row_dend = FALSE,
                              cluster_column_slices = FALSE)

lgd_icc <- Legend(title = "ICC", 
                  legend_gp = gpar(fill = color_icc),
                  at = c("within", "between"), 
                  labels = c("<0.5", ">0.5"),
                  direction = "horizontal",
                  title_position = "leftcenter",
                  border =  gpar(col = "black"),
                  nrow = 1,
                  labels_gp = gpar(col = "black"),
                  title_gp = gpar(col = "black"))

draw(ht_infection_genes, 
     heatmap_legend_side = "top", 
     annotation_legend_side = "top", 
     merge_legends = TRUE, 
     annotation_legend_list = list(lgd_icc),
     newpage = FALSE)



### Code for sample overview as in Extended Data Fig. 8e ------------------------------

df_infection <- read.table("results/infections/output/metadata_infections_healthy.txt", sep = "\t", header = TRUE)

my_colors <- c("before" = "black", "infected_Bacterial" = "#CC79A7", "infected_Viral" = "#E69F00", "after3" = "#56B4E9", "after6" = "#0072B2", "Control" = "grey")
my_fills <- c("before" = "white", "infected_Bacterial" = "#CC79A7", "infected_Viral" = "#E69F00", "after3" = "#56B4E9", "after6" = "#0072B2", "Control" = "grey")

df_plot <- df_infection %>%
  subset(.$group != "control") %>%
  add_column("group_plot" = ifelse(.$group == "infected", paste0("infected_", .$Infection_group), .$group))

df_plot$group_plot <- factor(df_plot$group_plot, levels = names(my_colors))

df_plot$Infection_type <- factor(df_plot$Infection_type, 
                                 levels = c("", "Viral", "Bacterial") %>% rev(),
                                 labels = c("no", "Viral", "Bacterial") %>% rev())

df_plot <- df_plot %>%
  .[order(.$Infection_type, .$Timepoint, decreasing = FALSE), ] 

df_plot$Individual_ID <- factor(df_plot$Individual_ID, levels = df_plot$Individual_ID %>% unique() %>% rev())


#Cohort A
df_plot$Annual_season <- factor(df_plot$Annual_season, levels = c("Winter", "Spring", "Summer", "Autumn"))

df_infection_A <- df_plot %>%
  subset(.$Season == "A")

samples_A <- ggplot(df_infection_A, aes(x = Annual_season, y = Individual_ID, color = group_plot, fill = group_plot)) +
  geom_point(size = 3.5, shape = 21) +
  scale_color_manual(values = my_colors, 
                     breaks = c("before", "infected_Bacterial", "infected_Viral", "after3", "after6", "Control"),
                     labels = c("Before infection", "Bacterial infection", "Viral infection", "After 3", "After 6", "Control"),
                     drop = FALSE) +
  scale_fill_manual(values = my_fills, 
                    breaks = c("before", "infected_Bacterial", "infected_Viral", "after3", "after6", "Control"),
                    labels = c("Before infection", "Bacterial infection", "Viral infection", "After 3", "After 6", "Control"),
                    drop = FALSE) +
  labs(y = "Individual", color = "Infection type", title = "Cohort A", fill = "Infection type", ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


#Cohort B
df_plot$Annual_season <- factor(df_plot$Annual_season, levels = c("Summer", "Autumn", "Winter", "Spring"))

df_infection_B<- df_plot %>%
  subset(.$Season == "B")

samples_B <- ggplot(df_infection_B, aes(x = Annual_season, y = Individual_ID, color = group_plot, fill = group_plot)) +
  geom_point(size = 3.5, shape = 21) +
  scale_color_manual(values = my_colors, 
                     breaks = c("before", "infected_Bacterial", "infected_Viral", "after3", "after6", "Control"),
                     labels = c("Before infection", "Bacterial infection", "Viral infection", "After 3", "After 6", "Control"),
                     drop = FALSE) +
  scale_fill_manual(values = my_fills, 
                    breaks = c("before", "infected_Bacterial", "infected_Viral", "after3", "after6", "Control"),
                    labels = c("Before infection", "Bacterial infection", "Viral infection", "After 3", "After 6", "Control"),
                    drop = FALSE) +
  labs(y = "Individual", color = "Infection type", title = "Cohort B", fill = "Infection type", ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_B_and_legend <- (guide_area()/samples_B) + plot_layout(nrow = 2, heights = c(0.82, 0.18), guides = 'collect')
plot(samples_A + p_B_and_legend + plot_layout(ncol = 2))





### Cell count and hsCRP barplots as in Extended Data Fig. 8f ------------------------------

df_infection <- read.table("results/infections/output/metadata_infections_healthy.txt", sep = "\t", header = TRUE) %>%
  subset(.$group != "before")

var_names <- c("Neutrophils...µL." = "Neutrophils", 
               "Lymphocytes...µL." = "Lymphocytes", 
               "Monocytes...µL." = "Monocytes", 
               "Thrombocytes..x1000.µL." = "Thrombocytes", 
               "hsCRP..mg.L." = "hsCRP", 
               "Viral" = "Viral",
               "Bacterial" = "Bacterial")

my_colors <- c("control" = "grey", "infected_Viral" = "#E69F00", "after3" = "#56B4E9", "after6" = "#0072B2", "infected_Bacterial" = "#CC79A7")

df_plot <- df_infection[, c("Individual_ID", "Sample_ID", "group", "Infection_group", "Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L.")] %>%
  melt(id.vars = c("Sample_ID", "group", "Infection_group", "Individual_ID")) %>%
  na.omit() %>%
  add_column("group_plot" = ifelse(.$group == "infected", paste0("infected_", .$Infection_group), .$group))

df_plot$my_order <- factor(df_plot$group, levels = c("control", "infected", "after3", "after6")) %>%
  as.numeric()

df_plot$variable <- factor(df_plot$variable, levels = c("Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L."))

df_test <- df_infection[, c("Sample_ID", "group", "Infection_group", "Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L.")] %>%
  melt(id.vars = c("Sample_ID", "group", "Infection_group")) %>%
  group_by(variable, Infection_group) %>%
  kruskal_test(value ~ group, data = .) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("max" = df_infection[, c("Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L.")] %>%
               apply(max, na.rm = TRUE, MARGIN = 2) %>%
               '*' (0.85) %>%
               rep(each = 2)) %>%
  add_column("group1" = 2,
             "group2" = 3)

p_cell_counts <- ggplot(df_plot, aes(x = my_order)) +
  geom_boxplot(aes(x = my_order, y = value, group = group_plot, color = group_plot, fill = group_plot), 
               width = 0.7, alpha = 0.5, linewidth = 0.3, outlier.size = 0.8) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(labels = c("Summer", "Winter", "Test", "After 6"), breaks = 1:4) +
  facet_grid(variable ~ Infection_group, scales = "free_y", labeller = as_labeller(var_names), switch = "y") +
  stat_pvalue_manual(data = df_test, label.x.npc = "middle", y.position = "max", label = "p.adj_signif", vjust = 0, tip.length = 0, bracket.shorten = 0.5, size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.placement = "outside",
        legend.position = "none",
        panel.border = element_rect(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())






### Code for GO dot plots as in Figure 5f and Extended Data Fig. 8g ----------------------------------

direction <- c("palevioletred1", "skyblue4", "tan4", "lavenderblush1", "brown2", "firebrick4", "lightsteelblue1")
GO_tables <- list()
for(d in 1:length(direction)){
  
  table_temp <- read.table(paste0("results/Coexpression/output/Module_", direction[d], "_GO.txt"), sep = "\t", header = TRUE)
  table_temp <- table_temp %>%
    mutate(Fisher.elim = Fisher.elim %>% str_remove("< ") %>% as.numeric()) %>%
    add_column("Ratio" = .$Significant/.$Annotated,
               "changed_scores" = -log10(.$Fisher.elim)) %>%
    mutate(Term = getTermsDefinition(whichTerms = GO.ID, "BP"))
  GO_tables <- append(GO_tables,
                      list(table_temp))
}
names(GO_tables) <- direction


#GO dotplot as in Figure 5f
GO_tables_figure_5f <- GO_tables[c("skyblue4", "lightsteelblue1", "tan4")]

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables_figure_5f,
                     direction = c("skyblue4", "lightsteelblue1", "tan4"),
                     score_name = score_name,
                     onto = "ALL",
                     reverse_bool = TRUE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 80)

df_plot <- df_plot %>%
  merge(df_module_list, by.x = "Direction", by.y = "Module_color")
df_plot$Module_name <- factor(df_plot$Module_name, levels = c("M12", "M15", "M21"))

p_GO <- ggplot(df_plot, mapping = aes(x = Module_name, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  scale_colour_gradient(high = "#990000", low = "#FF9999") +
  labs(color = expression(paste("-log"[10], "p")),
       size = paste0("Gene ratio")) +
  theme_bw() +
  theme(legend.title = element_text(size = 10),
        axis.title = element_blank(),
        legend.position = "right")


#GO dotplot as in Extended Data Fig. 7g
GO_tables_figure_E7g <- GO_tables[c("brown2", "firebrick4", "palevioletred1", "lavenderblush1")]

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables_figure_E7g,
                     direction = c("brown2", "firebrick4", "palevioletred1", "lavenderblush1"),
                     score_name = score_name,
                     onto = "ALL",
                     reverse_bool = TRUE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 80)

df_plot <- df_plot %>%
  merge(df_module_list, by.x = "Direction", by.y = "Module_color")
df_plot$Module_name <- factor(df_plot$Module_name, levels = c("M2", "M18", "M22", "M23"))

p_GO <- ggplot(df_plot, mapping = aes(x = Module_name, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  scale_colour_gradient(high = "#990000", low = "#FF9999") +
  labs(color = expression(paste("-log"[10], "p")),
       size = paste0("Gene ratio")) +
  theme_bw() +
  theme(legend.title = element_text(size = 10),
        axis.title = element_blank(),
        legend.position = "right")



