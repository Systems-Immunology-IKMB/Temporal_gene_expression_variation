## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(WGCNA)
library(DESeq2)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(reshape2)
library(lme4)
library(ICC)
library(circlize)
library(ComplexHeatmap)
library(rstatix)
library(ppcor)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


##Loading data ------------------------------------------------------------

gene_list <- read.csv("results/filtered_genes_samples25_fpm5.txt", sep = '\t')
gene_list <- gene_list %>%
  subset(.$gene_type == "protein_coding") %>%
  subset(! .$Chr %in% c("chrX", "chrY"))

count_data <- read.csv("count_files/merged_gene_counts_all_samples.txt", sep = '\t')
colnames(count_data) <- substr(colnames(count_data), 1, 6)
col_clinical_data <- read.csv("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = '\t')

col_clinical_data$Annual_season <- ifelse(col_clinical_data$Month %in% c("January", "February", "March"), "Winter", 
                                          ifelse(col_clinical_data$Month %in% c("July", "August", "September"), "Summer", 
                                                 ifelse(col_clinical_data$Month %in% c("April", "May"), "Spring", "Autumn")))

count_data <- count_data[as.character(gene_list$Gene_ID), as.character(col_clinical_data$Sample_ID)]
colnames(count_data) <- col_clinical_data$Description



## Run WGCNA ----------------

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ Individual_ID)
dds_counts <- estimateSizeFactors(dds_counts)
vst_counts <- vst(dds_counts)
vst_counts <- assay(vst_counts)
vst_counts <- t(vst_counts)

sample_tree <- hclust(dist(vst_counts), method = "average")

pdf("results/Coexpression/output/sample_clustering_to_detect_outliers.pdf", width=25)
par(cex =(0.2));
par(mar = c(0,4,2,0))
plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 0.025,cex.axis = 1.5, cex.main = 2)
abline(h = 95, col = "red")
dev.off()

clust <- cutreeStatic(sample_tree, cutHeight = 95, minSize = 10)
table(clust)

#determine parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(vst_counts, powerVector = powers, verbose = 5)
sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.9;

#Scale-free topology fit index as a function of soft-threshold power: Scale Independence Plot
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90, cool="red")

#Mean connectivity Plot
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Make gene tree
adj_mat <- adjacency(vst_counts, power=9)
TOM = TOMsimilarity(adj_mat)
dissTOM = 1-TOM
TOM_gene_tree <- hclust(as.dist(dissTOM), method = "average")

#Define and merge modules
module_labels <- cutreeDynamicTree(dendro=TOM_gene_tree, minModuleSize = 15, deepSplit = TRUE)
module_colors <- labels2colors(module_labels)
mergedColor <- mergeCloseModules(vst_counts, module_colors, cutHeight = .2)$color

#Save co-expression results
gene_module <- data.frame(Gene = colnames(vst_counts), ModuleLabel = module_labels, ModuleColor = mergedColor)
write.table(gene_module, "results/Coexpression/output/gene_module.txt", quote = FALSE, sep = '\t', row.names = FALSE)

MEsO <- moduleEigengenes(vst_counts, mergedColor)$eigengenes
MEs <- orderMEs(MEsO)
rownames(MEs) <- rownames(vst_counts)
write.table(MEs, "results/Coexpression/output/module_eigengenes.txt", quote = FALSE, sep = '\t')

module_size <- as.data.frame(table(gene_module$ModuleColor))
colnames(module_size) <- c("Module_color", "Number_of_genes")
module_size <- module_size[module_size$Module_color != "grey", ]
module_size <- module_size[order(module_size$Number_of_genes, decreasing = TRUE), ]
module_size$Module_name <- paste0("M", 1:nrow(module_size))
write.table(module_size, "results/Coexpression/output/Module_list.txt", sep = '\t', quote = FALSE)




## Test for association between module eigengenes and clinical parameters --------------------

#import gene co-expression results
res_module <- read.table("results/Coexpression/output/gene_module.txt", header = TRUE) %>%
  merge(annotation, by.x = "Gene", by.y = "Gene_ID")
df_eigengene <- read.table("results/Coexpression/output/module_eigengenes.txt", header = TRUE)
df_module_list <- read.table("results/Coexpression/output/module_list.txt", header = TRUE)


### Adjusted for donor ID -------------------------------

my_modules <- res_module$Module_color %>% unique() %>% setdiff("grey")

df_plot <- df_eigengene %>%
  add_column("Description" = rownames(.)) %>%
  melt(varnames = c("Description")) %>%
  merge(clinical_traits, by = "Description") %>%
  subset(.$variable %in% paste0("ME", my_modules))

#Compare fit of full and reduced model
my_variables <- colnames(clinical_traits)[3:length(colnames(clinical_traits))]
effect_table <- data.frame()
for(m in 1:length(my_modules)){
  for(v in 1:length(my_variables)){
    
    df_test <- df_plot %>%
      subset(.$variable == paste0("ME", my_modules[m])) %>%
      .[, c("Description", "variable", "value", "Individual_ID", my_variables[v])] %>%
      na.omit() 
    
    colnames(df_test)[5] <- "test_var"
    
    if(my_variables[v] %in% c("Neutrophils...µL.", "Lymphocytes...µL.", "Creatinine..mg.dL.", "Annual_season", "Monocytes...µL.", "Hemoglobin..g.dL.", "Age..years.", "Triglycerides..mg.dL.", "Eosinophils...µL.", "Gender", "Thrombocytes..x1000.µL.", "Time_of_sampling", "Albumin..g.L.", "BMI..kg.m2.", "Uric.acid..mg.dL.", "hsCRP..mg.L.")){
      df_test$test_var <- scale(df_test$test_var)
    }
    
    H_1 <- lme(fixed = value ~ test_var, random = ~ 1|Individual_ID, data = df_test, method = "ML")
    H_0 <- lme(fixed = value ~ 1, random = ~ 1|Individual_ID, data = df_test, method = "ML")
    res_table <- anova(H_0, H_1)
    
    effect_table <- rbind(effect_table, 
                          data.frame("Module" = my_modules[m],
                                     "Variable" = my_variables[v],
                                     "no_of_samples" = nrow(df_test),
                                     "coef" = H_1$coefficients$fixed[2],
                                     "L_ratio" = res_table$L.Ratio[2],
                                     "pvalue" = res_table$`p-value`[2]))
  }
}

effect_table <- effect_table %>%
  add_column("padj" = p.adjust(.$pvalue, method = "BH"))

write.table(effect_table, "results/Coexpression/protein_coding_genes/output/module_trait_association.txt", sep = "\t", row.names = FALSE, quote = FALSE)



### Adjusted for donor ID and cell composition -------------------------------

my_modules <- res_module$Module_color %>% unique() %>% setdiff("grey")

df_plot <- df_eigengene %>%
  add_column("Description" = rownames(.)) %>%
  melt(varnames = c("Description")) %>%
  merge(clinical_traits, by = "Description") %>%
  subset(.$variable %in% paste0("ME", my_modules))

#Compare fit of full and reduced model
my_variables <- colnames(clinical_traits)[3:length(colnames(clinical_traits))] %>%
  setdiff(c("Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL."))

effect_table <- data.frame()
m <- 1
v <- 1
for(m in 1:length(my_modules)){
  for(v in 1:length(my_variables)){
    
    df_test <- df_plot %>%
      subset(.$variable == paste0("ME", my_modules[m])) %>%
      .[, c("Description", "variable", "value", "Individual_ID", my_variables[v], "Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL.")] %>%
      na.omit() 
    
    colnames(df_test)[5] <- "test_var"
    
    df_test$test_var <- scale(df_test$test_var)
    df_test$Neutrophils...µL. <- scale(df_test$Neutrophils...µL.)
    df_test$Lymphocytes...µL. <- scale(df_test$Lymphocytes...µL.)
    df_test$Monocytes...µL. <- scale(df_test$Monocytes...µL.)
    df_test$Eosinophils...µL. <- scale(df_test$Eosinophils...µL.)
    df_test$Thrombocytes..x1000.µL. <- scale(df_test$Thrombocytes..x1000.µL.)
    
    H_1 <- lme(fixed = value ~ Neutrophils...µL. + Lymphocytes...µL. + Monocytes...µL. + Eosinophils...µL. + Thrombocytes..x1000.µL. + test_var, random = ~ 1|Individual_ID, data = df_test, method = "ML")
    H_0 <- lme(fixed = value ~ Neutrophils...µL. + Lymphocytes...µL. + Monocytes...µL. + Eosinophils...µL. + Thrombocytes..x1000.µL., random = ~ 1|Individual_ID, data = df_test, method = "ML")
    res_table <- anova(H_0, H_1)
    
    effect_table <- rbind(effect_table, 
                          data.frame("Module" = my_modules[m],
                                     "Variable" = my_variables[v],
                                     "no_of_samples" = nrow(df_test),
                                     "coef" = H_1$coefficients$fixed[7],
                                     "L_ratio" = res_table$L.Ratio[2],
                                     "pvalue" = res_table$`p-value`[2]))
  }
}

effect_table <- effect_table %>%
  add_column("padj" = p.adjust(.$pvalue, method = "BH"))

write.table(effect_table, "results/Coexpression/protein_coding_genes/output/module_trait_association_adjCellComp.txt", sep = "\t", row.names = FALSE, quote = FALSE)





## Compute ICC of module eigengenes and clinical data ------------------------------

#import gene co-expression results
res_module <- read.table("results/Coexpression/output/gene_module.txt", header = TRUE) %>%
  merge(annotation, by.x = "Gene", by.y = "Gene_ID")
df_eigengene <- read.table("results/Coexpression/output/module_eigengenes.txt", header = TRUE)
df_module_list <- read.table("results/Coexpression/output/module_list.txt", header = TRUE)

##compute ICC for clinical variables
var_names <- c("Neutrophils...µL." = "Neutrophils", 
               "Lymphocytes...µL." = "Lymphocytes", 
               "Monocytes...µL." = "Monocytes", 
               "Creatinine..mg.dL." = "Creatinine", 
               "Hemoglobin..g.dL." = "Hemoglobin", 
               "Age..years." = "Age", 
               "Triglycerides..mg.dL." = "Triglycerides", 
               "Eosinophils...µL." = "Eosinophils", 
               "Gender" = "Gender", 
               "Thrombocytes..x1000.µL." = "Thrombocytes", 
               "Albumin..g.L." = "Albumin", 
               "BMI..kg.m2." = "BMI", 
               "hsCRP..mg.L." = "hsCRP", 
               "Uric.acid..mg.dL." = "Uric acid", 
               "Annual_season" = "Annual Season", 
               "Time_of_sampling"  = "Time of sampling")

df_icc <- col_clinical_data[, c("Individual_ID", "Sample_ID", "Neutrophils...µL.", "Lymphocytes...µL.", "Creatinine..mg.dL.", "Monocytes...µL.", "Hemoglobin..g.dL.", "Triglycerides..mg.dL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL.", "Albumin..g.L.", "Uric.acid..mg.dL.", "hsCRP..mg.L.", "Time_of_sampling")] %>%
  melt(id.vars = c("Sample_ID", "Individual_ID")) %>%
  na.omit() %>%
  split(~ .$variable) %>%
  lapply(function(x){res <- ICCest(data = x, Individual_ID, value); return(data.frame("ICC" = res$ICC,
                                                                                      "lower_ICC" = res$LowerCI,
                                                                                      "upper_ICC" = res$UpperCI,
                                                                                      "var_within" = res$varw,
                                                                                      "var_between" = res$vara))}) %>%
  do.call(rbind, .) %>%
  add_column("Variable" = rownames(.), .before = 1)

df_icc$Name <- as_labeller(var_names)(df_icc$Variable) %>% unlist()

write.table(df_icc, "results/Coexpression/output/ICC_clinical_parameters.txt", quote = FALSE, sep = "\t", row.names = FALSE)


##compute ICC for module eigengenes
df_module_list <- df_module_list %>%
  add_column("Variable" = paste0("ME", .$Module_color))

df_icc_modules <- df_eigengene %>%
  add_column("Description" = rownames(.)) %>%
  melt(id.vars = "Description") %>%
  merge(col_clinical_data, by = "Description") %>%
  .[, c("Sample_ID", "Individual_ID", "variable", "value")] %>%
  split(~ .$variable) %>%
  lapply(function(x){res <- ICCest(data = x, Individual_ID, value); return(data.frame("ICC" = res$ICC,
                                                                                      "lower_ICC" = res$LowerCI,
                                                                                      "upper_ICC" = res$UpperCI,
                                                                                      "var_within" = res$varw,
                                                                                      "var_between" = res$vara))}) %>%
  do.call(rbind, .) %>%
  add_column("Variable" = rownames(.), .before = 1) %>%
  merge(df_module_list, by ="Variable")

write.table(df_icc_modules, "results/Coexpression/output/ICC_module_eigengenes.txt", quote = FALSE, sep = "\t", row.names = FALSE)





## Plot association heatmap as in Figure 5b ------------------------------

#load results
res_LMM_adjcellComp <- read.table("results/Coexpression/protein_coding_genes/output/module_trait_association_adjCellComp.txt", header = TRUE)
res_LMM <- read.table("results/Coexpression/protein_coding_genes/output/module_trait_association.txt", header = TRUE)
res_LMM_merged <- res_LMM_adjcellComp %>%
  rbind(res_LMM %>% subset(.$Variable %in% c("Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL."))) %>%
  add_column("padj_combined" = p.adjust(.$pvalue, method = "BH"))

#prepare heatmap input
my_variables_ordered <- c("Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL.", "Gender", "Age..years.", "BMI..kg.m2.", "Annual_season", "Time_of_sampling", "Creatinine..mg.dL.", "Hemoglobin..g.dL.", "Triglycerides..mg.dL.", "Albumin..g.L.", "Uric.acid..mg.dL.", "hsCRP..mg.L.")

res_effect <- res_LMM_merged %>%
  pivot_wider(names_from = "Variable", id_cols = "Module", values_from = "coef") %>%
  .[, c("Module", my_variables_ordered)] %>%
  as.data.frame() 
rownames(res_effect) <- res_effect$Module
res_effect <- res_effect[, -1] 

res_padj <- res_LMM_merged %>%
  subset(.$Variable %in% my_variables_filtered) %>%
  pivot_wider(names_from = "Variable", id_cols = "Module", values_from = "padj_combined") %>%
  .[, c("Module", my_variables_filtered)] %>%
  as.data.frame()
rownames(res_padj) <- res_padj$Module
res_padj <- res_padj[, -1] 

df_adjp <- res_padj %>%
  as.data.frame() %>%
  apply(function(x){case_when(x >= 0.05 ~ "",
                              x >= 0.01 ~ "*",
                              x >= 0.0001 ~ "**",
                              .default = "***")}, MARGIN = 1:2)

#change heatmap orientation
df_adjp <- df_adjp %>%
  t() %>%
  as.data.frame()

res_effect <- res_effect %>%
  t() %>%
  as.data.frame()

#add significance to the heatmap
add_pval <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%s", df_adjp[i, j]), x, y, gp = gpar(fontsize = my_axis_text_size, col = "black"))
}

#add ICC of clinical parameters and module eigengenes
df_icc_clin_param <- read.table("results/Coexpression/output/ICC_clinical_parameters.txt", sep = "\t", header = TRUE)
df_icc_modules <- read.table("results/Coexpression/output/ICC_module_eigengenes.txt", sep = "\t", header = TRUE)

my_icc_param <- df_icc_clin_param %>%
  add_column("my_color" = ifelse(.$ICC > 0.5, "#009900", "#CCCC00")) %>%
  .[match(res_effect %>% rownames(), .$Name), ] 

row_ha <- rowAnnotation("ICC" = anno_barplot(my_icc_param$ICC, 
                                             gp = gpar(fill = my_icc_param$my_color),
                                             width = unit(0.7, "cm"),
                                             border = TRUE,
                                             ylim = c(0, 1),
                                             axis_param = list(
                                               side = "bottom",
                                               at = c(0, 0.5, 1), 
                                               labels = c(" 0.0", "0.5", "1.0"),
                                               labels_rot = 270)),
                        annotation_name_rot = 0,
                        annotation_name_side = "bottom")

my_icc_module <- df_icc_modules %>%
  add_column("my_color" = ifelse(.$ICC > 0.5, "#009900", "#CCCC00")) %>%
  .[match(res_effect %>% colnames(), .$Module_name), ] 

col_ha <- columnAnnotation("ICC" = anno_barplot(my_icc_module$ICC, 
                                                gp = gpar(fill = my_icc_module$my_color),
                                                height = unit(0.7, "cm"),
                                                border = TRUE,
                                                ylim = c(0, 1),
                                                axis_param = list(
                                                  side = "right",
                                                  at = c(0, 0.5, 1), 
                                                  labels = c("0.0", "0.5", "1.0"),
                                                  labels_rot = 270)),
                           annotation_name_rot = 0,
                           annotation_name_side = "right")

#split depending on LMM model
row_split <- factor(c(rep("Cell counts\n(Adjusted\nfor donor)", 5), rep("Clinical parameters\n(Adjusted for\ndonor and cell counts)", 11)))

#plot heatmap and legend
Heatmap(res_effect, 
        name = "LMM\ncoefficient",
        col = colorRamp2(seq(from = -0.025, to = 0.025, by = 0.005), rev(brewer.pal(11, "RdYlBu"))),
        show_heatmap_legend = TRUE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        show_column_names = TRUE,
        column_names_side = "bottom",
        cell_fun = add_pval,
        top_annotation = col_ha, 
        right_annotation = row_ha,
        row_split = row_split,
        row_gap = unit(1, "mm"), 
        border = TRUE, 
        border_gp = gpar(col = "black"),
        heatmap_legend_param = list(direction = "vertical",
                                    title_position = "topleft"),
        width = ncol(res_effect)*unit(3, "mm"), 
        height = nrow(res_effect)*unit(3, "mm"),
        row_title_rot = 90,
        row_title_side = "left",
        cluster_column_slices = FALSE)

Legend(pch = c("*", "**", "***"),
       background = "white",
       border = "lightgrey",
       type = "points", 
       title = "p adj.",
       title_position = "topleft",
       labels_gp = gpar(col = "black"),
       title_gp = gpar(col = "black"),
       labels = c("<0.05", "<0.01", "<0.0001"))




## Plot association heatmap as in Extended Data Fig. 8a ------------------------------

#load results
res_LMM <- read.table("results/Coexpression/protein_coding_genes/output/module_trait_association.txt", header = TRUE)
res_LMM_merged <- res_LMM %>%
  add_column("padj_combined" = p.adjust(.$pvalue, method = "BH"))

#prepare heatmap input
my_variables_ordered <- c("Neutrophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Thrombocytes..x1000.µL.", "Gender", "Age..years.", "BMI..kg.m2.", "Annual_season", "Time_of_sampling", "Creatinine..mg.dL.", "Hemoglobin..g.dL.", "Triglycerides..mg.dL.", "Albumin..g.L.", "Uric.acid..mg.dL.", "hsCRP..mg.L.")

res_effect <- res_LMM_merged %>%
  pivot_wider(names_from = "Variable", id_cols = "Module", values_from = "coef") %>%
  .[, c("Module", my_variables_ordered)] %>%
  as.data.frame() 
rownames(res_effect) <- res_effect$Module
res_effect <- res_effect[, -1] 

res_padj <- res_LMM_merged %>%
  subset(.$Variable %in% my_variables_filtered) %>%
  pivot_wider(names_from = "Variable", id_cols = "Module", values_from = "padj_combined") %>%
  .[, c("Module", my_variables_filtered)] %>%
  as.data.frame()
rownames(res_padj) <- res_padj$Module
res_padj <- res_padj[, -1] 

df_adjp <- res_padj %>%
  as.data.frame() %>%
  apply(function(x){case_when(x >= 0.05 ~ "",
                              x >= 0.01 ~ "*",
                              x >= 0.0001 ~ "**",
                              .default = "***")}, MARGIN = 1:2)

#only show clinical parameters without cell counts
df_adjp <- df_adjp %>%
  .[, -(1:5)]

res_effect <- res_effect %>%
  .[, -(1:5)]

#add significance to the heatmap
add_pval <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%s", df_adjp[i, j]), x, y, gp = gpar(fontsize = my_axis_text_size, col = "black"))
}

#plot heatmap and legend
Heatmap(res_effect, 
        name = "LMM\ncoefficient",
        col = colorRamp2(seq(from = -0.025, to = 0.025, by = 0.005), rev(brewer.pal(11, "RdYlBu"))),
        show_heatmap_legend = TRUE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        show_column_names = TRUE,
        column_names_side = "bottom",
        cell_fun = add_pval,
        row_gap = unit(1, "mm"), 
        border = TRUE, 
        border_gp = gpar(col = "black"),
        heatmap_legend_param = list(direction = "vertical",
                                    title_position = "topleft"),
        width = ncol(res_effect)*unit(3, "mm"), 
        height = nrow(res_effect)*unit(3, "mm"),
        row_title_rot = 90,
        row_title_side = "left",
        cluster_column_slices = FALSE)

Legend(pch = c("*", "**", "***"),
       background = "white",
       border = "lightgrey",
       type = "points", 
       title = "p adj.",
       title_position = "topleft",
       labels_gp = gpar(col = "black"),
       title_gp = gpar(col = "black"),
       labels = c("<0.05", "<0.01", "<0.0001"))




### Code for GO dot plots as in Extended Data Fig. 8c ----------------------------------

direction <- c("M4", "M9", "M10", "M15", "M16")
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
score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = c("M4", "M9", "M10", "M15", "M16"),
                     score_name = score_name,
                     onto = "ALL",
                     reverse_bool = TRUE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 80)

ggplot(df_plot, mapping = aes(x = Module_name, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  scale_colour_gradient(high = "#990000", low = "#FF9999") +
  labs(color = expression(paste("-log"[10], "p")),
       size = paste0("Gene ratio")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = "top")




## Plot module ICC overview as in Extended Data Fig. 8d ------------------------------

df_icc_modules <- read.table("results/Coexpression/output/ICC_module_eigengenes.txt", sep = "\t", header = TRUE)
gene_icc_data <- read.csv("results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t')
res_module <- read.table("results/Coexpression/output/gene_module.txt", header = TRUE) %>%
  merge(annotation, by.x = "Gene", by.y = "Gene_ID")

module_gene_icc_data <- gene_icc_data %>%
  merge(res_module, by = "Gene_ID")
module_names <- unique(module_gene_icc_data$ModuleColor)

ggplot(module_gene_icc_data, aes(x=reorder(Module_name, ICC, FUN = median), y=ICC, fill = ModuleColor)) + 
  geom_boxplot() + 
  geom_point(data = df_icc_modules, mapping = aes(x=Module_name, y=ICC, color = "Module eigengene ICC"), alpha=0.8) + 
  geom_line(data = df_icc_modules, mapping = aes(x=Module_name, y=ICC, group=1, color = "Module eigengene ICC"), alpha=0.8) + 
  geom_hline(yintercept = 0.5, lty=2) + 
  xlab("Module") + ylab("ICC") + 
  scale_fill_manual(values = sort(module_names), guide = "none") + 
  scale_color_manual(name = "", values = c("Module eigengene ICC" = "darkred") ) +
  lims(y = c(0, 1)) +
  theme_bw() + 
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        legend.background = element_blank())


