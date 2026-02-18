## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape)
library(Seurat)
library(ggpubr)
library(ggmosaic)
library(genefilter)
library(geneplotter)
library(topGO)
library(DBI)
library(plyr)
library(clusterProfiler)
library(msigdbr)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")
source("results/enrichment_analysis_helper_functions.R")


## Loading data -----------------------------------

#scObject from Schulte-Schrepping et al. - Severe COVID-19 Is Marked by a Dysregulated Myeloid Cell Compartment (https://doi.org:10.1016/j.cell.2020.08.001)
single_cell_object <- readRDS("count_files/seurat_COVID19_freshWB-PBMC_cohort2_rhapsody_jonas_FG_2020-08-18.rds")
single_cell_object <- SetIdent(single_cell_object, value = single_cell_object@meta.data$cluster_labels_res.0.8)

#keep only the samples from healthy individuals with no comorbidities
Object_V1 <- subset(x = single_cell_object, 
                    cells = colnames(single_cell_object[, single_cell_object$comorbidities == "none" & single_cell_object$diagnosis == "control"]))
Object_V1 <- subset(x = Object_V1, idents = c("Mixed_cells", "CD34+ GATA2+ cells"), invert = TRUE)

metadata <- Object_V1@meta.data


## Compute cell type specificity score, high resolution ------------------------------------------

Object_V1$major_celltypes <- Object_V1$cluster_labels_res.0.8 %>%
  lapply(FUN = switch,
         "Neutrophils_1" = "Neutrophils",
         "Neutrophils_2" = "Neutrophils",
         "Neutrophils_3" = "Neutrophils",
         "Neutrophils_4" = "Neutrophils",
         "Immature Neutrophils_1" = "Immature Neutrophils",
         "Immature Neutrophils_2" = "Immature Neutrophils",
         "Eosinophils" = "Eosinophils",
         "CD14_Monocytes_1" = "CD14 Monocytes",
         "CD14_Monocytes_2" = "CD14 Monocytes",
         "CD14_Monocytes_3" = "CD14 Monocytes",
         "CD16_Monocytes" = "CD16 Monocytes",
         "CD8_T_cells" = "CD8+ T cells",
         "CD4_T_cells_1" = "CD4+ T cells",
         "CD4_T_cells_2" = "CD4+ T cells",
         "CD4_T_cells_3" = "CD4+ T cells",
         "Prol. cells" = "Prol. cells",
         "NK_cells" = "NK cells",
         "B_cells_1" = "B cells",
         "B_cells_2" = "B cells",
         "Plasmablast" = "Plasmablasts",
         "Megakaryocytes" = "Megakaryocytes",
         "mDC" = "DC",
         "pDC" = "DC") %>% unlist() %>% factor()

Object_V1$major_celltypes <- Object_V1$major_celltypes %>%
  factor(levels = levels(.)[c(2, 3, 7, 10, 8, 11, 13, 4, 5, 6, 1, 12, 9)])

Object_V1@meta.data %>%
  .$major_celltypes %>%
  table()
#CD14 Monocytes      CD16 Monocytes          Eosinophils          Neutrophils Immature Neutrophils             NK cells          Prol. cells 
#          1897                329                  145                 7737                  202                 1105                   44 
#CD4+ T cells         CD8+ T cells                   DC              B cells         Plasmablasts       Megakaryocytes 
#        4409                 1209                  155                  997                   37                  399 

Object_V1 <- SetIdent(Object_V1, value = Object_V1$major_celltypes)


#compute average proportionate gene expression for each cell type
cluster <- Object_V1$major_celltypes %>%
  droplevels() %>%
  levels()
j <- 1
for(j in 1:length(cluster)){
  
  print(cluster[j])
  Object_cluster <- subset(x = Object_V1, idents = cluster[j])
  print(dim(Object_cluster))
  
  counts <- GetAssayData(object = Object_cluster, slot = "data", assay = "RNA")
  
  temp <- rowSums(counts) > 0
  counts_filtered <- counts[temp, ]
  print(dim(counts_filtered))
  
  prop_table <- data.frame(matrix(nrow = nrow(counts_filtered)))
  my_step <- 2500
  last_number <- ceiling(ncol(counts_filtered)/my_step)
  i <- 1
  for(i in 1:last_number){
    
    a <- my_step*(i-1)+1
    b <- ifelse(i < last_number, i*my_step, (ncol(counts_filtered) %% my_step) + (i-1)*my_step)
    
    prop_table_temp <- counts_filtered %>%
      .[, a:b] %>%
      as.matrix() %>%
      prop.table(margin = 2) %>% 
      as.data.frame()
    
    prop_table <- cbind(prop_table, prop_table_temp)
  }
  
  average_gene_per_group <- prop_table %>%
    .[, -1] %>%
    rowMeans(.) %>%
    data.frame("Gene" = names(.), "mean" = .)
  
  write.table(average_gene_per_group, paste0("results/cell_type_specificity/output/average_highRes_celltypes/", cluster[j], ".txt"), row.names = FALSE, sep = "\t", quote = FALSE)
}

#merge and clean tables
df_average_group <- data.frame()
j <- 1
for(j in 1:length(cluster)){
  
  table_temp <- read.table(paste0("results/cell_type_specificity/output/average_highRes_celltypes/", cluster[j], ".txt"), sep = "\t", header = TRUE)
  colnames(table_temp)[2] <- paste0("mean_", cluster[j])
  ifelse(j == 1, 
         df_average_group <- table_temp,
         df_average_group <- df_average_group %>% merge(table_temp, bx = "Gene", all = TRUE)
  )
}
rownames(df_average_group) <- df_average_group$Gene
df_average_group[, 2:ncol(df_average_group)] <- df_average_group[, 2:ncol(df_average_group)] %>%
  apply(function(x){ifelse(is.na(x), 0, x)}, MARGIN = 1:2)


#compute gene expression quantiles for tau score
lower_limit <- 10^{-7}
compute_bins <- function(x){
  y <- cut(x, breaks = c(0, lower_limit, quantile(x[x>lower_limit], 1:10*0.1)), include.lowest = TRUE)
  y <- factor(y, levels = levels(y), labels = 0:10)
  y <- as.numeric(y) - 1
  return(y)
}
df_average_group[, 2:ncol(df_average_group)] <- df_average_group[, 2:ncol(df_average_group)] %>%
  apply(compute_bins, MARGIN = 2) 

#Function to compute the gene specificity score from the gene expression profile
compute_tau <- function(x){
  max_component <- max(x)
  tau <- x %>% `/` (max_component) %>% {sum(1 - .)} %>% `/` (length(cluster) - 1)
  return(tau)
}

#compute cell type specificity (tau score)
cell_type_specificity <- df_average_group[, 2:ncol(df_average_group)] %>%
  apply(FUN = compute_tau, MARGIN = 1) %>%
  data.frame("Gene_ID" = names(.), "cell_type_specificity" = .)

#extract the cell types to which each gene is specific
df_cell_type_specificity <- cell_type_specificity %>%
  merge(df_average_group %>% add_column("Gene_ID" = rownames(.)), by = "Gene_ID") %>%
  add_column("max" = apply(.[, 4:ncol(.)], max, MARGIN = 1))

df_cell_name <- df_cell_type_specificity[, 4:(ncol(df_cell_type_specificity)-1)] %>% 
  apply(function(x){return(x == df_cell_type_specificity$max)}, MARGIN = 2) %>%
  as.data.frame() %>%
  add_column("Gene_ID" = df_cell_type_specificity$Gene_ID, .before = 1)
colnames(df_cell_name)[2:ncol(df_cell_name)] <- colnames(df_cell_name)[2:ncol(df_cell_name)] %>%
  str_sub(6, nchar(.))

for(i in 2:ncol(df_cell_name)){
  my_rowname <- colnames(df_cell_name)[i]
  df_cell_name[, my_rowname] <- ifelse(df_cell_name[, my_rowname], my_rowname, NA)
}

df_cell_name_2 <- df_cell_name[2:ncol(df_cell_name)] %>%
  apply(function(x){paste0(na.omit(as.vector(x)), collapse = ",")}, MARGIN = 1) 
df_cell_name$Cell_type <- df_cell_name_2

df_cell_type_specificity <- df_cell_type_specificity %>%
  merge(df_cell_name, by = "Gene_ID")

write.table(df_cell_type_specificity, "results/cell_type_specificity/output/tau_highRes_celltypes.txt", row.names = FALSE, quote = FALSE, sep = "\t")







## Test for enrichment of high ICC genes in cell type markers as defined by cell type specificity ---------------

#identify marker genes and merge with ICC analysis results
gene_attributes <- read.table("info_files/filtered_genes_samples25_fpm5.txt", sep = "\t", header = TRUE)
annotation <- gene_attributes[, c("Gene_ID", "gene_name")] %>% unique()

cell_type_specificity <- read.table("results/cell_type_specificity/output/tau_highRes_celltypes.txt", sep = "\t", header = TRUE) %>%
  merge(annotation, by.x = "Gene_ID", by.y = "gene_name")

df_marker <- cell_type_specificity %>%
  subset(.$Gene_ID.y %in% annotation$Gene_ID %>% unique()) %>%
  subset(str_count(.$Cell_type, ",") == 0) %>%
  subset(.$max > 1) %>%
  subset(.$gene_specificity > 0.5) %>%
  .[order(.$gene_specificity, decreasing = TRUE), ] %>%
  .[, c("Gene_ID.y", "Gene_ID", "gene_specificity", "max", "Cell_type")] %>%
  dplyr::rename("Gene_ID" = "Gene_ID.y", "Symbol" = "Gene_ID", "expression_quantile" = "max") %>%
  unique()

icc_results <- read.csv("results/ICC_analysis/output/gene_icc_data_whole_dataset.txt", sep = '\t')
icc_specificity <- inner_join(icc_results, gene_specificity, by = c("Gene_ID" = "Gene_ID.y"))
icc_specificity$Group <- ifelse(icc_specificity$ICC < 0.5, "Within", "Between")
icc_specificity$is_marker <- ifelse(icc_specificity$Gene_ID %in% df_marker$Gene_ID, "YES", "NO")




# Code for mosaic plot as in Fig. 2e ---------------

mosaic_plot <- ggplot(icc_specificity) + 
  geom_mosaic(aes(x = ggmosaic::product(Group, is_marker), fill = Group)) +
  geom_text(data = layer_data(last_plot(), 1) %>% filter(.wt > 0),
            aes(x = (xmin + xmax) / 2,
                y = (ymin + ymax) / 2,
                label = .wt), 
            size = 1.7) + 
  scale_fill_manual(values = c("Within"='#CCCC00', "Between"='#009900'), name="") + 
  xlab("") + ylab("") + 
  scale_x_productlist(labels = c("Not cell type\nspecific", "Cell type\nspecific")) + 
  scale_y_productlist(labels = c("High inter-\nindividual variance", "High intra-\nindividual variance")) +
  theme_void() +
  theme(axis.text.y = element_text(angle = 90), 
        legend.position = "none",
        aspect.ratio = 1,
        legend.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"))


chisq.test(icc_specificity$Group, icc_specificity$is_marker)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  icc_specificity$Group and icc_specificity$is_marker
#X-squared = 382.92, df = 1, p-value < 2.2e-16




# Code for correlation scatter plots as in Extended Data Fig. 4a ---------------

my_colors <- c("CD4+ T cells" = "#EA811F",
               "CD8+ T cells" = "maroon",
               "DC"= "peru",
               "NK cells" = "#6A3D8A",
               "B cells" = "goldenrod1",
               'Plasmablasts' = 'goldenrod4',
               'Prol. cells'='#8fdbaf',
               "Megakaryocytes" = "darkolivegreen",
               'Eosinophils'='thistle1',
               "Neutrophils" = "orchid3",
               "CD14 Monocytes" = "paleturquoise3",
               "CD16 Monocytes" = "#577676",
               "no marker" = "grey",
               "Immature Neutrophils" = "#99CCFF")

my_names <- c("CD4+ T cells" = "CD4+ T cells",
              "CD8+ T cells" = "CD8+ T cells",
              "DC"= "Dendritic cells",
              "NK cells" = "NK cells",
              "B cells" = "B cells",
              'Plasmablasts' = 'Plasmablasts',
              'Prol. cells'='Prol. cells',
              "Megakaryocytes" = "Megakaryocytes",
              'Eosinophils'='Eosinophils',
              "Neutrophils" = "Neutrophils",
              "Immature Neutrophils" = "Immature Neutrophils",
              "CD14 Monocytes" = "CD14 Monocytes",
              "CD16 Monocytes" = "CD16 Monocytes",
              "no marker" = "No marker")

celltype <- c("CD4+ T cells", "CD8+ T cells", "B cells", "Plasmablasts", "NK cells", "DC", "Eosinophils", "Neutrophils", "CD14 Monocytes", "CD16 Monocytes", "Megakaryocytes", "Immature Neutrophils", "Prol. cells")

effect_table <- data.frame()
i <- 1
for(i in 1:length(celltype)) {
  
  effect_table_temp <- icc_specificity %>%
    subset(str_count(.$Cell_type, ",") == 0) %>%
    subset(str_detect(.$Cell_type, fixed(celltype[i]))) %>%
    add_column("my_group" = celltype[i]) %>%
    add_column("my_color" = ifelse(.$is_marker == "YES", celltype[i], "no marker"))
  
  effect_table <- rbind(effect_table, effect_table_temp)
}

effect_table$my_group <- factor(effect_table$my_group, levels = celltype)

cor_plot <- ggplot(effect_table, aes(x = cell_type_specificity, y = ICC)) + 
  geom_point(size = my_point_size/4, aes(col = my_color)) +
  scale_color_manual(values = my_colors, labels = as_labeller(my_names), breaks = c(celltype, "no marker")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"), limits = c(0, 1)) +
  geom_smooth(method = "glm", formula = y ~ x, col = "black") +
  stat_cor(method = "spearman", label.y = 0.95, size = 1.7) + 
  facet_wrap(. ~ my_group, nrow = 2, labeller = as_labeller(my_names)) +
  labs(x = "Gene specificity", y = "ICC", color = "Marker for") + 
  theme(legend.position = "right",
        strip.background = element_rect(color = "white", fill = "white"),
        panel.border = element_rect(color = "black"),
        panel.grid = element_blank())







# Code for dotplot as in Extended Data Fig. 4b ---------------

#select genes identified as cell type markers and with high ICC
my_features <- c("ITGA2B", "TUBB1", "IGHA2", "IGHG4", "CD72", "FCRLA", "CD19", "BLK", "CLEC4C", "LILRA4", "ZNF683", "CXCR6", "NELL2", "CCR4", "TCF19", "CD160", "NCR1", "LCN2", "MMP8", "IL5RA", "C1QA", "CES1")

dot_plot_marker_genes <- DotPlot(Object_V1, features = my_features,  dot.scale = 4)  + 
  RotatedAxis() +
  scale_color_gradientn(colours  = rev(brewer.pal(n = 9, name = "RdBu")), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  guides(color = guide_colourbar(title.position = "top", title = "Average\nexpression"),
         size = guide_legend(title.position = "top", ncol = 1, title = "Percent\nexpressed")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "right",
        panel.grid = element_blank()) +
  coord_flip()

df_boxplots <- df_marker %>%
  merge(icc_results, by = "Gene_ID") %>%
  subset(.$Symbol %in% my_features) %>%
  .[match(my_features, .$Symbol), ] %>%
  .[, c("Symbol", "ICC", "cell_type_specificity")] %>%
  melt()
df_boxplots$Symbol <- factor(df_boxplots$Symbol, levels = my_features)

boxplot_combined <- ggplot(df_boxplots, aes(x = Symbol, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge", color = "black", width = 0.7) +
  labs(y = "ICC") +
  scale_fill_manual(values = c("ICC" = '#009900', "cell_type_specificity" = "#044B1D"), labels = c("ICC", "Cell type\nspecificity")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "right",
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

dot_plot_marker_genes + boxplot_combined + plot_layout(ncol = 2, widths = c(0.8, 0.2))






## Code for GO dot plot as in Extended Data Fig. 4c ---------------

#GO enrichment analysis
mean_expression <- read.table("results/mean_gene_expression.txt", header = TRUE)
rownames(mean_expression) <- mean_expression$Gene_ID

celltype <- c("CD4+ T cells", "CD8+ T cells", "B cells", "Plasmablasts", "NK cells", "DC", "Eosinophils", "Neutrophils", "CD14 Monocytes", "CD16 Monocytes", "Megakaryocytes", "Immature Neutrophils", "Prol. cells")
direction <- c("between", "within")
for(c in 1:length(celltype)){
  
  icc_specificity_temp <- icc_specificity %>%
    subset(str_count(.$Cell_type, ",") == 0) %>%
    subset(str_detect(.$Cell_type, fixed(celltype[c]))) %>%
    subset(.$is_marker == "YES")
  
  for(d in 1:length(direction)){
    
    if(direction[d] == "between"){
      degs <- icc_specificity_temp %>%
        subset(.$is_marker == "YES" & .$Group == "Between") %>%
        .$Gene_ID
    }else if(direction[d] == "within"){
      degs <- icc_specificity_temp %>%
        subset(.$is_marker == "YES" & .$Group == "Within") %>%
        .$Gene_ID
    }
    
    base_mean <- mean_expression %>%
      subset(rownames(.) %in% icc_specificity_temp)
    overallBaseMean <- as.matrix(base_mean[, "mean_expression", drop = F])
    
    topGOResults <- run_GO_ORA(degs, 
                               overallBaseMean, 
                               selected_database = "org.Hs.eg.db", 
                               back_num = 10,
                               same_expression_level = FALSE,
                               annotation = annotation)
    
    go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
    
    file_filtered <- paste0("results/cell_type_specificity/output/ORA_GO_cellType_specific_terms/ORA_GO_", celltype[c], "_YES_", direction[d], ".csv")
    write.csv(go_results_filtered, file_filtered, row.names = FALSE)
  }
}

#plot GO results
celltype <- c("CD4+ T cells", "CD8+ T cells", "B cells", "NK cells", "Neutrophils")
direction <- c("between", "within")
GO_tables <- list()

for(c in 1:length(celltype)){
  for(d in 1:length(direction)){
    
    my_error <- try({read.csv(paste0("results/cell_type_specificity/output/ORA_GO_cellType_specific_terms/ORA_GO_", celltype[c], "_YES_", direction[d], ".csv"), 
                              sep = ",", 
                              header = TRUE)}, silent = TRUE)
    
    if(is(my_error, 'try-error')){
      GO_tables <- append(GO_tables,
                          list(NULL))
    }else{
      table_temp <- read.csv(paste0("results/cell_type_specificity/output/ORA_GO_cellType_specific_terms/ORA_GO_", celltype[c], "_YES_", direction[d], ".csv"), 
                             sep = ",", 
                             header = TRUE)
      GO_tables <- append(GO_tables,
                          list(table_temp))
    }
  }
}

names(GO_tables) <- paste(rep(celltype, each = 2), rep(direction, length(celltype)), sep = "_")

GO_tables <- GO_tables %>%
  Filter(Negate(is.null), .)

GO_tables[["B cells_between"]] <- GO_tables[["B cells_between"]][c(1, 4, 5), ]
GO_tables[["CD8+ T cells_between"]] <- GO_tables[["CD8+ T cells_between"]][c(1, 3, 4), ]
GO_tables[["CD4+ T cells_between"]] <- GO_tables[["CD4+ T cells_between"]][c(1, 4, 5), ]

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = names(GO_tables),
                     score_name = score_name,
                     onto = "ALL",
                     reverse_bool = TRUE,
                     showTerms = 3,
                     multLines = TRUE,
                     numChar = 80)

df_plot$Group1 <- df_plot$Direction %>% str_split("_") %>% do.call(rbind, .) %>% .[, 1]
df_plot$Group1 <- factor(df_plot$Group1, levels = df_plot$Group1 %>% unique())
df_plot$Group2 <- df_plot$Direction %>% str_split("_") %>% do.call(rbind, .) %>% .[, 2]
df_plot$Group2 <- ifelse(grepl("between", df_plot$Group2, fixed = TRUE), "High inter-\nindividual variance", "High intra-\nindividual variance")

GO_dotplot <- ggplot(df_plot, mapping = aes(x = Group2, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  facet_grid(~ Group1, switch = 'x', scales = "free_x") + 
  scale_colour_gradient(high="#990000", low="#FF9999", name = expression(-~log["10"]~italic(p))) + 
  scale_size(name = "Gene ratio", range = c(0, 3)) +
  guides(color = guide_colourbar(title.position = "top"),
         size = guide_legend(title.position = "top", ncol = 2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title = element_blank(),
        strip.placement = "outside", 
        strip.background = element_blank(), 
        panel.spacing = unit(0, "lines"))







