## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape)
library(Seurat)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")


## Loading data -----------------------------------

#scObject from Schulte-Schrepping et al. - Severe COVID-19 Is Marked by a Dysregulated Myeloid Cell Compartment (https://doi.org:10.1016/j.cell.2020.08.001)
single_cell_object <- readRDS("info_files/seurat_COVID19_freshWB-PBMC_cohort2_rhapsody_jonas_FG_2020-08-18.rds")
single_cell_object <- SetIdent(single_cell_object, value = single_cell_object@meta.data$cluster_labels_res.0.8)

#keep only the samples from healthy individuals with no comorbidities
Object_V1 <- subset(x = single_cell_object, 
                    cells = colnames(single_cell_object[, single_cell_object$comorbidities == "none" & single_cell_object$diagnosis == "control"]))
Object_V1 <- subset(x = Object_V1, idents = c("Mixed_cells", "CD34+ GATA2+ cells"), invert = TRUE)

metadata <- Object_V1@meta.data



## Compute cell type specificity score, low resolution ------------------------------------------

#merge cell types to match measured blood cell count types
Object_V1$major_celltypes <- Object_V1$cluster_labels_res.0.8 %>%
  lapply(FUN = switch,
         "Neutrophils_1" = "Neutrophils",
         "Neutrophils_2" = "Neutrophils",
         "Neutrophils_3" = "Neutrophils",
         "Neutrophils_4" = "Neutrophils",
         "Immature Neutrophils_1" = "Neutrophils",
         "Immature Neutrophils_2" = "Neutrophils",
         "Eosinophils" = "Eosinophils",
         "CD14_Monocytes_1" = "Monocytes",
         "CD14_Monocytes_2" = "Monocytes",
         "CD14_Monocytes_3" = "Monocytes",
         "CD16_Monocytes" = "Monocytes",
         "CD8_T_cells" = "Lymphocytes",
         "CD4_T_cells_1" = "Lymphocytes",
         "CD4_T_cells_2" = "Lymphocytes",
         "CD4_T_cells_3" = "Lymphocytes",
         "Prol. cells" = "Prol. cells",
         "NK_cells" = "Lymphocytes",
         "B_cells_1" = "Lymphocytes",
         "B_cells_2" = "Lymphocytes",
         "Plasmablast" = "Lymphocytes",
         "Megakaryocytes" = "Megakaryocytes",
         "mDC" = "DC",
         "pDC" = "DC") %>% unlist() %>% factor()

Object_V1@meta.data %>%
  .$major_celltypes %>%
  table()
# DC    Eosinophils    Lymphocytes Megakaryocytes      Monocytes    Neutrophils    Prol. cells 
#155            145           7757            399           2226           7939             44 

Object_V1 <- SetIdent(Object_V1, value = Object_V1$major_celltypes)


#compute average proportionate gene expression for each cell type
cluster <- Object_V1$major_celltypes %>%
  droplevels() %>%
  levels()
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
  
  write.table(average_gene_per_group, paste0("results/cell_type_specificity/output/average_lowRes_celltypes/", cluster[j], ".txt"), row.names = FALSE, sep = "\t", quote = FALSE)
}

#merge and clean tables
df_average_group <- data.frame()
j <- 1
for(j in 1:length(cluster)){
  
  table_temp <- read.table(paste0("results/cell_type_specificity/output/average_lowRes_celltypes/", cluster[j], ".txt"), sep = "\t", header = TRUE)
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

#Function to compute the cell type specificity score from the gene expression profile
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

df_cell_name_merged <- df_cell_name[2:ncol(df_cell_name)] %>%
  apply(function(x){paste0(na.omit(as.vector(x)), collapse = ",")}, MARGIN = 1) 
df_cell_name$Cell_type <- df_cell_name_merged

df_cell_type_specificity <- df_cell_type_specificity %>%
  merge(df_cell_name, by = "Gene_ID")

write.table(df_cell_type_specificity, "results/cell_type_specificity/output/tau_lowRes_celltypes.txt", row.names = FALSE, quote = FALSE, sep = "\t")




## Test for enrichment of high variance genes in cell type markers as defined by cell type specificity ---------------

tx2gene <- read.table("info_files/tx2gene.txt")
annotation <- tx2gene[, c("gene_id", "gene_name")] %>% 
  unique()

df_cell_type_specificity <- read.table("results/cell_type_specificity/output/tau_lowRes_celltypes.txt", sep = "\t", header = TRUE) %>%
  merge(annotation, by.x = "Gene_ID", by.y = "gene_name")

res_varPart <- read.table("results/Variance_partition/output/variance_partition.txt", header = TRUE)
res_varPart <- res_varPart %>%
  add_column("gene_id" = rownames(.))

my_colors <- c("#E69F00", "#0072B2")

all_genes <- res_varPart$gene_id %>% unique()
celltype <- c("Neutrophils...µL.", "Monocytes...µL.", "Eosinophils...µL.", "Lymphocytes...µL.")
celltype_marker <- c("Neutrophils", "Monocytes", "Eosinophils", "Lymphocytes")


#Fishers test to test if the genes which expression is significantly affected by cell type composition 
#are enriched with marker genes for the corresponding cell types
gene_table <- data.frame("gene_id" = all_genes) %>%
  merge(annotation, by = "gene_id")
df_enrichment <- data.frame()
df_enrichment_pval <- data.frame()

for(i in 1:length(celltype)) {
  
  highVar_genes <- res_varPart %>%
    .[order(.[, celltype[i]], decreasing = TRUE), ] %>%
    head(1000) %>%
    .$gene_id %>%
    unique()
  not_highVar_genes <- setdiff(all_genes, highVar_genes) %>% unique()
  
  marker_for_celltype <- df_cell_type_specificity %>%
    subset(.$Cell_type == celltype_marker[i] & .$gene_id %in% all_genes) %>%
    subset(.$max > 1) %>%
    subset(.$cell_type_specificity > 0.5) %>%
    .[order(.$cell_type_specificity, decreasing = TRUE), ] %>%
    .$gene_id %>% 
    unique()
  
  no_marker <- setdiff(all_genes, marker_for_celltype) %>% unique()
  
  highVar_marker <- intersect(marker_for_celltype, highVar_genes) %>% unique()
  not_highVar_marker <- intersect(marker_for_celltype, not_highVar_genes) %>% unique()
  highVar_no_marker <- intersect(highVar_genes, no_marker) %>% unique()
  not_highVar_no_marker <- intersect(not_highVar_genes, no_marker)
  
  
  fisher_matrix <- matrix(c(length(highVar_marker), length(not_highVar_marker), 
                            length(highVar_no_marker), length(not_highVar_no_marker)),
                          nrow = 2,
                          dimnames = list(c("highVar", "not highVar"),
                                          c("marker", "no marker")))
  
  res_fisher <- fisher.test(fisher_matrix, alternative = "two.sided")
  
  #compute the expected number of significant markers assuming independence
  #P(sig + marker) = P(sig)*P(marker)=
  prob <- (length(highVar_genes)/length(all_genes)) * (length(marker_for_celltype)/length(all_genes))
  expected_number_highVar_marker <- prob*length(all_genes)
  
  df_enrichment_pval <- rbind(df_enrichment_pval, 
                              data.frame("Cell_type" = celltype_marker[i],
                                         "p_value" = res_fisher$p.value,
                                         "Odds_ratio" = res_fisher$estimate,
                                         "conf_int_low" = res_fisher$conf.int[1],
                                         "conf_int_high" = res_fisher$conf.int[2],
                                         "highVar_marker" = length(highVar_marker),
                                         "expected_highVar_marker" = expected_number_highVar_marker,
                                         "total_number_of_markers" = length(marker_for_celltype)))
  
  df_enrichment_temp <- gene_table %>%
    add_column("Cell_type" = celltype_marker[i],
               "Variance" = ifelse(.$gene_id %in% highVar_genes, "high variance", "low variance"),
               "Marker" = ifelse(.$gene_id %in% marker_for_celltype, "marker", "no marker"))
  
  df_enrichment <- rbind(df_enrichment, df_enrichment_temp)
  
}

df_enrichment_pval$adj_p <- p.adjust(df_enrichment_pval$p_value, method = "BH")




## Code for barplot as in Extended Data Fig. 2e ---------------

df_enrichment$Cell_type <- factor(df_enrichment$Cell_type, levels = c("Neutrophils", "Monocytes", "Lymphocytes", "Eosinophils"))
df_enrichment$Variance <- factor(df_enrichment$Variance, levels = c("low variance", "high variance"))
df_enrichment$Marker <- factor(df_enrichment$Marker, levels = c("no marker", "marker"))

df_enrichment_pval <- df_enrichment_pval %>%
  add_significance(p.col = "adj_p", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  add_column("Marker" = "no marker",
             "Variance" = "low variance",
             "x_pos" = "high variance", 
             "y_pos" = 0.85,
             "my_label" = paste0("OR=", round(.$Odds_ratio, 2), "  \n(", .$p.adj_signif, ")"))

df_enrichment_pval$Cell_type <- factor(df_enrichment_pval$Cell_type, levels = c("Neutrophils", "Monocytes", "Lymphocytes", "Eosinophils"))

my_colors <- c("#E69F00", "#0072B2")

enrichment_plot <- ggplot(df_enrichment) +
  geom_bar(aes(x = Marker, fill = Variance), alpha = 1, position = "fill", linewidth = my_line_width) +
  scale_fill_manual(breaks = c("low variance", "high variance"), values = my_colors, labels = c("low", "high")) +
  scale_x_discrete(breaks = c("no marker", "marker"), labels = c("low", "high")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100)) +
  labs(fill = "Variance explained", x = "Cell type specificity", y = "Percentage") +
  geom_text(data = df_enrichment_pval, aes(x = Marker, y = y_pos, label = my_label), color = "black", size = my_axis_text_size/.pt) + 
  facet_wrap(~ Cell_type) +
  my_theme_ggplot +
  theme(legend.position = "top",
        strip.placement = "outside",
        axis.title.y = element_text(margin = margin(0,0,0,0.15, unit = "cm")),
        strip.background = element_blank(),
        panel.grid = element_blank())



         
