## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggmosaic)
library(tidyverse)
library(arrow)
library(data.table)
library(readr)
library(rstatix)

setwd("C:/Users/n.mishra/Temporal_gene_expression_variation")


## Calculating eQTL discovery rates  ------------------------------------------------------------

gtex_icc_eGene <- read.csv("results/eQTL_based_analysis/Whole_dataset/GTEx_v10/ICC_cis_eGene_data.txt", sep = '\t')

gtex_cis_discovery_rates <- gtex_icc_eGene %>%
  group_by(ICC_bin) %>%
  summarise(
    n_genes  = n(),
    n_eGenes = sum(eGene_binary),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    eQTL_rate      = n_eGenes / n_genes,
    eQTL_rate_low  = binom.test(n_eGenes, n_genes)$conf.int[1],
    eQTL_rate_high = binom.test(n_eGenes, n_genes)$conf.int[2], 
    dataset       = "GTEx v10 (cis)"
  ) %>%
  ungroup() %>%
  filter(!is.na(ICC_bin))

eqtl_cis_icc_eGene <- read.csv("results/eQTL_based_analysis/Whole_dataset/eQTLGen_cis/ICC_cis_eGene_data.txt", sep = '\t')

eqtl_cis_discovery_rates <- eqtl_cis_icc_eGene %>%
  group_by(ICC_bin) %>%
  summarise(
    n_genes  = n(),
    n_eGenes = sum(eGene_binary),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    eQTL_rate      = n_eGenes / n_genes,
    eQTL_rate_low  = binom.test(n_eGenes, n_genes)$conf.int[1],
    eQTL_rate_high = binom.test(n_eGenes, n_genes)$conf.int[2], 
    dataset       = "eQTLGen Phase I (cis)"
  ) %>%
  ungroup() %>%
  filter(!is.na(ICC_bin))

eqtl_trans_icc_eGene <- read.csv("results/eQTL_based_analysis/Whole_dataset/eQTLGen_trans/ICC_trans_eGene_data.txt", sep = '\t')

eqtl_trans_discovery_rates <- eqtl_trans_icc_eGene %>%
  group_by(ICC_bin) %>%
  summarise(
    n_genes  = n(),
    n_eGenes = sum(eGene_binary),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    eQTL_rate      = n_eGenes / n_genes,
    eQTL_rate_low  = binom.test(n_eGenes, n_genes)$conf.int[1],
    eQTL_rate_high = binom.test(n_eGenes, n_genes)$conf.int[2], 
    dataset       = "eQTLGen Phase I (trans)"
  ) %>%
  ungroup() %>%
  filter(!is.na(ICC_bin))

discovery_rates <- rbind(gtex_cis_discovery_rates, eqtl_cis_discovery_rates, eqtl_trans_discovery_rates)
discovery_rates$dataset <- factor(discovery_rates$dataset, levels = c("GTEx v10 (cis)", "eQTLGen Phase I (cis)", "eQTLGen Phase I (trans)"))

write.table(discovery_rates, "results/eQTL_based_analysis/Whole_dataset/eQTL_discovery_rates.txt", sep = '\t', row.names = FALSE, quote = FALSE)

dodge <- position_dodge(width = 0.9)

## ICC bin discovery rate plot as in Fig. 6d ---------------


p <- ggplot(discovery_rates,
            aes(x = dataset, y = eQTL_rate, col = ICC_bin, fill = ICC_bin)) +
  geom_col(alpha = 0.7, position = dodge) +
  geom_errorbar(aes(ymin = eQTL_rate_low,
                    ymax = eQTL_rate_high),
                width = 0.2, 
                position = dodge) +
  ylab("eQTL discovery rate (proportion of eGenes)") +
  xlab("") +
  scale_color_manual(values = c("#CCCC00", "#99CC00", "#66B200", "#009900")) + 
  scale_fill_manual(values = c("#CCCC00", "#99CC00", "#66B200", "#009900")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14)
  )
p

## Effect size comparison -----------------------------------------------

df1 <- gtex_icc_eGene %>% 
  dplyr::select(Gene_ID, ICC, ICC_bin, Group, eGene, slope, n_eQTL) %>%
  rename_with(~ paste0(.x, "_GTEx_cis"), -c(Gene_ID, ICC, ICC_bin, Group))
df2 <- eqtl_cis_icc_eGene %>%
  dplyr::select(Gene_ID, eGene, slope, n_eQTL) %>%
  rename_with(~ paste0(.x, "_eQTLGen_cis"), -Gene_ID)
df3 <- eqtl_trans_icc_eGene %>%
  dplyr::select(Gene_ID, eGene, slope, n_eQTL) %>%
  rename_with(~ paste0(.x, "_eQTLGen_trans"), -Gene_ID)

icc_eGene <- purrr::reduce(
  list(df1, df2, df3),
  left_join,
  by = "Gene_ID"
)

icc_slope <- icc_eGene %>%
  dplyr::select(Gene_ID, ICC, eGene_GTEx_cis, eGene_eQTLGen_cis, eGene_eQTLGen_trans, slope_GTEx_cis, slope_eQTLGen_cis, slope_eQTLGen_trans) %>%
  pivot_longer(
    cols = c(slope_GTEx_cis, slope_eQTLGen_cis, slope_eQTLGen_trans),
    names_to = "dataset",
    values_to = "slope"
  ) %>%
  mutate(
    eGene = case_when(
      dataset == "slope_GTEx_cis"       ~ eGene_GTEx_cis,
      dataset == "slope_eQTLGen_cis"    ~ eGene_eQTLGen_cis,
      dataset == "slope_eQTLGen_trans"  ~ eGene_eQTLGen_trans,
      TRUE ~ NA_character_
    ),
    dataset = recode(dataset,
                     slope_GTEx_cis      = "GTEx v10 (cis)",
                     slope_eQTLGen_cis   = "eQTLGen Phase I (cis)",
                     slope_eQTLGen_trans = "eQTLGen Phase I (trans)"
    )
  ) %>%
  filter(eGene == "Yes", !is.na(slope))

icc_slope$dataset <- factor(icc_slope$dataset, levels = c("GTEx v10 (cis)", "eQTLGen Phase I (cis)", "eQTLGen Phase I (trans)"))

## Scatter plot for ICC vs. slope as in Extended Data Fig. 9b ----------------------

icc_slope_plot <- ggplot(icc_slope, aes(x = ICC, y = abs(slope))) +
  geom_point(size = my_point_size/2, col = "#D1BEB0", alpha = 0.5) +
  geom_smooth(method = "glm", formula = y ~ x, col = "#383F51") +
  stat_cor(method = "spearman", cor.coef.name = "rho") +
  facet_wrap(~dataset, scales = "free_y") +
  xlab("ICC") + ylab("Absolute eQTL slope/z-score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14)
  )

icc_slope_plot

## Scatter plot for ICC vs. number of eQTLs as in Extended Data Fig. 9c ----------------------

icc_nqtl <- icc_eGene %>%
  dplyr::select(Gene_ID, ICC, eGene_GTEx_cis, eGene_eQTLGen_cis, eGene_eQTLGen_trans, n_eQTL_GTEx_cis, n_eQTL_eQTLGen_cis, n_eQTL_eQTLGen_trans) %>%
  pivot_longer(
    cols = c(n_eQTL_GTEx_cis, n_eQTL_eQTLGen_cis, n_eQTL_eQTLGen_trans),
    names_to = "dataset",
    values_to = "n_eQTL"
  ) %>%
  mutate(
    eGene = case_when(
      dataset == "n_eQTL_GTEx_cis"       ~ eGene_GTEx_cis,
      dataset == "n_eQTL_eQTLGen_cis"    ~ eGene_eQTLGen_cis,
      dataset == "n_eQTL_eQTLGen_trans"  ~ eGene_eQTLGen_trans,
      TRUE ~ NA_character_
    ),
    dataset = recode(dataset,
                     n_eQTL_GTEx_cis      = "GTEx v10 (cis)",
                     n_eQTL_eQTLGen_cis   = "eQTLGen Phase I (cis)",
                     n_eQTL_eQTLGen_trans = "eQTLGen Phase I (trans)"
    )
  ) %>%
  filter(eGene == "Yes", !is.na(n_eQTL))

icc_nqtl$dataset <- factor(icc_nqtl$dataset, levels = c("GTEx v10 (cis)", "eQTLGen Phase I (cis)", "eQTLGen Phase I (trans)"))


icc_nqtl_plot <- ggplot(icc_nqtl, aes(x = ICC, y = n_eQTL)) +
  geom_point(size = my_point_size/2, col = "#A49966", alpha = 0.5) +
  geom_smooth(method = "glm", formula = y ~ x, col = "#363020") +
  stat_cor(method = "spearman", cor.coef.name = "rho") +
  facet_wrap(~dataset, scales = "free_y") +
  xlab("ICC") + ylab("Number of eQTLs") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14)
  )

icc_nqtl_plot


## Enrichment of eGenes in low and high ICC groups and plotting mosaic plots as in Fig. 6a-c-----------------------

icc_eqtl_enrichment <- function(dataset_name,
                                  icc_eGene,
                                  out_base   = "results/eQTL_based_analysis/Whole_dataset") {
  out_dir <- file.path(out_base, dataset_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  icc_eGene$Group <- ifelse(icc_eGene$ICC < 0.5, "Within", "Between")
  
  # Chi-square + mosaic
  chisq_result <- chisq_test(table(icc_eGene$Group, icc_eGene$eGene))
  
  plot_data <- icc_eGene[, c("eGene", "Group")] %>%
    mutate(across(c(eGene, Group), as.factor))
  
  pdf(file.path(out_dir, "cis_eGene_variance_mosaicplot.pdf"))
  p <- ggplot(plot_data) +
    geom_mosaic(aes(x = product(Group, eGene), fill = Group)) +
    geom_text(
      data = layer_data(last_plot(), 1) %>% dplyr::filter(.wt > 0),
      aes(x = (xmin + xmax) / 2,
          y = (ymin + ymax) / 2,
          label = .wt)
    ) +
    scale_fill_manual(
      values = c("Within" = '#CCCC00', "Between" = '#009900'),
      name = ""
    ) +
    xlab("eGene") + ylab("High variance") +
    theme_void() +
    theme(
      axis.text  = element_text(size = 15),
      axis.title = element_text(size = 15),
      axis.title.y = element_text(angle = 90),
      legend.text  = element_text(size = 15),
      legend.position = "top",
      aspect.ratio = 1,
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  print(p)
  dev.off()
  
  write.table(chisq_result,
              file.path(out_dir, "cis_eGene_variance_chisq_test.txt"),
              sep = '\t', quote = FALSE)
  write.table(plot_data,
              file.path(out_dir, "cis_eGene_variance_mosaicplot_data.txt"),
              sep = '\t', quote = FALSE, row.names = FALSE)

}

icc_eqtl_enrichment("GTEx_cis", gtex_icc_eGene)
icc_eqtl_enrichment("eQTLGen_cis", eqtl_cis_icc_eGene)
icc_eqtl_enrichment("eQTLGen_trans", eqtl_trans_icc_eGene)


## Calculating eQTL replication rates  ------------------------------------------------------------

calc_replication_by_ICC <- function(eGene_discovery,   
                                    eGene_replication, 
                                    dataset_label      
) {
  # Define ICC bins
  breaks <- c(0, 0.25, 0.5, 0.75, 1)
  labels <- c("ICC_0-0.25", "ICC_0.25-0.5", "ICC_0.5-0.75", "ICC_0.75-1")
  
  # Universe = all genes in ICC table
  all_genes <- unique(icc_results$Gene_ID)
  
  # Prepare eGene tables with explicit No for missing genes
  disc_tbl <- data.frame(Gene_ID = all_genes) %>%
    left_join(eGene_discovery %>% dplyr::select(Gene_ID, eGene, ICC), by = "Gene_ID") %>%
    mutate(eGene_disc = ifelse(is.na(eGene), "No", eGene)) %>%
    dplyr::select(Gene_ID, eGene_disc)
  
  repl_tbl <- data.frame(Gene_ID = all_genes) %>%
    left_join(eGene_replication %>% dplyr::select(Gene_ID, eGene), by = "Gene_ID") %>%
    mutate(eGene_rep = ifelse(is.na(eGene), "No", eGene)) %>%
    dplyr::select(Gene_ID, eGene_rep)
  
  # Merge with ICC and define ICC bins
  rep_df <- disc_tbl %>%
    left_join(repl_tbl, by = "Gene_ID") %>%
    mutate(
      ICC_bin = cut(ICC,
                    breaks = breaks,
                    labels = labels,
                    include.lowest = TRUE,
                    right = TRUE)
    )
  
  # Keep only discovery eGenes for replication analysis
  rep_df_disc <- rep_df %>%
    filter(eGene_disc == "Yes", !is.na(ICC_bin))
  
  return(rep_df_disc)
}


# GTEx (discovery) → eQTLGen cis (replication)
gtex_vs_eqtlgen_cis <- calc_replication_by_ICC(
  eGene_discovery  = gtex_icc_eGene,          
  eGene_replication= eqtl_cis_icc_eGene,      
    dataset_label    = "GTEx cis → eQTLGen cis"
)


write.table(gtex_vs_eqtlgen_cis, "results/eQTL_based_analysis/Whole_dataset/GTEx_v10_vs_eQTLGen_cis/replication_GTEx_to_eQTLGen.txt", 
            quote = FALSE, row.names = FALSE, sep = '\t')

eqtlgen_vs_gtex_cis <- calc_replication_by_ICC(
  eGene_discovery  = eqtl_cis_icc_eGene,  
  eGene_replication= gtex_icc_eGene,      
  dataset_label    = "eQTLGen cis → GTEx cis"
)

write.table(eqtlgen_vs_gtex_cis, "results/eQTL_based_analysis/Whole_dataset/eQTLGen_cis_vs_GTEx_v10/replication_eQTLGen_to_GTEx.txt", 
            quote = FALSE, row.names = FALSE, sep = '\t')

gtex_to_eqtlgen_replication_rates <- gtex_vs_eqtlgen_cis %>%
  filter(eGene_disc == "Yes", !is.na(ICC_bin)) %>%
  group_by(ICC_bin) %>%
  summarise(
    n_discovery = n(),
    n_rep       = sum(eGene_rep == "Yes"),
    repl_rate   = n_rep / n_discovery,
    repl_rate_low = binom.test(n_rep, n_discovery)$conf.int[1],
    repl_rate_high = binom.test(n_rep, n_discovery)$conf.int[2],
    .groups     = "drop"
  ) %>%
  mutate(
    dataset   = "GTEx cis → eQTLGen cis"
  )

eqtlgen_to_gtex_replication_rates <- eqtlgen_vs_gtex_cis %>%
  filter(eGene_disc == "Yes", !is.na(ICC_bin)) %>%
  group_by(ICC_bin) %>%
  summarise(
    n_discovery = n(),
    n_rep       = sum(eGene_rep == "Yes"),
    repl_rate   = n_rep / n_discovery,
    repl_rate_low = binom.test(n_rep, n_discovery)$conf.int[1],
    repl_rate_high = binom.test(n_rep, n_discovery)$conf.int[2],
    .groups     = "drop"
  ) %>%
  mutate(
    dataset   = "eQTLGen cis → GTEx cis"
  )

replication_rates <- rbind(gtex_to_eqtlgen_replication_rates, eqtlgen_to_gtex_replication_rates)
replication_rates$dataset <- factor(replication_rates$dataset, levels = c("GTEx cis → eQTLGen cis", "eQTLGen cis → GTEx cis"))

dodge <- position_dodge(width = 0.9)

## ICC bin replication rate plot as in Extended Data Fig. 9a ---------------


replication_rates_plot <- ggplot(replication_rates,
                                 aes(x = dataset, y = repl_rate, col = ICC_bin, fill = ICC_bin)) +
  geom_col(alpha = 0.7, position = dodge) +
  geom_errorbar(aes(ymin = repl_rate_low,
                    ymax = repl_rate_high),
                width = 0.2,
                position = dodge) +
  ylab("eQTL replication rate\n(proportion of eGenes)") +
  xlab("ICC bin (discovery eGenes)") +
  scale_color_manual(values = c("#CCCC00", "#99CC00", "#66B200", "#009900"), guide = "none") + 
  scale_fill_manual(values = c("#CCCC00", "#99CC00", "#66B200", "#009900"), 
                    labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"), 
                    name = "ICC") +
  guides(fill  = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14)
  )

replication_rates_plot
