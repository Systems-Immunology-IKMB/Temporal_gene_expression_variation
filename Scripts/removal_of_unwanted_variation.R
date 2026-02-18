## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(DESeq2)
library(RUVSeq)
library(variancePartition)
library(reshape2)
library(patchwork)

setwd("C:/Users/f.kimmig/Temporal_gene_expression_variation")


## Loading data ------------------------------------------------------------

col_clinical_data <- read.table("info_files/SYSCID_metadata_merged_2023_07_21.txt", sep = "\t", header = TRUE)
col_clinical_data$Individual_ID <- factor(col_clinical_data$Individual_ID)
rownames(col_clinical_data) <- col_clinical_data$Sample_ID
col_clinical_data$Run <- factor(col_clinical_data$Run)

col_clinical_data$Annual_season <- ifelse(col_clinical_data$Month %in% c("January", "February", "March"), "Winter", 
                                 ifelse(col_clinical_data$Month %in% c("July", "August", "September"), "Summer",
                                        ifelse(col_clinical_data$Month %in% c("April", "May"), "Spring", "Autumn")))
col_clinical_data$Annual_season <- factor(col_clinical_data$Annual_season, levels = c("Winter", "Spring", "Summer", "Autumn"))

#import count data
count_data <- read.table("count_files/merged_gene_counts_all_samples.txt", sep = '\t')

colnames(count_data) <- substr(colnames(count_data), 1, 6)
col_clinical_data <- subset(col_clinical_data, col_clinical_data$Sample_ID %in% colnames(count_data))
count_data <- count_data[, as.character(col_clinical_data$Sample_ID)]
colnames(count_data) <- col_clinical_data$Sample_ID

#exclude genes on chromosome X and Y from this analysis
gene_attributes <- read.table("info_files/filtered_genes_samples25_fpm5.txt", sep = "\t", header = TRUE)
sex_chromosome_genes <- subset(gene_attributes, gene_attributes$Chr %in% c("chrX", "chrY"))$Gene_ID
gene_attributes <- gene_attributes %>%
  subset(!.$Gene_ID %in% sex_chromosome_genes)

count_data <- count_data[as.character(gene_attributes$Gene_ID), ]
rownames(count_data) <- gene_attributes$Gene_ID

dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_clinical_data, design = ~ 1)




## Adapt RUVg function to work with selected factors --------------------
#see also 'https://github.com/drisso/RUVSeq/blob/master/R/RUVg-methods.R'

#This functions takes a count matrix x, a list of negative control genes cIdx, the number k of total factors to consider and the selected factors that should be used
RUVg_factor <- function(x, cIdx, k, factors, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8, isLog=FALSE) {
  
  if(!isLog && !all(.isWholeNumber(x))) {
    warning(paste0("The expression matrix does not contain counts.\n",
                   "Please, pass a matrix of counts (not logged) or set isLog to TRUE to skip the log transformation"))
  }
  
  if(isLog) {
    Y <- t(x)
  } else {
    Y <- t(log(x+epsilon))
  }
  
  if (center) {
    Ycenter <- apply(Y, 2, function(x) scale(x, center = TRUE, scale=FALSE))
  } else {
    Ycenter <- Y
  }
  
  m <- nrow(Y)
  n <- ncol(Y)
  
  #Compute the singular-value decomposition 
  svdWa <- svd(Ycenter[, cIdx])
  k <- min(k, max(which(svdWa$d > tolerance)))
  W <- svdWa$u[, factors, drop = FALSE]
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  correctedY <- Y - W %*% alpha
  if(!isLog && all(.isWholeNumber(x))) {
    if(round) {
      correctedY <- round(exp(correctedY) - epsilon)
      correctedY[correctedY<0] <- 0
    } else {
      correctedY <- exp(correctedY) - epsilon
    }
  }
  colnames(W) <- paste("W", factors, sep="_")
  return(list(W = W, normalizedCounts = t(correctedY)))
}

#helper function
.isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  !is.na(x) & abs(x - round(x)) < tol
}





## Remove unwanted variation from RNA-Seq data using RUVSeq --------------------------------------

#house-keeping genes used as negative controls
control_genes <- gene_attributes %>%
  subset(.$gene_name %in% c("TFRC", "SDHA", "TBP", "ACTB", "PPIA", "GUSB", "YWHAZ", "HMBS", "GAPDH", "RPLP0", "B2M", "RPL13A"))

#determine variance partition in gene expression after accounting for no, all or individual factors detected by RUVSeq
file_names <- c(1:12, "DESeq2", "all")
for(i in 1:length(file_names)) {
  
  file_name <- file_names[i]
  
  if(file_name != "DESeq2"){
    
    if(file_name == "all"){
      set_g <- RUVg_factor(counts(dds_counts), control_genes$Gene_ID, factors = 1:12, k = 12)
    }else if(file_name %in% 1:12){
      set_g <- RUVg_factor(counts(dds_counts), control_genes$Gene_ID, factors = i, k = 12)
    }
    
    rownames(set_g$W) <- colnames(counts(dds_counts))
    
    normalized_counts_RUVg <- set_g$normalizedCounts %>%
      as.data.frame() %>%
      add_column("Gene_ID" = rownames(.), .before = 1)
    write.table(normalized_counts_RUVg, paste0("results/RUVSeq/output/normalized_counts_RUVg_", file_name, ".txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    
    unwanted_factors_RUVg <- set_g$W %>%
      as.data.frame() %>%
      add_column("Sample_ID" = rownames(.), .before = 1)
    write.table(unwanted_factors_RUVg,  paste0("results/RUVSeq/output/unwanted_factors_RUVg_", file_name, ".txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    
    dds_RUVg <- DESeqDataSetFromMatrix(countData = set_g$normalizedCounts, colData = col_clinical_data, design = ~ 1)
    
  }else if(file_name == "DESeq2"){
    dds_RUVg <- dds_counts
  }
  
  dds_RUVg <- estimateSizeFactors(dds_RUVg)
  quantvst <- vst(dds_RUVg)
  quantvst <- assay(quantvst)
  
  form <- ~ (1|Individual_ID) + (1|Run) + (1|Annual_season) + (1|Season) + Time_of_sampling + (1|Infection_type) + (1|Vaccine_last_month) + (1|allergies_season) + (1|allergies_food) + Age..years. + (1|Gender) + BMI..kg.m2. + Hemoglobin..g.dL. + Neutrophils...µL. + Eosinophils...µL. + Lymphocytes...µL. + Monocytes...µL. + Thrombocytes..x1000.µL. + hsCRP..mg.L. + Hemoglobine.A1c...... + Creatinine..mg.dL. + Uric.acid..mg.dL. + Albumin..g.L. + Triglycerides..mg.dL.
  
  scaled_clinical_data <- col_clinical_data
  scaled_clinical_data[, c("Age..years.", "BMI..kg.m2.", "Hemoglobin..g.dL.", "Neutrophils...µL.", "Eosinophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L.", "Hemoglobine.A1c......", "Creatinine..mg.dL.", "Uric.acid..mg.dL.", "Albumin..g.L.", "Gamma.globulins..g.L.", "Triglycerides..mg.dL.", "Cholesterol.LDL..measured...mg.dL.",    "Bilirubin.total..mg.dL.", "Gamma.GT..U.L.", "CK..U.L.", "Time_of_sampling")] <- scale(scaled_clinical_data[, c("Age..years.", "BMI..kg.m2.", "Hemoglobin..g.dL.", "Neutrophils...µL.", "Eosinophils...µL.", "Lymphocytes...µL.", "Monocytes...µL.", "Thrombocytes..x1000.µL.", "hsCRP..mg.L.", "Hemoglobine.A1c......", "Creatinine..mg.dL.", "Uric.acid..mg.dL.", "Albumin..g.L.", "Gamma.globulins..g.L.", "Triglycerides..mg.dL.", "Cholesterol.LDL..measured...mg.dL.", "Bilirubin.total..mg.dL.", "Gamma.GT..U.L.", "CK..U.L.", "Time_of_sampling")])
  
  varPart <- fitExtractVarPartModel(quantvst, form, scaled_clinical_data)
  varMean <- colMeans(varPart)
  varPart_sorted <- varPart[, order(varMean, decreasing = TRUE)]
  write.table(varPart_sorted, paste0("results/RUVSeq/output/variance_partition_RUVg_", file_name, ".txt"), sep = '\t', quote = FALSE)
  
  varMean <- as.data.frame(varMean)
  varMean$Parameter <- rownames(varMean)
  varMean <- varMean[order(varMean$varMean, decreasing = TRUE), ]
  write.table(varMean, paste0("results/RUVSeq/output/variance_partition_mean_RUVg_", file_name, ".txt"), sep = '\t', quote = FALSE, row.names = FALSE)
}




## Plot variance explained by sequencing run and cohort membership after accounting for different RUVSeq factors as in Extended Data Fig. 2c ---------

#load data
df_ruvseq <- data.frame()
files <- c(paste0(rep("results/RUVSeq/output/variance_partition_RUVg_", 12), 1:12, ".txt"), 
           "results/RUVSeq/output/variance_partition_DESeq2.txt",
           "results/RUVSeq/output/variance_partition_RUVg_all.txt")

for (i in 1:14) {
  
  varPart <- read.table(files[i])
  varPart <- add_column(varPart, "ID" = paste0("Factor ", i))
  melted_data <- melt(varPart, id.vars = "ID", measure.vars = c("Individual_ID", "Run", "Residuals", "Season"))
  df_ruvseq <- rbind(df_ruvseq, melted_data)
  
}

#prepare table and rename variables
df_ruvseq$ID <- factor(df_ruvseq$ID, levels = c("Factor 13", "Factor 14", paste0(rep("Factor ", 12), 1:12)))
df_ruvseq$variable <- factor(df_ruvseq$variable, levels = c("Individual_ID", "Run", "Season", "Residuals"))
df_ruvseq$value <- df_ruvseq$value * 100

df_ruvseq <- df_ruvseq %>% subset(.$variable %in% c("Run", "Season"))
df_ruvseq$variable <- ifelse(df_ruvseq$variable == "Run", "Sequencing run", "Cohort")

ruvseq_boxplot <- ggplot(df_ruvseq, aes(x = ID,  y = value)) +
  geom_boxplot(aes(fill = variable)) +
  labs(y = "Variance explained (%)\nafter accounting for different\ncombinations of RUV factors") +
  scale_x_discrete(labels = c("Factor 13" = "No factors", "Factor 14" = "All factors")) + 
  scale_y_continuous(trans = "log1p") +
  scale_fill_manual(values = c("Sequencing run" = "darkturquoise", "Cohort" = "green4")) +
  facet_wrap(.~variable) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid = element_blank())





