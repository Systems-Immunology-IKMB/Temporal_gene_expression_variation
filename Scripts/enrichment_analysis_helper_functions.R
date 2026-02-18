library(genefilter)
library(geneplotter)
library(topGO)
library(DBI)
library(plyr)

## Run GO enrichment analysis using genes with a similar expression level as universe --------------

#This function takes 
#a set of selected genes "degs", 
#a vector containing the average expression of all genes in the analysis "overallBaseMean",
#the database selected for mapping (usually "org.Hs.eg.db" or "org.Mm.eg.db"),
#a data.frame "annotation" with a gene to symbol mapping ("Gene_ID" and "gene_name" as column names),
#the number of background genes "back_num" to select for each gene in "degs",
#the minimal node size "min_node_size",
#a boolean same_expression_level indicating if ORA should be done with the universe consisting of genes with a similar expression level as the selected genes
run_GO_ORA <- function(degs, overallBaseMean, selected_database, annotation, same_expression_level = TRUE, back_num = 10, min_node_size = 5){
  
  #Identify genes with similar gene expression level
  sig_idx <- match(degs, rownames(overallBaseMean))
  backG <- c()
  
  for(i in sig_idx){
    ind <- genefinder(overallBaseMean, i, back_num, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
  }
  
  backG <- unique(backG)
  backG <- rownames(overallBaseMean)[backG]
  backG <- setdiff(backG, degs)
  
  multidensity(list(all = log2(overallBaseMean[,"mean_expression"]+1),
                    foreground = log2(overallBaseMean[degs, "mean_expression"]+1), 
                    background = log2(overallBaseMean[backG, "mean_expression"]+1)), 
               xlab ="log2 mean normalized counts", 
               main = "Matching for enrichment analysis")
  
  #GO enrichment analysis
  onts <- c("MF", "BP", "CC")
  
  geneIDs <- rownames(overallBaseMean)
  
  inUniverse <- geneIDs %in% c(degs, backG) 
  if(!same_expression_level){
    inUniverse <- geneIDs %in% geneIDs
  }
  
  inSelection <- geneIDs %in% degs 
  alg <- factor(as.integer(inSelection[inUniverse]))
  names(alg) <- geneIDs[inUniverse]
  
  tab <- as.list(onts)
  names(tab) <- onts
  i <- 2
  
  for(i in 1:3){
    
    tgd <- new("topGOdata", 
               ontology = onts[i], 
               allGenes = alg, 
               nodeSize = min_node_size,
               annot = annFUN.org, 
               mapping = selected_database, 
               ID = "ENSEMBL")
    
    resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
    
    if(length(nodes(graph(tgd))) < 200){
      table_temp <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                              Fisher.classic = resultTopGO.classic,
                              orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
    }else{
      table_temp <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                              Fisher.classic = resultTopGO.classic,
                              orderBy = "Fisher.elim" , topNodes = 200)
    }
    
    #add gene symbols
    my_genes <- table_temp$GO.ID %>%
      genesInTerm(object = tgd)
    
    table_temp$Genes <- my_genes %>%
      map_chr(.f = function(x){intersect(x, sigGenes(tgd)) %>% 
          paste(collapse = ", ") %>% 
          return()})
    
    table_temp$Symbols <- my_genes %>%
      map_chr(.f = function(x){annotation %>%
          subset(.$Gene_ID %in% intersect(x, sigGenes(tgd))) %>%
          .$gene_name %>% 
          paste(collapse = ", ") %>%
          return()})
    
    table_temp$Fisher.elim <- as.numeric(table_temp$Fisher.elim)
    table_temp$Fisher.classic <- as.numeric(table_temp$Fisher.classic)
    table_temp$Term <- getTermsDefinition(whichTerms = table_temp$GO.ID, onts[i])
    
    tab[[i]] <- table_temp
  }
  
  tab$MF$ont <- "MF"
  tab$BP$ont <- "BP"
  tab$CC$ont <- "CC"
  topGOResults <- rbind.fill(tab)
  
  return(topGOResults)
}




## Filter GO results for significance, minimum number of genes, ontology, unique gene sets and maximum number of terms -------------

filter_go_results <- function(go_results, top_num, onto = "ALL"){
  
  go_results <- go_results %>%
    subset(.$Fisher.elim < 0.05) %>%
    subset(.$Significant > 1)
  
  if(onto != "ALL"){
    go_results <- subset(go_results, go_results$ont == onto)
  }
  
  if (nrow(go_results) > 0) {
    go_results <- go_results[order(go_results$Fisher.elim), ]
    exclude_indices <- c()
    unique_gene_sets <- c()
    for(i in 1:length(go_results$GO.ID)){
      if(as.vector(go_results$Genes[i]) %in% unique_gene_sets){
        exclude_indices <- c(exclude_indices, i)
      }else{
        unique_gene_sets <- c(unique_gene_sets, as.vector(go_results$Genes[i]))
      }
    }
    
    if(length(exclude_indices) > 0){
      go_results_filtered <- go_results[-exclude_indices, ]
    }else{
      go_results_filtered <- go_results
    }
    
    if(length(go_results_filtered$GO.ID) > top_num){
      top <- go_results_filtered[1:top_num,]
    }else{
      top <- go_results_filtered
    }
    
    top$Fisher.elim
    top$p <- -1*log10(as.numeric(top$Fisher.elim))
    top$Ratio <- top$Significant / top$Annotated
    
    return(top)
  }
}


## Get full GO term names -----------------

getTermsDefinition <- function(whichTerms, ontology, numChar = 1000, multipLines = FALSE) {
  
  qTerms <- paste(paste("'", whichTerms, "'", sep = ""), collapse = ",")
  retVal <- dbGetQuery(GO_dbconn(), paste("SELECT term, go_id FROM go_term WHERE ontology IN",
                                          "('", ontology, "') AND go_id IN (", qTerms, ");", sep = ""))
  
  termsNames <- retVal$term
  names(termsNames) <- retVal$go_id
  
  if(!multipLines)
    shortNames <- paste(substr(termsNames, 1, numChar),
                        ifelse(nchar(termsNames) > numChar, '...', ''), sep = '')
  else
    shortNames <- sapply(termsNames,
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\\\n"))
                         })
  
  names(shortNames) <- names(termsNames)
  
  #return NAs for the terms that are not found in the DB and make sure the 'names' attribute is as specified
  shortNames <- shortNames[whichTerms]
  names(shortNames) <- whichTerms
  
  return(shortNames)
}




## Get data.frame object that can be used for GO term dotplot -------------------------

get_GO_df <- function(GO_tables, direction, score_name, reverse_bool = FALSE, onto = "BP", showTerms = 10, multLines = TRUE, numChar = 60){
  
  #get top terms 
  table <- GO_tables
  score_index <- which(names(table[[1]]) == score_name)
  i <- 1
  df_terms <- c()
  for(i in 1:length(GO_tables)) {
    colnames(table[[i]])[score_index] <- "Scores"
    
    if(onto != "ALL"){
      table[[i]] <- table[[i]] %>%
        subset(.$ont == onto)
    }
    
    #Order data frame by P-value
    idx <- order(table[[i]]$Scores, decreasing = FALSE)
    table[[i]] <- table[[i]][idx, ]
    
    ifelse(is.numeric(showTerms), 
           table[[i]] <- table[[i]][1:showTerms, ],
           table[[i]] <- table[[i]] %>% subset(.$GO.ID %in% showTerms))
    df_terms <- union(df_terms, table[[i]]$GO.ID)
  }
  
  df <- data.frame()
  for(i in 1:length(GO_tables)) {
    df_temp <- GO_tables[[i]][GO_tables[[i]]$GO.ID %in% df_terms, ]
    df_temp <- df_temp %>%
      add_column(Direction = i)
    
    df <- rbind(df, df_temp)
  }
  
  if(reverse_bool){
    my_levels <- table %>%
      lapply(function(x){x$Term %>% rev()}) %>%
      c(recursive = TRUE) %>%
      unique()
  }else if(!reverse_bool){
    my_levels <- table %>%
      lapply(function(x){x$Term}) %>%
      c(recursive = TRUE) %>%
      unique() %>%
      rev()
  }
  
  df$Term <- factor(df$Term, levels = my_levels)
  colnames(df)[7] <- "Scores"
  df$Direction <- factor(df$Direction)
  levels(df$Direction) <- direction
  
  #reduce length of term names 
  termsNames <- levels(df$Term)
  
  if(multLines == FALSE) {
    shortNames <- paste(substr(as.character(termsNames), 1, numChar),
                        ifelse(nchar(as.character(termsNames)) > numChar, '...', ''), sep = '')
  } else {
    shortNames <- sapply(as.character(termsNames),
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\n"))
                         })
  }
  
  levels(df$Term) <- shortNames

  df <- df %>% 
    add_column(changed_scores = -1*log10(df$Scores))
  
  return(df)
}


