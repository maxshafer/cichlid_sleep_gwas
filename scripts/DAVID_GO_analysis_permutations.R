library(Matrix)
library(dplyr)
library(ggplot2)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(data.table)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(patchwork)
library(here)
library(PhenoExam)
library(stringr)

splitNames <- function(x) {
  if (length(x) == 3) {
    names(x) <- c("human", "mouse", "zebrafish")
  }
  return(x)
}

listList <- function(x) {
  list <- list()
  for (i in 1:length(x)){
    list[[i]] <- lapply(x[[i]], function(y) length(y) > 0)
  }
  return(list)
}

source("scripts/variants_functions.R")


## Load in results from orthology, it includes lists for human, mouse, and zebrafish, for PC1/PC2/TR
## Also have gene names and ensembl IDs

permuted.genes <- readRDS(file = "Permute_SNPs_perchr_1000_permutations_800_SNPs_orthologs.rds")
ortholog_df <- readRDS("GTF_human_orthologs_all.rds")

## OK, for each of these data.frames, run DAVID analysis on the gene lists
ortholog_df <- ortholog_df[ortholog_df$human_ortholog != "",]
# ortholog_df <- ortholog_df[!duplicated(ortholog_df$human_ortholog),]

dbi_out <- select(org.Hs.eg.db, ortholog_df$human_ortholog, "ENTREZID", "SYMBOL")

ortholog_df$human_ENTREZID <- dbi_out$ENTREZID[match(ortholog_df$human_ortholog, dbi_out$SYMBOL)]

annoCharts <- lapply(permuted.genes, function(x) {
  markers <- ortholog_df$human_ENTREZID[match(unique(x$human_ortholog), ortholog_df$human_ortholog)]
  # markers <- select(org.Hs.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
  david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
  setAnnotationCategories(david, c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY", "DISGENET", "GAD_DISEASE", "UP_TISSUE"))
  result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = "permutated_snps", listType = "Gene")
  result <- addList(david, dbi_out$ENTREZID, idType = "ENTREZ_GENE_ID", listName = "backgroundorthos", listType = "Background")
  out <- getFunctionalAnnotationChart(david, threshold = 1)
  return(out)
})
## Temp save out, so I don't have to run this again
saveRDS(annoCharts, file = "DAVID-GO_analysis_Permute_SNPs_perchr_1000_permutations_800_SNPs_orthoBKGRD.rds")

### Use PhenoExam to get MGI and others
pop.totals <- t(data.frame(MGD = 17895, PSYGENET = 19262, HPO = 19248, CGI = 19198, CLINGEN = 19198, CRB = 19274, CTD = 19636, GENOMICS_ENGLAND = 19230, ORPHANET = 19262, UNIPROT = 19204))

pheno_exam <- lapply(permuted.genes, function(x) {
  PE_human <- PhenoEnrichGenes(genes= unique(x$human_ortholog), database = getdbnames())
  
  PE_human <- PE_human$alldata
  PE_human$term_id <- str_replace(PE_human$term_id, ":", "-")
  PE_human$Term <- paste(PE_human$term_id, PE_human$term_name, sep = ":")
  colnames(PE_human) <- c("term_id", "term_name", "Category", "Bonferroni", "FDR", "Pop.Hits", "Count", "overlap_ratio", "PValue", "Genes", "Term")
  PE_human$List.Total <- length(unique(x$human_ortholog))
  
  PE_human$Pop.Total <- pop.totals[,1][match(PE_human$Category, rownames(pop.totals))]
  PE_human$Fold.Enrichment <- (PE_human$Count / PE_human$List.Total) / (PE_human$Pop.Hits / PE_human$Pop.Total)
  PE_human$Benjamini <- NA
  PE_human$X. <- PE_human$Count / PE_human$Pop.Hits
  
  PE_human <- PE_human[, c("Category","Term","Count","X.","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")]
  return(PE_human)
})


#### Combine and save out
permuted_go_analysis <- lapply(seq_along(permuted.genes), function(x) rbind(annoCharts[[x]], pheno_exam[[x]]))

saveRDS(permuted_go_analysis, file = "DAVID-GO_analysis_Permute_SNPs_perchr_1000_permutations_800_SNPs.rds")


## Read back in to analyse
permuted_go_analysis <- readRDS("DAVID-GO_analysis_Permute_SNPs_perchr_1000_permutations_800_SNPs.rds")

annoCharts.2 <- readRDS("DAVID-GO_analysis_perchr_60sp_1e-05_percentile.rds")



permute_out <- lapply(annoCharts.2, function(ac2) {
  
  test_go <- ac2$human[ac2$human$PValue < 0.05,]
  
  ## This returns the matching pvalues for each permutation
  result <- apply(test_go, 1, function(x) {
    
    test <- unlist(lapply(permuted_go_analysis, function(y) {
      y <- y[y$Category %in% x[1] & y$Term %in% x[2], "PValue"][1]
      return(y)
    }))
    
    return(test)
  })
  
  # This then tests if they are less than the pvalue from our test
  # Each row is a permutation, each column a category
  result_pvalues_t <- t(result)
  
  new_result <- lapply(seq_along(1:nrow(result_pvalues_t)), function(x) {
    return(table(result_pvalues_t[x,] < test_go[x,"PValue"]))
  })
  
  
  test_go$permute_result <- unlist(lapply(new_result, function(x) prop.table(x)[1]))
  
  return(test_go)
  
})

saveRDS(permute_out, file = "Permutation_corrected_GO_analysis.rds")

permute_out_fdr <- lapply(annoCharts.2, function(ac2) {
  
  test_go <- ac2$human[ac2$human$FDR < 0.05,]
  
  ## This returns the matching pvalues for each permutation
  result <- apply(test_go, 1, function(x) {
    
    test <- unlist(lapply(permuted_go_analysis, function(y) {
      y <- y[y$Category %in% x[1] & y$Term %in% x[2], "PValue"][1]
      return(y)
    }))
    
    return(test)
  })
  
  # This then tests if they are less than the pvalue from our test
  # Each row is a permutation, each column a category
  result_pvalues_t <- t(result)
  
  new_result <- lapply(seq_along(1:nrow(result_pvalues_t)), function(x) {
    return(table(result_pvalues_t[x,] < test_go[x,"PValue"]))
  })
  
  
  test_go$permute_result <- unlist(lapply(new_result, function(x) prop.table(x)[1]))
  
  return(test_go)
  
})




## Somehow the majority of the random GO tests always have lower p-values, like 8/9/10 out of 10 permuted gene sets
## This doesn't make sense to me, and it is even when controlling for the background
## I think there is something wrong with how the permuted sets were chosen? I don't have any other explanation
## Unless there is something different about the first few
## OK I figured it out, my code to check if the pvalues were lower wasn't working correctly (I should say, the function didn't seem to be working properly)
## I still don't know how to use a different background set in PhenoExam (I could do it online, but that is impossible for 1000x gene sets)


