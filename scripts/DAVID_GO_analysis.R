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

genes <- readRDS(file = "orthos/orthologs_perchr_1e-05_percentile.rds")

genes_new <- genes$gene_names


## OK, for each of these data.frames, run DAVID analysis on the gene lists (zebrafish)

annoCharts <- list(list(), list(), list())
for (i in 1:length(genes_new)) {
  
  markers <- genes_new[[i]]$human
  markers <- select(org.Hs.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
  david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
  setAnnotationCategories(david, c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY", "DISGENET", "GAD_DISEASE", "UP_TISSUE"))
  result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(genes_new)[i], listType = "Gene")
  annoCharts[[i]][[1]] <- getFunctionalAnnotationChart(david, threshold = 0.05)
  
  markers <- genes_new[[i]]$mouse
  markers <- select(org.Mm.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
  david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
  setAnnotationCategories(david, c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY", "DISGENET", "GAD_DISEASE", "UP_TISSUE"))
  result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(genes_new)[i], listType = "Gene")
  annoCharts[[i]][[2]] <- getFunctionalAnnotationChart(david, threshold = 0.05)
  
  markers <- genes_new[[i]]$zebrafish
  markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
  david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
  setAnnotationCategories(david, c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY", "DISGENET", "GAD_DISEASE", "UP_TISSUE"))
  result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(genes_new)[i], listType = "Gene")
  annoCharts[[i]][[3]] <- getFunctionalAnnotationChart(david, threshold = 0.05)
  
}


# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts) <- names(genes_new)
annoCharts <- lapply(annoCharts, function(x) splitNames(x))

list.index <- listList(annoCharts)
annoCharts.2 <- lapply(seq_along(annoCharts), function(x) annoCharts[[x]][unlist(list.index[[x]])])
names(annoCharts.2) <- names(genes_new)
annoCharts.2 <- lapply(annoCharts.2, function(x) splitNames(x))


## Read in and append MGI phenotypes in orthos/modPhEA

modPhEA <- list(list(list()), list(list()), list(list()))
names(modPhEA) <- c("pc1", "pc2", "tr")

modPhEA$pc1[[1]] <- readModPhEA(file = "orthos/modPhEA/modPhEA_PC1_oreochromis-human_vs_mouse.txt", pvalue = "FishersExact")
modPhEA$pc1[[2]] <- readModPhEA(file = "orthos/modPhEA/modPhEA_PC1_oreochromis-mouse_vs_mouse.txt", pvalue = "FishersExact")
names(modPhEA$pc1) <- c("human", "mouse")
modPhEA$pc2[[1]] <- readModPhEA(file = "orthos/modPhEA/modPhEA_PC2_oreochromis-human_vs_mouse.txt", pvalue = "FishersExact")
modPhEA$pc2[[2]] <- readModPhEA(file = "orthos/modPhEA/modPhEA_PC2_oreochromis-mouse_vs_mouse.txt", pvalue = "FishersExact")
names(modPhEA$pc2) <- c("human", "mouse")
modPhEA$tr[[1]] <- readModPhEA(file = "orthos/modPhEA/modPhEA_TR_oreochromis-human_vs_mouse.txt", pvalue = "FishersExact")
modPhEA$tr[[2]] <- readModPhEA(file = "orthos/modPhEA/modPhEA_TR_oreochromis-mouse_vs_mouse.txt", pvalue = "FishersExact")
names(modPhEA$tr) <- c("human", "mouse")


## Append to annocharts
## This works, but need to combine with annocharts in some way? As it is only some of them?

annoCharts.2$pc1$human <- rbind(annoCharts.2$pc1$human, modPhEA$pc1$human)
annoCharts.2$pc1$mouse <- rbind(annoCharts.2$pc1$mouse, modPhEA$pc1$mouse)
annoCharts.2$pc2$human <- rbind(annoCharts.2$pc2$human, modPhEA$pc2$human)
annoCharts.2$pc2$mouse <- rbind(annoCharts.2$pc2$mouse, modPhEA$pc2$mouse)
annoCharts.2$tr$human <- rbind(annoCharts.2$tr$human, modPhEA$tr$human)
annoCharts.2$tr$mouse <- rbind(annoCharts.2$tr$mouse, modPhEA$tr$mouse)

saveRDS(annoCharts.2, file = "DAVID-GO_analysis_perchr_1e-05_percentile.rds")










