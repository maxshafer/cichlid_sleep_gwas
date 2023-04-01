library(Matrix)
library(dplyr)
library(ggplot2)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(data.table)
library(org.Mm.eg.db)
library(patchwork)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/orthos")

## Read in the results from annotating SNPs
## Using a pvalue cutoff (of 0.0001) to ID significant SNPs
## Returns ~1k genes per comparison (only 216 for peaks)

genes1 <- readRDS(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.00025", "genes_ALL_SNPS.rds", sep = "_"))

## OK, I need to subset to unique hits
## Do I also care if they are the closest? Maybe not

genes <- lapply(genes1, function(x) x[x$feature == "gene",])
genes <- lapply(genes, function(x) x[x$closest == "yes",])
genes <- lapply(genes, function(x) {
  index <- grepl("LOC", x$gene)
  x$NCBI.GeneID <- NA
  x$NCBI.GeneID[index] <- substr(x$gene[index], start = 4, stop = 15)
  x$gene.name <- NA
  x$gene.name[!(index)] <- x$gene[!(index)]
  return(x)
})


new_gtf <- readRDS(file = "../orthos_gsea/gtf2_oreochromis_orthologs.rds")

## Add them back to the genes dataframes
## This selects the first entry in the ortho database as the ortho

genes_new <- list()
for (i in seq_along(genes)) {
  genes_new[[i]] <- genes[[i]]
  genes_new[[i]]$zeb_ortholog <- new_gtf$zebrafish_ortholog_symbol[match(genes_new[[i]]$Dbxref, new_gtf$Dbxref_sub)]
  genes_new[[i]]$mouse_ortholog <- new_gtf$mouse_ortholog_symbol[match(genes_new[[i]]$Dbxref, new_gtf$Dbxref_sub)]
}

names(genes_new) <- c("dn55", "dn58", "dn60", "tr", "peak")






splitNames <- function(x) {
  if (length(x) == 2) {
    names(x) <- c("danio", "mus")
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

hCluster <- function(x = go_analysis, measure.var = "Benj.value", category = c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY")) {
  casted <- acast(x[x$Category %in% category,], formula = L1 + L2 ~ Term, value.var = measure.var)
  casted[is.na(casted)] <- 0
  hr <- hclust(as.dist(1-cor(t(casted), method="pearson")), method="complete")
  hc <- hclust(as.dist(1-cor(casted, method="spearman")), method="complete")
  clustered <- x[x$Category %in% category,]
  clustered$Term <- factor(clustered$Term, levels = hc[[4]][hc[[3]]])
  clustered$L3 <- factor(clustered$L3, levels = hr[[4]][hr[[3]]])
  return(clustered)
}


## OK, for each of these data.frames, run DAVID analysis on the gene lists (zebrafish)

annoCharts <- list()
for (i in 1:length(genes_new)) {
  
  markers <- genes_new[[i]]$zeb_ortholog

  markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
  david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
  
  setAnnotationCategories(david, c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY"))
  
  result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(genes_new)[i], listType = "Gene")
  annoCharts[[i]] <- list()
  annoCharts[[i]][[1]] <- getFunctionalAnnotationChart(david)
  
  markers <- genes_new[[i]]$mouse_ortholog
  
  markers <- select(org.Mm.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
  david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
  
  setAnnotationCategories(david, c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY"))
  
  result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(genes_new)[i], listType = "Gene")
  annoCharts[[i]][[2]] <- getFunctionalAnnotationChart(david)
}


# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts) <- names(genes_new)
annoCharts <- lapply(annoCharts, function(x) splitNames(x))

list.index <- listList(annoCharts)
annoCharts.2 <- lapply(seq_along(annoCharts), function(x) annoCharts[[x]][unlist(list.index[[x]])])
names(annoCharts.2) <- names(genes_new)
annoCharts.2 <- lapply(annoCharts.2, function(x) splitNames(x))

saveRDS(annoCharts.2, file = "GO_analysis_all_pvalue_0.00025.rds")

# # If used subcluster markers, but grouped by subtype
annoCharts.2 <- readRDS("GO_analysis_all_pvalue_0.00025.rds")

# Reshape for plotting
measure.vars <- matrix(c("Bonferroni", "Benjamini", "FDR", "Count", "Fold.Enrichment", "X.", "Bonferroni", "Benjamini", "FDR", "Count", "FE", "X.", "Bonf.value", "Benj.value", "FDR.value", "Counts.value", "FoldE.value", "X.value"), nrow=6, ncol = 3)
melted <- apply(measure.vars, 1, function(x) reshape2::melt(annoCharts.2, id.vars = c("Category", "Term"), measure.vars = x[1], variable.name = x[2], value.name = x[3]))


go_analysis <- melted[[1]]
go_analysis$Benj.value <- melted[[2]]$Benj.value
go_analysis$FDR.value <- melted[[3]]$FDR.value
go_analysis$Count <- melted[[4]]$Counts.value
go_analysis$FoldE <- log(melted[[5]]$FoldE.value)
go_analysis$Perc <- melted[[6]]$X.value
go_analysis$FoldE[go_analysis$L2 == "mus"] <- go_analysis$FoldE[go_analysis$L2 == "mus"]*-1
go_analysis$L3 <- paste(go_analysis$L1, go_analysis$L2, sep = "_")

# Cluster with function, and plot
# GOTERM_BP_DIRECT
# UP_KEYWORDS
# GOTERM_MF_DIRECT
# GOTERM_CC_DIRECT
# KEGG_PATHWAY
# UP_SEQ_FEATURE
# REACTOME_PATHWAY
# BBID

go_analysis$L1 <- factor(go_analysis$L1, levels = names(genes_new))

go.plot <- ggplot(data = hCluster(x = go_analysis[go_analysis$L2 == "danio",], measure.var = "FDR.value", category = "KEGG_PATHWAY"), aes(L1, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point(alpha = 0.75) + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.title.y = element_blank()) + coord_flip()
go.plot <- ggplot(data = hCluster(x = go_analysis[go_analysis$L2 == "mus",], measure.var = "FDR.value", category = "KEGG_PATHWAY"), aes(L1, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point(alpha = 0.75) + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.title.y = element_blank()) + coord_flip()

go.plot.all <- ggplot(data = hCluster(x = go_analysis[go_analysis$Bonf.value < 0.05 & go_analysis$L2 == "danio",], measure.var = "Bonf.value"), aes(L1, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point(alpha = 0.75) + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.title.y = element_blank()) + coord_flip()


dev.new()
go.plot + plot_layout(width = unit(200, "mm"), height = unit(75, "mm"))

go.plot.all + plot_layout(width = unit(220, "mm"), height = unit(45, "mm"))


## To map back from an enriched category

test <- str_split(annoCharts.2[[4]]$danio[grep("dre04012:ErbB signaling pathway", annoCharts.2[[4]]$danio$Term),"Genes"], ", ")[[1]]

test <- select(org.Dr.eg.db, test, "SYMBOL", "ENTREZID")[,2]

# Find them in the og list to id oreochromis genes
View(genes_new[[4]][genes_new[[4]]$zeb_ortholog %in% test,])
