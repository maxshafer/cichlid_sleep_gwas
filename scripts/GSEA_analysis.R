








## Add them back to the genes dataframes
## This selects the first entry in the ortho database as the ortho

genes1 <- readRDS("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_10kb_windows_minPvalue.rds") ## Can also do meanPvalue

## OK, I need to subset to unique hits
## Maybe take the lowest pvalue per gene

genes <- lapply(genes1, function(x) x[x$feature == "gene",])
genes <- lapply(genes, function(x) x[x$closest == "yes",])
genes <- lapply(genes, function(x) x %>% group_by(Dbxref) %>% slice_min(order_by = ps))
genes <- lapply(genes, function(x) {
  index <- grepl("LOC", x$gene)
  x$NCBI.GeneID <- NA
  x$NCBI.GeneID[index] <- substr(x$gene[index], start = 4, stop = 15)
  x$gene.name <- NA
  x$gene.name[!(index)] <- x$gene[!(index)]
  return(x)
})

genes_new <- list()
for (i in seq_along(genes)) {
  
  genes_new[[i]] <- genes[[i]]
  genes_new[[i]]$zeb_ortholog <- new_gtf$zebrafish_ortholog_symbol[match(genes_new[[i]]$Dbxref, new_gtf$Dbxref_sub)]
}


## The below works, but should be run on ortholog lists (at least oreochromis -> zebrafish)


### Run some GSEA?
library(org.Dr.eg.db)
library(fgsea)
library(msigdbr)

## Load the go lists using package msigdbr, which is really nice https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
## This has go lists for zebrafish, which are associated with gene_symbols, entrez IDs, and ensembl IDs, as well as human versions (symbol, entrez, and ensembl) (presumably could use mouse to get mouse and human)

all_gene_sets = msigdbr(species = "Danio rerio")

## Make these into a list based on the gs_name
all_gene_sets_list <- split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)

## Or can subset to a specific kind of GO term
## gs_cat and gs_subcat (where "CP:KEGG" is)
kegg <- all_gene_sets[all_gene_sets$gs_subcat == "CP:KEGG",]
kegg_list <- split(x = kegg$gene_symbol, f = kegg$gs_name)

## Make a sorted gene list
gene_list <- genes_new[[4]]$ps
names(gene_list) <- as.character(genes_new[[4]]$zeb_ortholog)
gene_list <- gene_list[!(is.na(names(gene_list)))]

gene_list <- gene_list[!(duplicated(names(gene_list)))]
# Need to log transform the list (and inverse it) so that the results make sense
gene_list <- log(gene_list)*-1
gene_list <- gene_list[order(gene_list, decreasing = T)]


# Then can run gsea analysis using fgsea, the pathway (list), and gene_list
fgseaRes <- fgsea(pathways = kegg_list, gene_list, eps = 0)
# View(fgseaRes)

fgseaRes_all <- fgsea(pathways = all_gene_sets_list, gene_list, minSize=15, maxSize=500, scoreType = "pos")
# View(fgseaRes_all[fgseaRes_all$padj < 0.05,])

## Can look at top ones like this to find the names for plotting some example graphs
head(fgseaRes_all$pathway[order(fgseaRes_all$padj, decreasing = FALSE)])
head(fgseaRes$pathway[order(fgseaRes$padj, decreasing = FALSE)])

plotEnrichment(all_gene_sets_list[["KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY"]], gene_list) + labs(title="KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY")


### For gsea, the pathways object needs to be a named list of IDs that matches the names of the gene_list
### YES! THIS WORKS!






