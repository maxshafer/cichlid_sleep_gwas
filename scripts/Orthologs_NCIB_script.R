library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)
library(fdrtool)
library(stringr)
library(rjson)
library(clipr)

library(clusterProfiler)
library(AnnotationDbi)
library(org.Dr.eg.db)

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R")
# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/orthos")


################################################################################################################
################ Generate lists of genes for finding orthologs #################################################
################################################################################################################

## Read in the results from annotating SNPs
## Using a pvalue cutoff (of 0.0001) to ID significant SNPs
## Returns ~1k genes per comparison (only 216 for peaks)

genes1 <- readRDS(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "genes_ALL_SNPS.rds", sep = "_"))

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

## Ok, for each of these lists, need to output a file without header with the IDs (Dbxref)

write.table(genes[[1]]$Dbxref, file = "perchr_dn55_genes.csv", col.names = F, row.names = F, quote = F)
write.table(genes[[2]]$Dbxref, file = "perchr_dn58_genes.csv", col.names = F, row.names = F, quote = F)
write.table(genes[[3]]$Dbxref, file = "perchr_dn60_genes.csv", col.names = F, row.names = F, quote = F)
write.table(genes[[4]]$Dbxref, file = "perchr_tr_genes.csv", col.names = F, row.names = F, quote = F)
write.table(genes[[5]]$Dbxref, file = "perchr_peak_genes.csv", col.names = F, row.names = F, quote = F)


################################################################################################################
################ Run terminal commands to download NCBI orthologs ##############################################
################################################################################################################

## Need to run these system commands for all of the lists (within the directories)
## Makes named folders for each Niloticus gene and puts the ortho inside (as a csv)
## Runs in terminal, not R (without single-quotes)

## Mouse taxon ID 10090
## Zebrafish taxon ID 7955
## Tilapia taxon ID 8128

'cat perchr_tr_genes.csv | while read GENE; do
datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 10090 "${GENE}" --filename "tr/mouse/${GENE}".zip;
done'

# Then to unzip all the .zip files in a folder, then the zip files can be deleted
'for f in *.zip; do unzip "$f" -d "${f%.zip}"; done'

## This one needs to recursively find all the paths
## Then load and extract the ortholog from the csv
## This output has duplicate entries (because of transcript variants), which are reconciled when combining with orthoDB

# For zeb
directories_zeb <- c("dn55/zeb/", "dn58/zeb/", "dn60/zeb/", "tr/zeb/", "peak/zeb/")
ncbi_orthologs_zeb <- lapply(directories_zeb, function(x) generateNCBIorthologs(directory = x))
names(ncbi_orthologs_zeb) <- c("dn55", "dn58", "dn60", "tr", "peak")
ncbi_orthologs_zeb <- lapply(ncbi_orthologs_zeb, function(x) Reduce(rbind, x))

# For mouse
directories_mouse <- c("dn55/mouse/", "dn58/mouse/", "dn60/mouse/", "tr/mouse/", "peak/mouse/")
ncbi_orthologs_mouse <- lapply(directories_mouse, function(x) generateNCBIorthologs(directory = x))
names(ncbi_orthologs_mouse) <- c("dn55", "dn58", "dn60", "tr", "peak")
ncbi_orthologs_mouse <- lapply(ncbi_orthologs_mouse, function(x) Reduce(rbind, x))


################################################################################################################
################ Query orthoDB database and ID orthologs from output ###########################################
################################################################################################################

## This one needs to run a bunch of system commands to query the orthoDB database
## Takes a Dbxref ID as input, and outputs a list of vectors (matching orthologs)
## These are better then the above, but miss some and should be supplemented
## This can be l/apply'd across a list of gene_ids

## This creates some unnecessary structure, but works
## Need to rework this function so I can run it for mouse

orthoDB_orthologs_zeb <- lapply(genes, function(x) lapply(x[,"Dbxref"], function(y) generateOrthoDBorthologs(dbxrefs = y, species = "zebrafish")))
names(orthoDB_orthologs_zeb) <- c("dn55", "dn58", "dn60", "tr", "peak")

orthoDB_orthologs_mouse <- lapply(genes, function(x) lapply(x[,"Dbxref"], function(y) generateOrthoDBorthologs(dbxrefs = y, species = "mouse")))
names(orthoDB_orthologs_mouse) <- c("dn55", "dn58", "dn60", "tr", "peak")

# Unlist the unnecessary structure by using a non-recursive unlist
orthoDB_orthologs_zeb <- lapply(orthoDB_orthologs_zeb, function(x) unlist(x, recursive = F))
orthoDB_orthologs_mouse <- lapply(orthoDB_orthologs_mouse, function(x) unlist(x, recursive = F))


################################################################################################################
################ Reconcile and combine the two database results together #######################################
################################################################################################################

### OK, so the ncbi ortho  and orthoDB pipelines retrieve slightly different results (some missing orthos)
### There are X that ncbi retreives that orthoDB doesn't, and X vice versa - combined is more
## I can add the ncbi_orthologs to the orthoDB_orthologs list

orthos_combined_zeb <- orthoDB_orthologs_zeb
for (j in 1:length(orthos_combined_zeb)) {
  for (i in 1:nrow(ncbi_orthologs_zeb[[j]])) {
    if (!(ncbi_orthologs_zeb[[j]]$gene_id_oreo[i] %in% names(orthos_combined_zeb[[j]]))) {
      orthos_combined_zeb[[j]][[ncbi_orthologs_zeb[[j]]$gene_id_oreo[i]]] <- list(ncbi_orthologs_zeb[[j]]$gene_symbol[i])
    }
  }
}

orthos_combined_mouse <- orthoDB_orthologs_mouse
for (j in 1:length(orthos_combined_mouse)) {
  for (i in 1:nrow(ncbi_orthologs_mouse[[j]])) {
    if (!(ncbi_orthologs_mouse[[j]]$gene_id_oreo[i] %in% names(orthos_combined_mouse[[j]]))) {
      orthos_combined_mouse[[j]][[ncbi_orthologs_mouse[[j]]$gene_id_oreo[i]]] <- list(ncbi_orthologs_mouse[[j]]$gene_symbol[i])
    }
  }
}

## Save out the list

saveRDS(orthos_combined_zeb, file = "oreochromis_orthologs_zebrafish.rds")
saveRDS(orthos_combined_mouse, file = "oreochromis_orthologs_mouse.rds")



## Add them back to the genes dataframes
## This selects the first entry in the ortho database as the ortho

genes_new <- list()
for (i in seq_along(genes)) {
  
  zeb_ortholog <- unlist(lapply(orthos_combined_zeb[[i]], function(x) {
    if (length(x) > 0) {
      return(x[[1]][[1]][1])
    } else {
      return(NA)
    }}))
  zeb_ortholog <- data.frame(Dbxref = names(zeb_ortholog), zeb_ortholog = zeb_ortholog)
  
  mouse_ortholog <- unlist(lapply(orthos_combined_mouse[[i]], function(x) {
    if (length(x) > 0) {
      return(x[[1]][[1]][1])
    } else {
      return(NA)
    }}))
  mouse_ortholog <- data.frame(Dbxref = names(mouse_ortholog), mouse_ortholog = mouse_ortholog)
  genes_new[[i]] <- genes[[i]]
  genes_new[[i]]$zeb_ortholog <- zeb_ortholog$zeb_ortholog[match(genes_new[[i]]$Dbxref, zeb_ortholog$Dbxref)]
  genes_new[[i]]$mouse_ortholog <- mouse_ortholog$mouse_ortholog[match(genes_new[[i]]$Dbxref, mouse_ortholog$Dbxref)]
}


### Run some GSEA?
library(org.Mm.eg.db)

gene_list <- genes_new[[1]]$ps
names(gene_list) <- as.character(genes_new[[1]]$zeb_ortholog)
gene_list <- gene_list[order(gene_list, decreasing = T)]
gene_list <- gene_list[!(is.na(names(gene_list)))]


gse <- gseGO(gene_list, ont = "BP", keyType = "SYMBOL", OrgDb = org.Dr.eg.db, scoreType = "pos", eps = 1e-300)






# OK so now I have a list of danio IDs that correspond to the oreochromis ones, and I'm only missing ncRNAs or things that don't have an ortholog in zebrafish!
# I should also be able to pull mouse orthologs (directly from niloticus)

missing_ids <- genes[[1]]$Dbxref[!(genes[[1]]$Dbxref %in% names(orthos_combined_zeb))]


## One issue that comes to mind - if a gene has two orthologs, they are likely both in the same pathway, and therefore it increases the chances of enrichment of pathways
## There is also a chance that only one of the orthologs is in a pathway, but not sure - I suppose thats why this isn't common practice...
# There are also duplicates, which I think doesn't affect it?
write_clip(unique(unlist(orthos_combined_zeb$tr)))

write_clip(c(unique(unlist(orthos_combined_zeb$dn55)), unique(unlist(orthos_combined_zeb$dn58)), unique(unlist(orthos_combined_zeb$dn60))))

## Can also backtrack to ID the oreo gene id that is associated with a zeb/mouse gene
## The location of the SNP is in the genes1 list of dfs
unlist(danio_orthos_combined$dn60)[grep("gnrhr4", unlist(danio_orthos_combined$dn60))]
















