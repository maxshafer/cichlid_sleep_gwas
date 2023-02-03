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


# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas")

## Read in the results from annotating SNPs
## Using a pvalue cutoff (of 0.0001) to ID significant SNPs
## Returns ~1k genes per comparison (only 216 for peaks)

genes1 <- readRDS(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "genes_ALL_SNPS.rds", sep = "_"))

# genes2 <- lapply(genes1, function(x) x[x$distance_to_gene == 0,])
# lapply(genes2, function(x) length(unique(x[x$distance_to_gene == 0, "location"])))
# 
# ## Even fewer in exons at this point (don't need here)
# exon_genes_2 <- readRDS(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "SNPSinExons.csv", sep = "_"))

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

# ## Ok, for each of these lists, need to output a file without header with the IDs (Dbxref)
# 
# write.table(genes[[1]]$Dbxref, file = "perchr_dn55_genes.csv", col.names = F, row.names = F, quote = F)
# write.table(genes[[2]]$Dbxref, file = "perchr_dn58_genes.csv", col.names = F, row.names = F, quote = F)
# write.table(genes[[3]]$Dbxref, file = "perchr_dn60_genes.csv", col.names = F, row.names = F, quote = F)
# write.table(genes[[4]]$Dbxref, file = "perchr_tr_genes.csv", col.names = F, row.names = F, quote = F)
# write.table(genes[[5]]$Dbxref, file = "perchr_peak_genes.csv", col.names = F, row.names = F, quote = F)

## Run through each


## OK this for loop works a bit better, makes named folders for each Niloticus gene and puts the ortho inside

## Mouse 10090
## Zebrafish 7955
## Oreochromis niloticus 8128

'cat perchr_tr_genes.csv | while read GENE; do
datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 10090 "${GENE}" --filename "tr/mouse/${GENE}".zip;
done'

'for f in *.zip; do unzip "$f" -d "${f%.zip}"; done'

## OK, for each of the sets, I need to write functions for the following

## This one needs to recursively find all the paths
## Then load and extract the ortholog from the csv
## This output has duplicate entries (because of transcript variants)

# # For zeb
# directories <- c("dn55/zeb/", "dn58/zeb/", "dn60/zeb/", "tr/zeb/", "peak/zeb/")

# For mouse
directories <- c("dn55/mouse/", "dn58/mouse/", "dn60/mouse/", "tr/mouse/", "peak/mouse/")


ncbi_orthologs <- lapply(directories, function(x) generateNCBIorthologs(directory = x))
names(ncbi_orthologs) <- c("dn55", "dn58", "dn60", "tr", "peak")

## This now has duplicates (because it returns all of the transcript variants)
ncbi_orthologs <- lapply(ncbi_orthologs, function(x) Reduce(rbind, x))


## This one needs to run a bunch of system commands to query the orthoDB database
## Takes a Dbxref ID as input, and outputs a list of vectors (matching orthologs)
## These are better then the above, but miss some and should be supplemented
## This can be l/apply'd across a list of gene_ids

## lapply, lapply loop?
## This creates some unnecessary structure, but works

orthoDB_orthologs <- lapply(genes, function(x) lapply(x[,"Dbxref"], function(y) generateOrthoDBorthologs(dbxrefs = y)))
names(orthoDB_orthologs) <- c("dn55", "dn58", "dn60", "tr", "peak")

# Unlist the unnecessary structure by using a non-recursive unlist
orthoDB_orthologs <- lapply(orthoDB_orthologs, function(x) unlist(x, recursive = F))




### OK, so the ncbi ortho  and orthoDB pipelines retrieve slightly different results (some missing orthos)
### There are X that ncbi retreives that orthoDB doesn't, and X vice versa - combined is more
## I can add the ncbi_orthologs to the orthoDB_orthologs list


danio_orthos_combined <- orthoDB_orthologs
for (j in 1:length(danio_orthos_combined)) {
  
  for (i in 1:nrow(ncbi_orthologs[[j]])) {
    if (!(ncbi_orthologs[[j]]$gene_id_oreo[i] %in% names(danio_orthos_combined[[j]]))) {
      danio_orthos_combined[[j]][[ncbi_orthologs[[j]]$gene_id_oreo[i]]] <- list(ncbi_orthologs[[j]]$gene_symbol[i])
    }
  }
  
}


## Save out the list

# saveRDS(danio_orthos_combined, file = "oreochromis_orthologs_zebrafish.rds")
saveRDS(danio_orthos_combined, file = "oreochromis_orthologs_mouse.rds")

# OK so now I have a list of danio IDs that correspond to the oreochromis ones, and I'm only missing ncRNAs or things that don't have an ortholog in zebrafish!
# I should also be able to pull mouse orthologs (directly from niloticus)

missing_ids <- genes[[1]]$Dbxref[!(genes[[1]]$Dbxref %in% names(danio_orthos_combined))]


## One issue that comes to mind - if a gene has two orthologs, they are likely both in the same pathway, and therefore it increases the chances of enrichment of pathways
## There is also a chance that only one of the orthologs is in a pathway, but not sure - I suppose thats why this isn't common practice...
# There are also duplicates, which I think doesn't affect it?
write_clip(unique(unlist(danio_orthos_combined$tr)))

write_clip(c(unique(unlist(danio_orthos_combined$dn55)), unique(unlist(danio_orthos_combined$dn58)), unique(unlist(danio_orthos_combined$dn60))))

## Can also backtrack to ID the oreo gene id that is associated with a zeb/mouse gene
## The location of the SNP is in the genes1 list of dfs
unlist(danio_orthos_combined$dn60)[grep("gnrhr4", unlist(danio_orthos_combined$dn60))]
















