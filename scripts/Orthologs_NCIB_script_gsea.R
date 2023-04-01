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
library(tictoc)

library(clusterProfiler)
library(AnnotationDbi)
library(org.Dr.eg.db)

source("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R")
# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/orthos_gsea")
# setwd("~/Documents/orthos_gsea/")

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]

## Ok, for each of these lists, need to output a file without header with the IDs (Dbxref)

index <- grepl("LOC", gtf2$Name)
gtf2$NCBI.GeneID <- NA
gtf2$NCBI.GeneID[index] <- substr(gtf2$gene[index], start = 4, stop = 15)
gtf2$gene.name <- NA
gtf2$gene.name[!(index)] <- gtf2$gene[!(index)]


gtf2$Dbxref_sub <- substr(gtf2$Dbxref, start = 8, stop = 20)


write.table(unique(gtf2$Dbxref_sub), file = "unique_gtf2_genes.csv", col.names = F, row.names = F, quote = F)

################################################################################################################
################ Run terminal commands to download NCBI orthologs ##############################################
################################################################################################################

## I should probably just run the following for all possible genes, make a data.frame, then use that to get orthos!

## Makes named folders for each Niloticus gene and puts the ortho inside (as a csv)
## Runs in terminal, not R (without single-quotes)

## Mouse taxon ID 10090
## Zebrafish taxon ID 7955
## Tilapia taxon ID 8128

'cat unique_gtf2_genes.csv | while read GENE; do
datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 7955 "${GENE}" --filename "zebrafish/${GENE}".zip;
done'

'cat unique_gtf2_genes.csv | while read GENE; do
datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 10090 "${GENE}" --filename "mouse/${GENE}".zip;
done'

# Then to unzip all the .zip files in a folder, then the zip files can be deleted
'for f in *.zip; do unzip "$f" -d "${f%.zip}"; done'

## This one needs to recursively find all the paths
## Then load and extract the ortholog from the csv
## This output has duplicate entries (because of transcript variants), which are reconciled when combining with orthoDB

# For zeb
ncbi_orthologs_zeb <- generateNCBIorthologs(directory = "zebrafish/")
ncbi_orthologs_zeb <- Reduce(rbind, ncbi_orthologs_zeb)
saveRDS(ncbi_orthologs_zeb, file = "ncbi_orthologs_zeb.rds")

# For mouse
ncbi_orthologs_mouse <- generateNCBIorthologs(directory = "mouse/")
ncbi_orthologs_mouse <- Reduce(rbind, ncbi_orthologs_mouse)
saveRDS(ncbi_orthologs_mouse, file = "ncbi_orthologs_mouse.rds")

################################################################################################################
################ Query orthoDB database and ID orthologs from output ###########################################
################################################################################################################

## This one needs to run a bunch of system commands to query the orthoDB database
## Takes a Dbxref ID as input, and outputs a list of vectors (matching orthologs)
## These are better then the above, but miss some and should be supplemented
## This can be l/apply'd across a list of gene_ids

## This creates some unnecessary structure, but works
## Need to rework this function so I can run it for mouse
# then Unlist the unnecessary structure by using a non-recursive unlist

orthoDB_orthologs_zeb <- lapply(seq_along(gtf2$Dbxref_sub), function(x) generateOrthoDBorthologs(dbxrefs = gtf2$Dbxref_sub[x], gene_ids = gtf2$gene[x], species = "zebrafish"))
orthoDB_orthologs_zeb <- unlist(orthoDB_orthologs_zeb, recursive = F)
saveRDS(orthoDB_orthologs_zeb, file = "orthoDB_orthologs_zeb.rds")


orthoDB_orthologs_mouse <- readRDS(file = "orthoDB_orthologs_mouse")
# orthoDB_orthologs_mouse <- list()
# orthoDB_orthologs_mouse[length(orthoDB_orthologs_mouse):(length(orthoDB_orthologs_mouse)+1000)] <- lapply(length(orthoDB_orthologs_mouse):(length(orthoDB_orthologs_mouse)+1000), function(x) generateOrthoDBorthologs(dbxrefs = gtf2$Dbxref_sub[x], gene_ids = gtf2$gene[x], species = "mouse"))
orthoDB_orthologs_mouse[length(orthoDB_orthologs_mouse):(length(gtf2$Dbxref))] <- lapply(length(orthoDB_orthologs_mouse):(length(gtf2$Dbxref)), function(x) generateOrthoDBorthologs(dbxrefs = gtf2$Dbxref_sub[x], gene_ids = gtf2$gene[x], species = "mouse"))
# Save it out to reload and continue, as long as I don't unlist it recursively, I can keep going like this
saveRDS(orthoDB_orthologs_mouse, file = "orthoDB_orthologs_mouse")


orthoDB_orthologs_mouse <- unlist(orthoDB_orthologs_mouse, recursive = F)
saveRDS(orthoDB_orthologs_mouse, file = "orthoDB_orthologs_mouse")

# ################################################################################################################
# ################ Reconcile and combine the two database results together #######################################
# ################################################################################################################
# 
# ### OK, so the ncbi ortho  and orthoDB pipelines retrieve slightly different results (some missing orthos)
# ### There are X that ncbi retreives that orthoDB doesn't, and X vice versa - combined is more
# ## I can add the ncbi_orthologs to the orthoDB_orthologs list
# 
# ## If ncbi_orthos_zeb$gene_id_oreo is not in names(orthos_combined_zeb), then I want to add it
# 
# orthos_combined_zeb <- orthoDB_orthologs_zeb
# 
# for (i in 1:nrow(ncbi_orthologs_zeb)) {
#   if (!(ncbi_orthologs_zeb$gene_id_oreo[i] %in% names(orthos_combined_zeb))) {
#     orthos_combined_zeb[[ncbi_orthologs_zeb$gene_id_oreo[i]]] <- list(ncbi_orthologs_zeb$gene_symbol[i])
#   }
# }
# 
# orthos_combined_mouse <- orthoDB_orthologs_mouse
# 
# for (i in 1:nrow(ncbi_orthologs_mouse)) {
#   if (!(ncbi_orthologs_mouse$gene_id_oreo[i] %in% names(orthos_combined_mouse))) {
#     orthos_combined_mouse[[ncbi_orthologs_mouse$gene_id_oreo[i]]] <- list(ncbi_orthologs_mouse$gene_symbol[i])
#   }
# }
# 
# ## Save out the list
# 
# saveRDS(orthos_combined_zeb, file = "oreochromis_orthologs_zebrafish.rds")
# saveRDS(orthos_combined_mouse, file = "oreochromis_orthologs_mouse.rds")
# 
# orthos_combined_zeb <- readRDS(file = "oreochromis_orthologs_zebrafish.rds")
# orthos_combined_mouse <- readRDS(file = "oreochromis_orthologs_mouse.rds")


################################################################################################################
################ Reconcile and combine the two database results together #######################################
################################################################################################################

## I think the best actually would be to take the ncbi orthos first, then supplement with the orthoDBs, since they are less certain (more orthos)
## Also need to decided how to pick the orthoDB ortholog to use...

ncbi_orthologs_zeb <- readRDS("ncbi_orthologs_zeb.rds")
ncbi_orthologs_mouse <- readRDS("ncbi_orthologs_mouse.rds")

orthoDB_orthologs_zeb <- readRDS("orthoDB_orthologs_zeb.rds")
orthoDB_orthologs_mouse <- readRDS("orthoDB_orthologs_mouse.rds")

## Maybe I can start with the gtf file, and match first the ncbi results, then the orthoDB results

new_gtf <- gtf2

# Match by new_gtf$Dbxref_sub

new_gtf$zebrafish_ortholog_symbol <- ncbi_orthologs_zeb$gene_symbol[match(new_gtf$Dbxref_sub, ncbi_orthologs_zeb$gene_id_oreo)]
new_gtf$zebrafish_ortholog_id <- ncbi_orthologs_zeb$gene_id[match(new_gtf$Dbxref_sub, ncbi_orthologs_zeb$gene_id_oreo)]

new_gtf$mouse_ortholog_symbol <- ncbi_orthologs_mouse$gene_symbol[match(new_gtf$Dbxref_sub, ncbi_orthologs_mouse$gene_id_oreo)]
new_gtf$mouse_ortholog_id <- ncbi_orthologs_mouse$gene_id[match(new_gtf$Dbxref_sub, ncbi_orthologs_mouse$gene_id_oreo)]


## Add the orthoDB stuff??

new_gtf2 <- new_gtf

index <- is.na(new_gtf2$zebrafish_ortholog_id)
new_gtf2$zebrafish_ortholog_symbol[index] <- unlist(lapply(new_gtf2$Dbxref_sub[index], function(x) { 
  if (length(grep(x, names(orthoDB_orthologs_zeb))) %in% c(1)) {
    return(orthoDB_orthologs_zeb[[grep(x, names(orthoDB_orthologs_zeb))]][1]) # returning the first entry
  } else {
    return(NA)
  }
  })) # I need to pull the orthoDB results, based on the dbxref id number

index <- is.na(new_gtf$mouse_ortholog_id)
new_gtf2$mouse_ortholog_symbol[index] <- unlist(lapply(new_gtf2$Dbxref_sub[index], function(x) { 
  if (length(grep(x, names(orthoDB_orthologs_mouse))) %in% c(1)) {
    return(orthoDB_orthologs_mouse[[grep(x, names(orthoDB_orthologs_mouse))]][1]) # returning the first entry
  } else {
    return(NA)
  }
}))
## Save out the list

saveRDS(new_gtf, file = "gtf2_oreochromis_orthologs.rds")


## The problem is that for genes without ncbi based orthologs, orthoDB tends to return more 1:2, and 1:many orthologs (makes sense)
## 

# These are all the orthos
table(lengths(orthoDB_orthologs_zeb))

# With ncbi orthos
table(  lengths(    orthoDB_orthologs_zeb[names(orthoDB_orthologs_zeb) %in% as.character(new_gtf$Dbxref_sub)[!index]]         )     )

# Without ncbi orthos
table(  lengths(    orthoDB_orthologs_zeb[names(orthoDB_orthologs_zeb) %in% as.character(new_gtf$Dbxref_sub)[index]]         )     )






## There seem to be quite a few without zebrafish gene_symbols, but with an oreochromis gene symbol (897 to be exact)
## Can I sub these in (maybe in many cases they will match a real zeb gene symbol)
table(is.na(new_gtf$zebrafish_ortholog_symbol), is.na(new_gtf$gene.name)) # x by y

## 369 after adding orthoDB database
table(is.na(new_gtf2$zebrafish_ortholog_symbol), is.na(new_gtf2$gene.name)) # x by y








