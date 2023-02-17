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
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/orthos_gsea")

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]

## Ok, for each of these lists, need to output a file without header with the IDs (Dbxref)

index <- grepl("LOC", gtf2$Name)
gtf2$NCBI.GeneID <- NA
gtf2$NCBI.GeneID[index] <- substr(x$gene[index], start = 4, stop = 15)
gtf2$gene.name <- NA
gtf2$gene.name[!(index)] <- x$gene[!(index)]
return(x)




write.table(unique(substr(gtf2$Dbxref, start = 8, stop = 20)), file = "unique_gtf2_genes.csv", col.names = F, row.names = F, quote = F)


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
datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 10090 "${GENE}" --filename "mouse/${GENE}".zip;
datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 7955 "${GENE}" --filename "zebrafish/${GENE}".zip;
done'

# Then to unzip all the .zip files in a folder, then the zip files can be deleted
'for f in *.zip; do unzip "$f" -d "${f%.zip}"; done'

## This one needs to recursively find all the paths
## Then load and extract the ortholog from the csv
## This output has duplicate entries (because of transcript variants), which are reconciled when combining with orthoDB

# For zeb
ncbi_orthologs_zeb <- generateNCBIorthologs(directory = "zebrafish/")

# For mouse
ncbi_orthologs_mouse <- generateNCBIorthologs(directory = "mouse/")


################################################################################################################
################ Query orthoDB database and ID orthologs from output ###########################################
################################################################################################################

## This one needs to run a bunch of system commands to query the orthoDB database
## Takes a Dbxref ID as input, and outputs a list of vectors (matching orthologs)
## These are better then the above, but miss some and should be supplemented
## This can be l/apply'd across a list of gene_ids

## This creates some unnecessary structure, but works
## Need to rework this function so I can run it for mouse
tic()
orthoDB_orthologs_zeb <- lapply(unique(substr(gtf2$Dbxref, start = 8, stop = 20)), function(y) generateOrthoDBorthologs(dbxrefs = y, species = "zebrafish"))
toc()

orthoDB_orthologs_mouse <- lapply(unique(substr(gtf2$Dbxref, start = 8, stop = 20)), function(y) generateOrthoDBorthologs(dbxrefs = y, species = "mouse"))

# Unlist the unnecessary structure by using a non-recursive unlist
orthoDB_orthologs_zeb <- unlist(orthoDB_orthologs_zeb, recursive = F)
orthoDB_orthologs_mouse <- unlist(orthoDB_orthologs_mouse, recursive = F)


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

orthos_combined_zeb <- readRDS(file = "oreochromis_orthologs_zebrafish.rds")
orthos_combined_mouse <- readRDS(file = "oreochromis_orthologs_mouse.rds")

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
library(fgsea)

gene_list <- genes_new[[4]]$ps
names(gene_list) <- as.character(genes_new[[4]]$mouse_ortholog)
gene_list <- gene_list[order(gene_list, decreasing = T)]
gene_list <- gene_list[!(is.na(names(gene_list)))]

gene_list <- gene_list[!(duplicated(names(gene_list)))]




##### OK the below might work, but it requires ENTREZIDs for the orthologs, of which I only have gene symbols

### AHH I think the reason this doesn't work is that the gene lists are too short - ugh, that means quite a bit more processing for this orthology...

# One can also start from .rnk and .gmt files as in original GSEA:

rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package="fgsea")
gmt.file <- "~/Downloads/mh.all.v2022.1.Mm.symbols.gmt"
gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")

# Loading ranks:
ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)
str(ranks)

str(ranks)
##  Named num [1:12000] -63.3 -49.7 -43.6 -41.5 -33.3 ...
##  - attr(*, "names")= chr [1:12000] "170942" "109711" "18124" "12775" ...

# Loading pathways:

pathways <- gmtPathways(gmt.file)
str(head(pathways))


# And running fgsea:
  
fgseaRes <- fgsea(pathways, gene_list, minSize=15, maxSize=500)
head(fgseaRes)
















