library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)
library(fdrtool)


# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/sra_reads_nobackup/alleleFreqs")

## Read in the results from annotating SNPs
## Using a pvalue cutoff (of 0.0001) to ID significant SNPs
## Returns ~1k genes per comparison (only 216 for peaks)

genes1 <- readRDS(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "genes_ALL_SNPS.rds", sep = "_"))

genes2 <- lapply(genes1, function(x) x[x$distance_to_gene == 0,])
lapply(genes2, function(x) length(unique(x[x$distance_to_gene == 0, "location"])))

## Even fewer in exons at this point (don't need here)
exon_genes_2 <- readRDS(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "SNPSinExons.csv", sep = "_"))



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


## OK, now I need to convert these into human gene names
## This is 


ncbi_enseml <- read.csv("../../genome/mart_export_NCBI_Ensembl.txt", sep = "\t")
oreo_ortho <- read.csv("../../genome/mart_export_Oreochromis_Orthologs.txt", sep = "\t")
zeb_ortho <- read.csv("../../genome/mart_export_Zeb_orthos_107.txt", sep = "\t")
oreo_genes <- read.csv("../../genome/ncbi_genes_oreochromis_niloticus.tsv", sep = "\t")

lapply(genes, function(x) table(x$gene %in% oreo_ortho$Gene.name | x$gene %in% oreo_ortho))

## OK, goal is to have all of these match something
## I likely have the wrong version from ensembl?

lapply(genes, function(x) table(x$Dbxref %in% oreo_genes$NCBI.GeneID))


## OK, go from Dbxref (NCBI.ID) to Ensembl, then from there to mouse?

ensembl.ids <- ncbi_enseml$Gene.stable.ID[match(genes[[1]]$Dbxref, ncbi_enseml$NCBI.gene.ID)]

## This is from the newest versino of ensembl release - about half of the NCBI IDs aren't in this list (ugh)
new_ensembl <- read.csv("../../genome/mart_export_2023_01_09.txt")
lapply(genes, function(x) table(x$Dbxref %in% new_ensembl$NCBI.gene..formerly.Entrezgene..ID))



## Terminal commmand to use ncbi's command line tool
## This uses a program that I downloaded from here: https://www.ncbi.nlm.nih.gov/datasets/docs/v1/how-tos/genes/download-ortholog-dataset/

# datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 10090 --inputfile test_genes_3.csv --filename ort2.zip


## OK this for loop works a bit better, makes named folders for each Niloticus gene and puts the ortho inside

## Mouse 10090
## Zebrafish 7955
cat test_genes_3.csv | while read GENE; do
datasets download ortholog gene-id --exclude-gene --exclude-protein --exclude-rna --taxon-filter 7955 "${GENE}" --filename "${GENE}".zip;
done



write.csv(genes[[1]]$Dbxref, file = "~/Downloads/test_genes.csv", head = FALSE)


mouse.test <- oreo_ortho$Mouse.gene.stable.ID[match(ncbi_enseml$Gene.stable.ID[first_match], oreo_ortho$Gene.stable.ID)]


## Which ones don't match?

gene_test[!(gene_test$NCBI.GeneID[!(is.na(gene_test$NCBI.GeneID))] %in% oreo_genes$NCBI.GeneID),]
View(gene_test[!(gene_test$gene.name[!(is.na(gene_test$gene.name))] %in% oreo_genes$Symbol),])








