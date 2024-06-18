## I want this all to be run across lists, and be hidden inside of functions
## Mostly so I can run this across the permuted SNP lists
library(Rgb)
library(here)
library(stringr)

library(ggrepel)
library(data.table)
library(dplyr)
library(readr)
library(patchwork)
library(fdrtool)
library(rjson)
library(clipr)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Dr.eg.db)

source("scripts/variants_functions.R")

## Load in things I need

## Load in gtf
gtf <- read.gtf("genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf")
gtf <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]

# ## Load in permutated ids and filter_snps
# 
# ## Load annotated snps
# filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))
# names(filter_snps) <- c("pc1_56-species", "pc2_56-species", "total_rest_56-species", "pc1_60-species", "pc2_60-species", "total_rest_60-species")
# 
# permuted_snps <- readRDS(here("sra_reads_nobackup/combined_ann/Permute_SNPs_perchr_1000_permutations_800_SNPs"))
# 
# permuted_genes <- lapply(permuted_snps, function(x) x$ANN_GENE)

## OK, so I have as input dataframes with gene ids
## To find ncbi orthologs I need a function that does:
## 1) Find dbxref ids from the ncbi_gene_table
## 2) Outputs an sh file with commands to download ncbi orthologs for each id
## 3) Runs that command, then unzips the results
## 4) Loads all the files and makes a list of orthologs (potentially named with the oreo dbxref)
## I should be able to do this in a function, using system()

idDBXREFs <- function(gene_list = gene_list, gtf = gtf) {
  gof <- gene_list$ANN_GENE
  gof <- unlist(strsplit(gof, "-"))
  gof <- unique(gof[gof %in% gtf$gene])
  
  ncbi_gof_id <- gtf$db_xref[match(gof, gtf$gene)]
  ncbi_gof_id <- substr(ncbi_gof_id, 8, 30)
  
  return(ncbi_gof_id)
}

## Should this work per dbxref? Or per list of dbxrefs? I think the former would be cleaner?
## OK this doesn't seem to work for one by one, I should run the steps separately, but for a list
ncbiOrthologsID <- function(dbxref_id = dbxref, species = "human", directory = here("orthos", "permutations")) {
  
  # Run the ncbi command
  system(paste('./orthos/datasets download gene gene-id ', dbxref_id,' --ortholog human --filename ', "'", file.path(directory, dbxref_id, fsep = "/"), ".zip'", sep = ""))
  # Unzip the result
  #system(paste('for f in ', directory, '/*.zip; do unzip "$f" -d "${f%.zip}"; done', sep = ""))
  system(paste("unzip -o -d ", directory, "/", dbxref_id, "/ ", directory, "/", dbxref_id, ".zip", sep = ""))
  # Load and extract the ortholog(s)
  
  
  # Remove the folder (should be temporary)
  # unlink(paste(directory, "ncbi_dataset", sep = "/"), recursive = T)
  unlink(paste(directory, "/", dbxref_id, ".zip", sep = ""))
}



## Run for all unique gene names in the gtf

unique_dbxrefs <- substr(unique(gtf$db_xref), 8, 30)

## Did this actually run, and finish? I know that it was hanging on something when I came in
## There is definitely an issue with the fact that I am not randomly sampling hte human genome, therefore all teh samples have some enrichments
## My data also suffers from this, b/c of the human ortholog problem. The background should be the human orthos for oreochromis genes
lapply(unique_dbxrefs, function(x) ncbiOrthologsID(dbxref_id = x, species = "human", directory = here("orthos", "permutations")))

out <- generateNCBIorthologsJSONL(directory = here("orthos", "permutations"))
out <- unlist(out)
names(out) <- list.dirs(here("orthos", "permutations"), recursive = F, full.names = F)


################################################################################################################
################ Use the list to ID orthologs for each permutation #############################################
################################################################################################################

ortholog_df <- data.frame(oreochromis_gene = gtf$gene, ncbi_human_ortholog = out[match(gtf$gene, gtf$gene[match(paste("GeneID:", names(out), sep = ""), gtf$db_xref)])], dbxref = str_sub(gtf$db_xref[match(gtf$gene, gtf$gene)], start = 8, end = 30))

# permuted_genes_df <- lapply(permuted_genes, function(x) {
#   gof <- unlist(strsplit(x, "-"))
#   output <- out[match(gof, gtf$gene[match(paste("GeneID:", names(out), sep = ""), gtf$db_xref)])]
#   df <- data.frame(oreochromis_gene = gof, ncbi_human_ortholog = output, dbxref = str_sub(gtf$db_xref[match(gof, gtf$gene)], start = 8, end = 30))
# })


################################################################################################################
################ Query ensembl using matching gene names #######################################################
################################################################################################################

# Load in ensembl indices for IDing orthologs
# These were downloaded in Feb 2024, using Ensembl 111
# ensembl_btw <- list(list(), list(), list())
# names(ensembl_btw) <- c("human", "mouse", "zebrafish")
ensembl_btw <- read.csv("ensembl_data/ensembl-O_niloticus_UMD_NMBU-Human.txt", sep = "\t")
# ensembl_btw$mouse <- read.csv("ensembl_data/ensembl-O_niloticus_UMD_NMBU-Mouse.txt", sep = "\t")
# ensembl_btw$zebrafish <- read.csv("ensembl_data/ensembl-O_niloticus_UMD_NMBU-Zebrafish.txt", sep = "\t")

ensembl_ncbi_id <- read.csv("ensembl_data/O_niloticus_UMD_NMBU-ensembl_ncbi-id.txt")
ensembl_ncbi_ac <- read.csv("ensembl_data/O_niloticus_UMD_NMBU-ensembl_ncbi-accession.txt")




permuted_genes_df <- lapply(ortholog_df, function(x) {
  
  # Find Oreochromis ENSEMBL ID by matching ncbi gene ID in both names and accessions
  x$ensembl_names_1 <- ensembl_ncbi_ac$Gene.name[match(x$oreochromis_gene, ensembl_ncbi_ac$NCBI.gene..formerly.Entrezgene..accession)]
  x$ensembl_names_2 <- ensembl_ncbi_id$Gene.name[match(x$dbxref, ensembl_ncbi_id$NCBI.gene..formerly.Entrezgene..ID)]
  
  x$ensembl_names <- apply(x, 1, function(i) {
    if(is.na(i[4]) | is.na(i[5])) {
      if(is.na(i[4]) & is.na(i[5])) {
        return("")
      } else {
        out <- unlist(c(i[4], i[5]))
        out <- out[!(is.na(out))]
        out <- unique(out)
        return(out)
      }
    } else {
      return(i[4])
    }
  })
  
  
  # Find Oreochromis ENSEMBL ID by matching gene name in both names and accessions
  x$ensembl_ids_1 <- ensembl_ncbi_ac$Gene.stable.ID[match(x$oreochromis_gene, ensembl_ncbi_ac$NCBI.gene..formerly.Entrezgene..accession)]
  x$ensembl_ids_2 <- ensembl_ncbi_id$Gene.stable.ID[match(x$dbxref, ensembl_ncbi_id$NCBI.gene..formerly.Entrezgene..ID)]
  
  x$ensembl_ids <- apply(x, 1, function(i) {
    if(is.na(i[7]) | is.na(i[8])) {
      if(is.na(i[7]) & is.na(i[8])) {
        return("")
      } else {
        out <- unlist(c(i[7], i[8]))
        out <- out[!(is.na(out))]
        out <- unique(out)
        return(out)
      }
    } else {
      return(i[7])
    }
  })
  
  ## use Oreochromis ENSEMBL IDs and names to find matching ENSEMBL IDs and names for human
  
  x$ensembl_names_human_1 <- ensembl_btw$Human.gene.name[match(x$ensembl_ids, ensembl_btw$Gene.stable.ID)]
  x$ensembl_names_human_2 <- ensembl_btw$Human.gene.name[match(x$ensembl_names, ensembl_btw$Gene.name)]
  
  x$ensembl_names_human <- apply(x, 1, function(i) {
    if(is.na(i[10]) | is.na(i[11])) {
      if(is.na(i[10]) & is.na(i[11])) {
        return("")
      } else {
        out <- unlist(c(i[10], i[11]))
        out <- out[!(is.na(out))]
        out <- unique(out)
        return(out)
      }
    } else {
      return(i[11])
    }
  })
  
  x$human_ortholog <- apply(x, 1, function(i) {
    if(is.na(i[2]) | is.na(i[12])) {
      if(is.na(i[2]) & is.na(i[12])) {
        return("")
      } else {
        out <- unlist(c(i[2], i[12]))
        out <- out[!(is.na(out))]
        out <- unique(out)
        return(out)
      }
    } else {
      return(i[2])
    }
  })
  
  # # This is for pulling ensembl IDs, which might be useful (and can also typically be queried in GO analyses)
  # x$ensembl_ids_human_1 <- ensembl_btw$Human.gene.stable.ID[match(x$ensembl_ids, ensembl_btw$human$Gene.stable.ID)]
  # x$ensembl_ids_human_2 <- ensembl_btw$Human.gene.stable.ID[match(x$ensembl_names, ensembl_btw$human$Gene.name)]

  return(x[,c(1,3,13)])
})


################################################################################################################
################ Combine and save out ##########################################################################
################################################################################################################

saveRDS(ortholog_df, file = "GTF_human_orthologs_all.rds")
saveRDS(permuted_genes_df, file = "Permute_SNPs_perchr_1000_permutations_800_SNPs_orthologs.rds")




