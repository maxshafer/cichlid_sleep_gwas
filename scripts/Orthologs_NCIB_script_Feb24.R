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
library(here)

library(clusterProfiler)
library(AnnotationDbi)
library(org.Dr.eg.db)

source("scripts/variants_functions.R")
# This should be run in the folder with the allele_frequncies

################################################################################################################
################ Generate lists of genes for finding orthologs #################################################
################################################################################################################

## Load gene table to convert to newest ids
ncbi_gene_table <- read.csv(here("genome","NCBI_gene-table_Oreochromis-niloticus_UMD_NMBU.tsv"), sep = "\t")

## Load in gtf
gtf <- read.gtf("genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf")
gtf <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]

## Load annotated snps
filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))
names(filter_snps) <- c("pc1_56-species", "pc2_56-species", "total_rest_56-species", "pc1_60-species", "pc2_60-species", "total_rest_60-species")

## Extract gene names
gene_ids <- lapply(filter_snps, function(x) {
  
  gof <- x$ANN_GENE
  gof <- unlist(strsplit(gof, "-"))
  gof <- unique(gof[gof %in% gtf$gene])
  
  
  ncbi_gof_id <- gtf$db_xref[match(gof, gtf$gene)]
  ncbi_gof_id <- substr(ncbi_gof_id, 8, 30)
  
  ncbi_gof <- ncbi_gene_table$Symbol[match(ncbi_gof_id, ncbi_gene_table$NCBI.GeneID)]
  
  return(ncbi_gof)
  
})

dbxrefs <- lapply(filter_snps, function(x) {
  
  gof <- x$ANN_GENE
  gof <- unlist(strsplit(gof, "-"))
  gof <- unique(gof[gof %in% gtf$gene])
  
  
  ncbi_gof_id <- gtf$db_xref[match(gof, gtf$gene)]
  ncbi_gof_id <- substr(ncbi_gof_id, 8, 30)
  
  # ncbi_gof <- ncbi_gene_table$Symbol[match(ncbi_gof_id, ncbi_gene_table$NCBI.GeneID)]
  
  return(ncbi_gof_id)
  
})

## Ok, for each of these lists, need to output a file without header with the IDs (Dbxref)

write.table(dbxrefs[[1]], file = here("orthos", "perchr_56sp_pc1_genes.csv"), col.names = F, row.names = F, quote = F)
write.table(dbxrefs[[2]], file = here("orthos", "perchr_56sp_pc2_genes.csv"), col.names = F, row.names = F, quote = F)
write.table(dbxrefs[[3]], file = here("orthos", "perchr_56sp_tr_genes.csv"), col.names = F, row.names = F, quote = F)

write.table(dbxrefs[[4]], file = here("orthos", "perchr_60sp_pc1_genes.csv"), col.names = F, row.names = F, quote = F)
write.table(dbxrefs[[5]], file = here("orthos", "perchr_60sp_pc2_genes.csv"), col.names = F, row.names = F, quote = F)
write.table(dbxrefs[[6]], file = here("orthos", "perchr_60sp_tr_genes.csv"), col.names = F, row.names = F, quote = F)

################################################################################################################
################ Run terminal commands to download NCBI orthologs ##############################################
################################################################################################################

## Need to run these system commands for all of the lists (within the directories)
## Makes named folders for each Niloticus gene and puts the ortho inside

## Mouse taxon ID 10090
## Zebrafish taxon ID 7955
## Tilapia taxon ID 8128

## Can just paste a bunch of commands?

human_pc1 <- paste('./datasets download gene gene-id ', dbxrefs[[4]],' --ortholog human --filename ', "'",file.path("pc1", "human", dbxrefs[[4]], fsep = "/"), ".zip'", sep = "")
human_pc2 <- paste('./datasets download gene gene-id ', dbxrefs[[5]],' --ortholog human --filename ', "'",file.path("pc2", "human", dbxrefs[[5]], fsep = "/"), ".zip'", sep = "")
human_tr <- paste('./datasets download gene gene-id ', dbxrefs[[6]],' --ortholog human --filename ', "'",file.path("tr", "human", dbxrefs[[6]], fsep = "/"), ".zip'", sep = "")

mouse_pc1 <- paste('./datasets download gene gene-id ', dbxrefs[[4]],' --ortholog 10090 --filename ', "'",file.path("pc1", "mouse", dbxrefs[[4]], fsep = "/"), ".zip'", sep = "")
mouse_pc2 <- paste('./datasets download gene gene-id ', dbxrefs[[5]],' --ortholog 10090 --filename ', "'",file.path("pc2", "mouse", dbxrefs[[5]], fsep = "/"), ".zip'", sep = "")
mouse_tr <- paste('./datasets download gene gene-id ', dbxrefs[[6]],' --ortholog 10090 --filename ', "'",file.path("tr", "mouse", dbxrefs[[6]], fsep = "/"), ".zip'", sep = "")

zeb_pc1 <- paste('./datasets download gene gene-id ', dbxrefs[[4]],' --ortholog zebrafish --filename ', "'",file.path("pc1", "zebrafish", dbxrefs[[4]], fsep = "/"), ".zip'", sep = "")
zeb_pc2 <- paste('./datasets download gene gene-id ', dbxrefs[[5]],' --ortholog zebrafish --filename ', "'",file.path("pc2", "zebrafish", dbxrefs[[5]], fsep = "/"), ".zip'", sep = "")
zeb_tr <- paste('./datasets download gene gene-id ', dbxrefs[[6]],' --ortholog zebrafish --filename ', "'",file.path("tr", "zebrafish", dbxrefs[[6]], fsep = "/"), ".zip'", sep = "")

write.table(c(human_pc1, human_pc2, human_tr), file = here("orthos", "human_script.csv"), col.names = F, row.names = F, quote = F)
write.table(c(zeb_pc1, zeb_pc2, zeb_tr), file = here("orthos", "zebrafish_script.csv"), col.names = F, row.names = F, quote = F)
write.table(c(mouse_pc1, mouse_pc2, mouse_tr), file = here("orthos", "mouse_script.csv"), col.names = F, row.names = F, quote = F)

## Then just run it as a sh in the terminal

# Then to unzip all the .zip files in a folder, then the zip files can be deleted
'for f in *.zip; do unzip "$f" -d "${f%.zip}"; done'

## This one needs to recursively find all the paths
## Then load and extract the ortholog from the csv
## This output has duplicate entries (because of transcript variants), which are reconciled when combining with orthoDB

ncbi_orthologs <- list(list(), list(), list())
names(ncbi_orthologs) <- c("pc1", "pc2", "tr")
# For PC1
directories_pc1 <- c("orthos/pc1/human/", "orthos/pc1/mouse/", "orthos/pc1/zebrafish/")
ncbi_orthologs[[1]] <- lapply(directories_pc1, function(x) generateNCBIorthologsJSONL(directory = x))
names(ncbi_orthologs[[1]]) <- c("human", "mouse", "zebrafish")
ncbi_orthologs[[1]] <- lapply(ncbi_orthologs[[1]], function(x) unlist(x))

# For mouse
directories_pc2 <- c("orthos/pc2/human/", "orthos/pc2/mouse/", "orthos/pc2/zebrafish/")
ncbi_orthologs[[2]] <- lapply(directories_pc2, function(x) generateNCBIorthologsJSONL(directory = x))
names(ncbi_orthologs[[2]]) <- c("human", "mouse", "zebrafish")
ncbi_orthologs[[2]] <- lapply(ncbi_orthologs[[2]], function(x) unlist(x))

# For zeb
directories_tr <- c("orthos/tr/human/", "orthos/tr/mouse/", "orthos/tr/zebrafish/")
ncbi_orthologs[[3]] <- lapply(directories_tr, function(x) generateNCBIorthologsJSONL(directory = x))
names(ncbi_orthologs[[3]]) <- c("human", "mouse", "zebrafish")
ncbi_orthologs[[3]] <- lapply(ncbi_orthologs[[3]], function(x) unlist(x))


################################################################################################################
################ Query ensembl using matching gene names #######################################################
################################################################################################################

# Load in ensembl indices for IDing orthologs
# These were downloaded in Feb 2024, using Ensembl 111
ensembl_btw <- list(list(), list(), list())
names(ensembl_btw) <- c("human", "mouse", "zebrafish")
ensembl_btw$human <- read.csv("ensembl_data/ensembl-O_niloticus_UMD_NMBU-Human.txt", sep = "\t")
ensembl_btw$mouse <- read.csv("ensembl_data/ensembl-O_niloticus_UMD_NMBU-Mouse.txt", sep = "\t")
ensembl_btw$zebrafish <- read.csv("ensembl_data/ensembl-O_niloticus_UMD_NMBU-Zebrafish.txt", sep = "\t")

ensembl_ncbi_id <- read.csv("ensembl_data/O_niloticus_UMD_NMBU-ensembl_ncbi-id.txt")
ensembl_ncbi_ac <- read.csv("ensembl_data/O_niloticus_UMD_NMBU-ensembl_ncbi-accession.txt")

ensembl_orthologs <- lapply(c(2:4), function(x) {
  gof <- unique(filter_snps[[x]]$ANN_GENE)
  gof <- unlist(strsplit(gof, "-"))
  idof <- gtf$db_xref[match(gof, gtf$gene)]
  idof <- str_sub(idof, start = 8, end = 30)
  
  # Find Oreochromis ENSEMBL ID by matching ncbi gene ID in both names and accessions
  ensembl_names <- unique(c(ensembl_ncbi_ac$Gene.name[match(gof, ensembl_ncbi_ac$NCBI.gene..formerly.Entrezgene..accession)], 
                            ensembl_ncbi_id$Gene.name[match(gof, ensembl_ncbi_id$NCBI.gene..formerly.Entrezgene..ID)]))
  
  # Find Oreochromis ENSEMBL ID by matching gene name in both names and accessions
  ensembl_ids <- unique(c(ensembl_ncbi_ac$Gene.stable.ID[match(idof, ensembl_ncbi_ac$NCBI.gene..formerly.Entrezgene..accession)],
                          ensembl_ncbi_id$Gene.stable.ID[match(idof, ensembl_ncbi_id$NCBI.gene..formerly.Entrezgene..ID)]))
  
  ## Create nested output structure
  out <- list(list(list(), list()), list(list(), list()), list(list(), list()))
  names(out) <- c("human", "mouse", "zebrafish")

  
  ## use Oreochromis ENSEMBL IDs and names to find matching ENSEMBL IDs and names for human
  out$human[[1]] <- unique(c(ensembl_btw$human$Human.gene.name[match(ensembl_ids, ensembl_btw$human$Gene.stable.ID)],
                             ensembl_btw$human$Human.gene.name[match(ensembl_names, ensembl_btw$human$Gene.name)]))
  out$human[[2]] <- unique(c(ensembl_btw$human$Human.gene.stable.ID[match(ensembl_ids, ensembl_btw$human$Gene.stable.ID)],
                             ensembl_btw$human$Human.gene.stable.ID[match(ensembl_names, ensembl_btw$human$Gene.name)]))
  
  ## use Oreochromis ENSEMBL IDs and names to find matching ENSEMBL IDs and names for mouse
  out$mouse[[1]] <- unique(c(ensembl_btw$mouse$Mouse.gene.name[match(ensembl_ids, ensembl_btw$mouse$Gene.stable.ID)],
                             ensembl_btw$mouse$Mouse.gene.name[match(ensembl_names, ensembl_btw$mouse$Gene.name)]))
  out$mouse[[2]] <- unique(c(ensembl_btw$mouse$Mouse.gene.stable.ID[match(ensembl_ids, ensembl_btw$mouse$Gene.stable.ID)],
                             ensembl_btw$mouse$Mouse.gene.stable.ID[match(ensembl_names, ensembl_btw$mouse$Gene.name)]))
  
  ## use Oreochromis ENSEMBL IDs and names to find matching ENSEMBL IDs and names for mouse
  out$zebrafish[[1]] <- unique(c(ensembl_btw$zebrafish$Zebrafish.gene.name[match(ensembl_ids, ensembl_btw$zebrafish$Gene.stable.ID)],
                                 ensembl_btw$zebrafish$Zebrafish.gene.name[match(ensembl_names, ensembl_btw$zebrafish$Gene.name)]))
  out$zebrafish[[1]] <- unique(c(ensembl_btw$zebrafish$Zebrafish.gene.stable.ID[match(ensembl_ids, ensembl_btw$zebrafish$Gene.stable.ID)],
                                 ensembl_btw$zebrafish$Zebrafish.gene.stable.ID[match(ensembl_names, ensembl_btw$zebrafish$Gene.name)]))
  
  return(out)
})
  

################################################################################################################
################ Combine and save out ##########################################################################
################################################################################################################

## First combine gene names across ensembl and ncbi
## Then use this to map back from gene names to ensembl IDs and combine with those

genes_of_interest_names <- lapply(seq_along(ncbi_orthologs), function(x) {
  out <- lapply(seq_along(ncbi_orthologs[[x]]), function(y) {
    out2 <- unique(c(ncbi_orthologs[[x]][[y]], ensembl_orthologs[[x]][[y]][[1]]))
    return(out2)
  })
  
  names(out) <- c("human", "mouse", "zebrafish")
  return(out)
})

names(genes_of_interest_names) <- c("pc1", "pc2", "tr")

str(genes_of_interest_names, max.level = 2)
# List of 3
# $ pc1:List of 3
# ..$ human    : chr [1:601] "TSHR" "DMRT2" "TCF21" "PCNX2" ...
# ..$ mouse    : chr [1:597] "Tshr" "Dmrt2" "Tcf21" "Pcnx2" ...
# ..$ zebrafish: chr [1:1132] "tshr" "dmrt2a" "tnfb" "tcf21" ...
# $ pc2:List of 3
# ..$ human    : chr [1:585] "METTL24" "OLIG3" "CFAP299" "DAZAP1" ...
# ..$ mouse    : chr [1:583] "Mettl24" "Olig3" "Cfap299" "Dazap1" ...
# ..$ zebrafish: chr [1:1121] "mettl24" "olig3" "cfap299" "dazap1" ...
# $ tr :List of 3
# ..$ human    : chr [1:553] "SASH1" "SRBD1" "NT5DC1" "TENM4" ...
# ..$ mouse    : chr [1:553] "Sash1" "Srbd1" "Nt5dc1" "Tenm4" ...
# ..$ zebrafish: chr [1:1048] "sash1a" "mpzl3" "srbd1" "nt5dc1" ...


## Now I need to use these IDs to match ENSEMBL IDs (particularly for those found with ncbi), and combine with the list of ENSEMBL IDs

genes_of_interest_ids <- lapply(seq_along(ncbi_orthologs), function(x) {
  
  out <- lapply(seq_along(ncbi_orthologs[[x]]), function(y) {
        out2 <- unique(ensembl_btw[[x]][match(ncbi_orthologs[[x]][[y]], ensembl_btw[[x]][,4]),3])
        out2 <- unique(c(out2, ensembl_orthologs[[x]][[y]][[2]]))
    return(out2)
  })
  
  names(out) <- c("human", "mouse", "zebrafish")
  return(out)
})

names(genes_of_interest_ids) <- c("pc1", "pc2", "tr")

str(genes_of_interest_ids, max.level = 2)
# List of 3
# $ pc1:List of 3
# ..$ human    : chr [1:575] "ENSG00000165409" "ENSG00000173253" "ENSG00000118526" "ENSG00000135749" ...
# ..$ mouse    : chr [1:470] NA "ENSMUSG00000043822" "ENSMUSG00000010755" "ENSMUSG00000070348" ...
# ..$ zebrafish:List of 1
# $ pc2:List of 3
# ..$ human    : chr [1:456] NA "ENSG00000162105" "" "ENSG00000110092" ...
# ..$ mouse    : chr [1:559] "ENSMUSG00000045555" "ENSMUSG00000045591" "ENSMUSG00000057816" "ENSMUSG00000069565" ...
# ..$ zebrafish:List of 1
# $ tr :List of 3
# ..$ human    : chr [1:429] NA "ENSDARG00000092855" "ENSDARG00000088100" "ENSDARG00000041372" ...
# ..$ mouse    : chr [1:414] NA "ENSMUSG00000053046" "ENSMUSG00000027253" "ENSMUSG00000032340" ...
# ..$ zebrafish:List of 471

genes_of_interest <- list(genes_of_interest_names, genes_of_interest_ids)
names(genes_of_interest) <- c("gene_names", "gene_ids")
saveRDS(genes_of_interest, file = "orthos/orthologs_perchr_60sp_1e-05_percentile.rds")

## These lists then get put into Gene Analytics and download, or David
write_clip(genes_of_interest_names$pc1$human)
write_clip(genes_of_interest_names$pc2$human)
write_clip(genes_of_interest_names$tr$human)


## These lits can be put into the MGI pheno enrichment tool
write_clip(genes_of_interest$gene_ids$pc1$human) #M
write_clip(genes_of_interest$gene_ids$pc2$human) #M
write_clip(genes_of_interest$gene_ids$tr$human) #M

write_clip(genes_of_interest$gene_ids$pc1$mouse)
write_clip(genes_of_interest$gene_ids$pc2$mouse)
write_clip(genes_of_interest$gene_ids$tr$mouse) #M

# saved in /orthos/modPhEA/


## The below doesn't work anymore, there is some difference in the query to orthodb that has changed, and I'm not sure how to get the same info again
## Either way, combining ncbi orthology and ensembl orthology (by matching gene names) gives me a ton of orthologs (60%, 57%, and 57%)


# ################################################################################################################
# ################ Query orthoDB database and ID orthologs from output ###########################################
# ################################################################################################################
# 
# ## This one needs to run a bunch of system commands to query the orthoDB database
# ## Takes a Dbxref ID as input, and outputs a list of vectors (matching orthologs)
# ## These are better then the above, but miss some and should be supplemented
# ## This can be l/apply'd across a list of gene_ids
# 
# ## This creates some unnecessary structure, but works
# ## Need to rework this function so I can run it for mouse
# 
# orthoDB_orthologs_zeb <- lapply(genes, function(x) lapply(x[,"Dbxref"], function(y) generateOrthoDBorthologs(dbxrefs = y, species = "zebrafish")))
# names(orthoDB_orthologs_zeb) <- c("dn55", "dn58", "dn60", "tr", "peak")
# 
# orthoDB_orthologs_mouse <- lapply(genes, function(x) lapply(x[,"Dbxref"], function(y) generateOrthoDBorthologs(dbxrefs = y, species = "mouse")))
# names(orthoDB_orthologs_mouse) <- c("dn55", "dn58", "dn60", "tr", "peak")
# 
# # Unlist the unnecessary structure by using a non-recursive unlist
# orthoDB_orthologs_zeb <- lapply(orthoDB_orthologs_zeb, function(x) unlist(x, recursive = F))
# orthoDB_orthologs_mouse <- lapply(orthoDB_orthologs_mouse, function(x) unlist(x, recursive = F))
# 
# 
# ################################################################################################################
# ################ Reconcile and combine the two database results together #######################################
# ################################################################################################################
# 
# ### OK, so the ncbi ortho  and orthoDB pipelines retrieve slightly different results (some missing orthos)
# ### There are X that ncbi retreives that orthoDB doesn't, and X vice versa - combined is more
# ## I can add the ncbi_orthologs to the orthoDB_orthologs list
# 
# orthos_combined_zeb <- orthoDB_orthologs_zeb
# for (j in 1:length(orthos_combined_zeb)) {
#   for (i in 1:nrow(ncbi_orthologs_zeb[[j]])) {
#     if (!(ncbi_orthologs_zeb[[j]]$gene_id_oreo[i] %in% names(orthos_combined_zeb[[j]]))) {
#       orthos_combined_zeb[[j]][[ncbi_orthologs_zeb[[j]]$gene_id_oreo[i]]] <- list(ncbi_orthologs_zeb[[j]]$gene_symbol[i])
#     }
#   }
# }
# 
# orthos_combined_mouse <- orthoDB_orthologs_mouse
# for (j in 1:length(orthos_combined_mouse)) {
#   for (i in 1:nrow(ncbi_orthologs_mouse[[j]])) {
#     if (!(ncbi_orthologs_mouse[[j]]$gene_id_oreo[i] %in% names(orthos_combined_mouse[[j]]))) {
#       orthos_combined_mouse[[j]][[ncbi_orthologs_mouse[[j]]$gene_id_oreo[i]]] <- list(ncbi_orthologs_mouse[[j]]$gene_symbol[i])
#     }
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
# 
# ## Add them back to the genes dataframes
# ## This selects the first entry in the ortho database as the ortho
# 
# genes_new <- list()
# for (i in seq_along(genes)) {
#   
#   zeb_ortholog <- unlist(lapply(orthos_combined_zeb[[i]], function(x) {
#     if (length(x) > 0) {
#       return(x[[1]][[1]][1])
#     } else {
#       return(NA)
#     }}))
#   zeb_ortholog <- data.frame(Dbxref = names(zeb_ortholog), zeb_ortholog = zeb_ortholog)
#   
#   mouse_ortholog <- unlist(lapply(orthos_combined_mouse[[i]], function(x) {
#     if (length(x) > 0) {
#       return(x[[1]][[1]][1])
#     } else {
#       return(NA)
#     }}))
#   mouse_ortholog <- data.frame(Dbxref = names(mouse_ortholog), mouse_ortholog = mouse_ortholog)
#   genes_new[[i]] <- genes[[i]]
#   genes_new[[i]]$zeb_ortholog <- zeb_ortholog$zeb_ortholog[match(genes_new[[i]]$Dbxref, zeb_ortholog$Dbxref)]
#   genes_new[[i]]$mouse_ortholog <- mouse_ortholog$mouse_ortholog[match(genes_new[[i]]$Dbxref, mouse_ortholog$Dbxref)]
# }
# 
# 
# ### Run some GSEA?
# library(org.Mm.eg.db)
# library(fgsea)
# 
# gene_list <- genes_new[[4]]$ps
# names(gene_list) <- as.character(genes_new[[4]]$mouse_ortholog)
# gene_list <- gene_list[order(gene_list, decreasing = T)]
# gene_list <- gene_list[!(is.na(names(gene_list)))]
# 
# gene_list <- gene_list[!(duplicated(names(gene_list)))]
# 
# 
# 
# 
# ##### OK the below might work, but it requires ENTREZIDs for the orthologs, of which I only have gene symbols
# 
# ### AHH I think the reason this doesn't work is that the gene lists are too short - ugh, that means quite a bit more processing for this orthology...
# 
# # One can also start from .rnk and .gmt files as in original GSEA:
# 
# rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package="fgsea")
# gmt.file <- "~/Downloads/mh.all.v2022.1.Mm.symbols.gmt"
# gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")
# 
# # Loading ranks:
# ranks <- read.table(rnk.file,
#                     header=TRUE, colClasses = c("character", "numeric"))
# ranks <- setNames(ranks$t, ranks$ID)
# str(ranks)
# 
# str(ranks)
# ##  Named num [1:12000] -63.3 -49.7 -43.6 -41.5 -33.3 ...
# ##  - attr(*, "names")= chr [1:12000] "170942" "109711" "18124" "12775" ...
# 
# # Loading pathways:
# 
# pathways <- gmtPathways(gmt.file)
# str(head(pathways))
# 
# 
# # And running fgsea:
#   
# fgseaRes <- fgsea(pathways, gene_list, minSize=15, maxSize=500)
# head(fgseaRes)
















