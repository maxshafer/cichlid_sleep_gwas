library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)
library(fdrtool)
library(DescTools)

# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas")

source('/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R')

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]
gtf3 <- gtf[gtf$feature == "mRNA",]

## Load in perchr.list, and summarise it by 10k windows
## I'm pretty sure this is PGLS-GWAS pvalue (can do reg pvalue also)

perchr.list <- readRDS("sra_reads_nobackup/per_chromosome_list.rds")

perchr.list.10k <- lapply(perchr.list, function(x) x %>% group_by(chr, ints = cut_width(coord, width = 10000, boundary = 0)) %>% summarise(start = min(coord), end = max(coord), mean_pvalue = mean(pvalue), pvalue = min(pvalue)) )

rm(perchr.list)

saveRDS(perchr.list.10k, file = "sra_reads_nobackup/per_chromosome_list_10k.rds")

## Read back in 
perchr.list.10k <- readRDS(file = "sra_reads_nobackup/per_chromosome_list_10k.rds")

## I want to run the above per chromosome (which is faster)

## First make lists of gtf dataframes for each chromosome
gtf2_chrs <- unique(gtf2$seqname)[1:23]
gtf_df <- data.frame(gtf_chrs = gtf2_chrs, vcf_chrs = as.character(substr(gtf2_chrs, 1, 9)))
list_gtf2_sub <- lapply(gtf2_chrs, function(x) gtf2[gtf2$seqname %in% x,])
names(list_gtf2_sub) <- gtf2_chrs

list_gtf3_sub <- lapply(gtf2_chrs, function(x) gtf3[gtf3$seqname %in% x,])
names(list_gtf3_sub) <- gtf2_chrs

## This now runs getGenesEns using chromosome specific data frames, which increases speed by a lot
## I'm also only running this once per 10kb window, which reduces the # of comparisons
## Still retains many of the top genes, so I think this is a good alternative and should work for GSEA (once I ID orthologs)

output_mean <- lapply(perchr.list.10k, function(perchr) {
  perchr <- perchr[,c(1,2,3,4,5)] # 6 is min pvalue, 5 is mean
  per_chrs <- unique(perchr$chr)
  perchr <- lapply(per_chrs, function(l) perchr[perchr$chr %in% l,])
  names(perchr) <- per_chrs
  out <- lapply(per_chrs[2:24], function(k) {
    put <- getGenesEns(results.table = perchr[[grep(k, names(perchr))]], gene_table = list_gtf2_sub[[gtf_df$gtf_chrs[match(k, gtf_df$vcf_chrs)]]], gene_table_2 = list_gtf3_sub[[gtf_df$gtf_chrs[match(k, gtf_df$vcf_chrs)]]], window = 25000, windowed = TRUE, combined = FALSE)
  })
  out <- Reduce(rbind, out)
})

names(output_mean) <- c("dn55", "dn58", "dn60", "tr", "peak")

output_min <- lapply(perchr.list.10k, function(perchr) {
  perchr <- perchr[,c(1,2,3,4,6)] # 6 is min pvalue, 5 is mean
  per_chrs <- unique(perchr$chr)
  perchr <- lapply(per_chrs, function(l) perchr[perchr$chr %in% l,])
  names(perchr) <- per_chrs
  out <- lapply(per_chrs[2:24], function(k) {
    put <- getGenesEns(results.table = perchr[[grep(k, names(perchr))]], gene_table = list_gtf2_sub[[gtf_df$gtf_chrs[match(k, gtf_df$vcf_chrs)]]], gene_table_2 = list_gtf3_sub[[gtf_df$gtf_chrs[match(k, gtf_df$vcf_chrs)]]], window = 25000, windowed = TRUE, combined = FALSE)
  })
  out <- Reduce(rbind, out)
})

names(output_min) <- c("dn55", "dn58", "dn60", "tr", "peak")

## OK, this seems to work quite well, and returns genes that look similar (at least for total rest)

saveRDS(output_mean, file = "/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_10kb_windows_meanPvalue.rds")
saveRDS(output_min, file = "/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_10kb_windows_minPvalue.rds")


################################################################################################################
################ Generate lists of genes for finding orthologs #################################################
################################################################################################################

# output_mean <- readRDS("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_10kb_windows_meanPvalue.rds")
# output_min <- readRDS("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_10kb_windows_minPvalue.rds")
## Read in the results from annotating SNPs
## Using a pvalue cutoff (of 0.0001) to ID significant SNPs
## Returns ~1k genes per comparison (only 216 for peaks)

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






## The below works, but should be run on ortholog lists (at least oreochromis -> zebrafish)


### Run some GSEA?
library(org.Dr.eg.db)
library(fgsea)
library(msigdbr)

gene_list <- genes[[4]]$ps
names(gene_list) <- as.character(genes[[4]]$gene)
gene_list <- gene_list[order(gene_list, decreasing = T)]
gene_list <- gene_list[!(is.na(names(gene_list)))]

gene_list <- gene_list[!(duplicated(names(gene_list)))]
# Need to log transform the list (and inverse it) so that the results make sense
gene_list <- log(gene_list)*-1

## Load the go lists using package msigdbr, which is really nice https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
## This has go lists for zebrafish, which are associated with gene_symbols, entrez IDs, and ensembl IDs, as well as human versions (symbol, entrez, and ensembl) (presumably could use mouse to get mouse and human)

all_gene_sets = msigdbr(species = "Danio rerio")

## Make these into a list based on the gs_name
all_gene_sets_list <- split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)

## Or can subset to a specific kind of GO term
## gs_cat and gs_subcat (where "CP:KEGG" is)
kegg <- all_gene_sets[all_gene_sets$gs_subcat == "CP:KEGG",]
kegg_list <- split(x = kegg$gene_symbol, f = kegg$gs_name)

# Then can run gsea analysis using fgsea, the pathway (list), and gene_list
fgseaRes <- fgsea(pathways = kegg_list, gene_list, minSize=15, maxSize=500, scoreType = "pos")
View(fgseaRes)

fgseaRes_all <- fgsea(pathways = all_gene_sets_list, gene_list, minSize=15, maxSize=500, scoreType = "pos")
View(fgseaRes_all)

## These plot the results, by taking a specific list and the gene_list
plotEnrichment(kegg_list[["KEGG_AXON_GUIDANCE"]], gene_list) + labs(title="KEGG_AXON_GUIDANCE")
plotEnrichment(kegg_list[["KEGG_ERBB_SIGNALING_PATHWAY"]], gene_list) + labs(title="KEGG_ERBB_SIGNALING_PATHWAY")

plotEnrichment(all_gene_sets_list[["DACOSTA_UV_RESPONSE_VIA_ERCC3_DN"]], gene_list) + labs(title="DACOSTA_UV_RESPONSE_VIA_ERCC3_DN")

## Can look at top ones like this to find the names for plotting some example graphs
head(fgseaRes_all$pathway[order(fgseaRes_all$padj, decreasing = FALSE)])
head(fgseaRes$pathway[order(fgseaRes$padj, decreasing = FALSE)])

plotEnrichment(kegg_list[["KEGG_FOCAL_ADHESION"]], gene_list) + labs(title="KEGG_FOCAL_ADHESION")

plotEnrichment(all_gene_sets_list[["GOBP_ANATOMICAL_STRUCTURE_FORMATION_INVOLVED_IN_MORPHOGENESIS"]], gene_list) + labs(title="GOBP_ANATOMICAL_STRUCTURE_FORMATION_INVOLVED_IN_MORPHOGENESIS")


### For gsea, the pathways object needs to be a named list of IDs that matches the names of the gene_list
### YES! THIS WORKS!




