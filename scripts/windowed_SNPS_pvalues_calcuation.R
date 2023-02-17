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



