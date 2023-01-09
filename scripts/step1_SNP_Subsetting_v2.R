library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)
library(fdrtool)

# list = ls()
# rm(list = list[!(list %in% c("allele_freq", "gtf", "gtf2", "gtf3", "gtf4"))])

# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/sra_reads_nobackup/alleleFreqs")

## Make a script to read in all of the allelefreq files, combine them, plot the relationships, then subset by p-values
## This can then be modified, and re-run if we want other plots produced!

source('/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R')

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene",]
gtf3 <- gtf[gtf$feature == "mRNA",]
gtf4 <- gtf[gtf$feature == "CDS",]
# rm(gtf)

# # This can be run from the command line with trailing arguments, or within IDE by setting them manually
# # Load arguments
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 3) {
#   stop("must have 3 arguments for Rscript command")
# }

# Setting arguments for testing script
# Need to make sure to specify the same phenotype file as was used in the gwas (if testing fewer species)
args <- c("../../pheno_data/combined_cichlid_data_2022-08-04.csv", 0.00005, "day_night_dif_55-species", "pgls")

Phenofile <- args[1]
cutoff <- as.numeric(args[2])
pheno_to_test <- as.character(args[3])
statistical_test <- args[4]


################################################################################################################################################################################################
#### Load the data  ############################################################################################################################################################################
################################################################################################################################################################################################

# Retrieve dataset names
# Should have files per chromosome, per trait, and per test (pi or PGLSpi)
if (pheno_to_test == "peak") {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
  datasets <- datasets[!(grepl("amplitude", datasets))]
  datasets <- datasets[grep(".txt", datasets)]
} 
if (pheno_to_test == "day_night_dif" | pheno_to_test == "day_night_dif_spd") {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
  datasets <- datasets[!(grepl("-species", datasets))]
  datasets <- datasets[grep(".txt", datasets)]
}  else {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
  datasets <- datasets[grep(".txt", datasets)]
}

if (statistical_test == "pgls") {
  datasets <- datasets[grep("_PGLSpiVals", datasets)]
}

if (statistical_test == "gwas") {
  datasets <- datasets[grep("_piVals", datasets)]
}

names <- substr(datasets, start = 52, stop = 60)


# The above are all of the GWAS results per chromosome
# Read in and combine into one dataframe for plotting/subsetting

perchr <- lapply(datasets, function(x) fread(x, showProgress = F, drop = "V1"))
names(perchr) <- names
perchr <- Reduce(rbind, perchr)
colnames(perchr) <- c("chr", "coord", "pvalue")


## Run the above for all comparisons, and than use a pvalue cutoff
## Any multiple testing correction methods seem to fail because of the vast number of tests (44 million) that were run
## For example, bonferroni multiples the pvalue by the # of tests (44,000,000), which means nothing is significant (or even lower than 0.5) after correction

perchr.list <- list()
perchr.list[[1]] <- perchr.dn55
perchr.list[[2]] <- perchr.dn58
perchr.list[[3]] <- perchr.dn60
perchr.list[[4]] <- perchr.tr
perchr.list[[5]] <- perchr.peak


perchr.list <- lapply(perchr.list, function(x) {
  x$location <- paste(x$chr, x$coord, sep = ":")
  return(x)
})


saveRDS(perchr.list, file = "../per_chromosome_list.rds")

# lapply(perchr.list, function(x) dim(x[x$pvalue < 0.05,]))

perchr.list.cut <- lapply(perchr.list, function(x) x[x$pvalue < 0.05,])

saveRDS(perchr.list.cut, file = "../per_chromosome_pvalue_0.05_list.rds")

perchr <- readRDS("../per_chromosome_list.rds")

fdr_list <- lapply(perchr, function(x) fdrtool(x$pvalue, statistic = "pvalue"))
fdr_list_pct0 <- lapply(perchr, function(x) fdrtool(x$pvalue, statistic = "pvalue", cutoff.method = "pct0"))


################################################################################################################################################################################################
#### Find SNPS of interest and annotate genes and exons  #######################################################################################################################################
################################################################################################################################################################################################

### Now for each list, pull the genes
### But only those that pass the stricter cutoff
perchr.list.cut <- lapply(perchr.list, function(x) x[x$pvalue < 0.0001,])
saveRDS(perchr.list.cut, file = "../per_chromosome_pvalue_0.0001_list.rds")

# Annotate all the snps that pass the cutoff for more Systemies biology
snps_of_interest <- list()
snps_of_interest <- lapply(perchr.list.cut, function(x) x$location)

# Annotate SNPs with close genes and save out csv
genes <- lapply(seq_along(snps_of_interest), function(x) getGenesEns(results.table = perchr.list.cut[[x]], gene_table = gtf, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE))

genes1 <- lapply(genes, function(x) {
  x$location <- paste(x$CHROM, x$POS, sep = ":")
  x <- x[x$feature %in% c("gene", "mRNA", "exon", "lnc_RNA", "mRNA", "pseudogene", "transcript", "tRNA", "V_gene_segment"),]
  return(x)
})

lapply(genes1, function(x) length(unique(x[x$distance_to_gene == 0, "location"]))) # the number of unique SNPs that have a 'hit'

saveRDS(genes1, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "genes_ALL_SNPS.rds", sep = "_"))


################################################################################################################################################################################################
#### Calculate % of snps in various categories  ################################################################################################################################################
################################################################################################################################################################################################

genes2 <- lapply(genes1, function(x) x[x$distance_to_gene == 0,])
lapply(genes2, function(x) length(unique(x[x$distance_to_gene == 0, "location"])))
# [1] 2399 # out of 3406 are found within 0bp of a feature, the rest should be intergenic (1007), some of which are within the cutoff distance to a gene, and some that are not


## The below needs to be rewritten

# SNP_gene_table <- as.data.frame.matrix(table(genes2$location, genes2$feature))
# 
# exonic_snps <- row.names(SNP_gene_table[SNP_gene_table$exon > 0 & SNP_gene_table$gene > 0, ]) # 507 exonic
# intronic_snps <- row.names(SNP_gene_table[SNP_gene_table$exon == 0 & SNP_gene_table$gene > 0, ]) # 1890 in gene but not exon, so intronic?
# pseudo_snps <- row.names(SNP_gene_table[SNP_gene_table$gene == 0 & SNP_gene_table$pseudogene > 0, ]) # 1890 in gene but not exon, so intronic?
# intergenic_snps <- setdiff(snps_of_interest, genes$location) # 106
# intergenic_close_snps <- setdiff(genes$location, genes2$location) # 901
# 
# df <- data.frame(category = c("exonic_snps", "intronic_snps", "pseudo_snps", "intergenic_snps", "intergenic_close_snps"), values = unlist(lapply(list(exonic_snps, intronic_snps, pseudo_snps, intergenic_snps, intergenic_close_snps), function(x) length(x))))
# df <- df %>% arrange(desc(category)) %>% mutate(prop = values / sum(df$values) *100) %>% mutate(ypos = cumsum(prop)- 0.5*prop)
# df <- df %>% mutate(values2 = values / sum(df$values) * 100)
# pie.chart <- ggplot(df, aes(x = "", y = values2, fill = category)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) + theme_void() + theme(legend.position = "none") + geom_text(aes(y = ypos, label = category)) #+ theme_classic() + theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
# 
# png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "PieChart_plot.png", sep = "_"), width = 5, height = 5, units = "in", res = 750)
# pie.chart
# dev.off()


################################################################################################################################################################################################
#### Look at the best hits, extract exonic non-synonomous changes  #############################################################################################################################
################################################################################################################################################################################################

exon_genes <- lapply(genes2, function(x) x[x$feature == "exon",])

exon_genes <- lapply(exon_genes, function(x) x[!(duplicated(x$location)),])

## Find out if they are in wobble or whichever position
## I could definitely annotate the SNPs in exons as being either 1st, 2nd or 3rd/wobble position, find ones in 'CDS', if frame == +, then add the frame value, if - then subtract it, if .33 then first, .66 then 2nd, if .0 then wobble. (frame is how many nucleotides ahead the frame is)
## Have to pull the fasta file to get the AA change, which could be done, but then, probably faster to download the partial vcf file (with just the SNP and anno data), and load/query it ::shrugemoji::

# Take each row of exon_genes, find the CDS from gtf that matches, and calculate the aa position

exon_genes_2 <- lapply(exon_genes, function(x) wobbleAA(results.table = x, gene_table = gtf4, window = 0))

saveRDS(exon_genes_2, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "SNPSinExons.csv", sep = "_"))

## The rest, plotting the manhattan plots, and saving the allele frequencies, is done in a separate script.
