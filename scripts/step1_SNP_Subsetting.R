library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)

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

# This can be run from the command line with trailing arguments, or within IDE by setting them manually
# Load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("must have 3 arguments for Rscript command")
}

# # Setting arguments for testing script
# # Need to make sure to specify the same phenotype file as was used in the gwas (if testing fewer species)
# args <- c("../../pheno_data/combined_cichlid_data_2022-08-04.csv", 0.00005, "total_rest")

Phenofile <- args[1]
cutoff <- as.numeric(args[2])
pheno_to_test <- as.character(args[3])

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

datasets.gwas <- datasets[grep("_piVals", datasets)]
datasets.pgls <- datasets[grep("_PGLSpiVals", datasets)]

gwas.names <- substr(datasets.gwas, start = 52, stop = 60)
pgls.names <- substr(datasets.pgls, start = 52, stop = 60)

# The above are all of the GWAS results per chromosome
# Read in and combine into one dataframe for plotting/subsetting

gwas.perchr <- lapply(datasets.gwas, function(x) fread(x, showProgress = F, drop = "V1"))
gwas.perchr <- lapply(gwas.perchr, function(x) x[x$ps < cutoff,]) # this keeps only the top ~10% of snps
names(gwas.perchr) <- gwas.names
gwas.perchr <- Reduce(rbind, gwas.perchr)
colnames(gwas.perchr) <- c("chr", "coord", "ps")
gwas.perchr$location <- paste(gwas.perchr$chr, gwas.perchr$coord, sep = ":")

pgls.perchr <- lapply(datasets.pgls, function(x) fread(x, showProgress = F, drop = "V1"))
pgls.perchr <- lapply(pgls.perchr, function(x) {
  colnames(x) <- c("chr", "coord", "PGLSps")
  return(x)
})
pgls.perchr <- lapply(pgls.perchr, function(x) x[x$PGLSps < cutoff,]) # this keeps only the top ~10% of snps
names(pgls.perchr) <- pgls.names
pgls.perchr <- Reduce(rbind, pgls.perchr)
pgls.perchr$location <- paste(pgls.perchr$chr, pgls.perchr$coord, sep = ":")

# # Make sure to use the intersect of what is left
# locations <- intersect(gwas.perchr$location, pgls.perchr$location)
# gwas.perchr <- gwas.perchr[gwas.perchr$location %in% locations,]
# pgls.perchr <- pgls.perchr[pgls.perchr$location %in% locations,]

################################################################################################################################################################################################
#### Find SNPS of interest and annotate genes and exons  #######################################################################################################################################
################################################################################################################################################################################################

# We only care about SNPs that had low p-values in both experiments
snps <- intersect(gwas.perchr$location, pgls.perchr$location)
# snps <- intersect( slice_min(gwas.perschr, ps, n=100000) %>% select(location),  slice_min(pgls.perchr, PGLSps, n=100000) %>% select(location) ) # this would also work if I want a certain number of snps

# Make one dataframe and remove the extras to free up memory
perchr <- merge(gwas.perchr[gwas.perchr$location %in% snps,], pgls.perchr[pgls.perchr$location %in% snps,])

rm(gwas.perchr, pgls.perchr)

# Annotate all the snps that pass the cutoff for more Systemies biology
snps_of_interest <- perchr$location

# Annotate SNPs with close genes and save out csv
genes <- getGenesEns(results.table = perchr[perchr$location %in% snps_of_interest,], gene_table = gtf, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE)
genes$location <- paste(genes$CHROM, genes$POS, sep = ":")
genes$comb_value_log_p <- -1*log(genes$ps + genes$PGLSps)

genes1 <- genes[genes$feature %in% c("gene", "mRNA", "exon", "lnc_RNA", "mRNA", "pseudogene", "transcript", "tRNA", "V_gene_segment"),]
length(unique(genes1[genes1$distance_to_gene == 0, "location"])) # the number of unique SNPs that have a 'hit'

write.csv(genes1, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "genes_ALL_SNPS.csv", sep = "_"))

################################################################################################################################################################################################
#### Calculate % of snps in various categories  ################################################################################################################################################
################################################################################################################################################################################################

genes2 <- genes1[genes1$distance_to_gene == 0,]
length(unique(genes2[genes2$distance_to_gene == 0, "location"]))
# [1] 2399 # out of 3406 are found within 0bp of a feature, the rest should be intergenic (1007), some of which are within the cutoff distance to a gene, and some that are not

SNP_gene_table <- as.data.frame.matrix(table(genes2$location, genes2$feature))

exonic_snps <- row.names(SNP_gene_table[SNP_gene_table$exon > 0 & SNP_gene_table$gene > 0, ]) # 507 exonic
intronic_snps <- row.names(SNP_gene_table[SNP_gene_table$exon == 0 & SNP_gene_table$gene > 0, ]) # 1890 in gene but not exon, so intronic?
pseudo_snps <- row.names(SNP_gene_table[SNP_gene_table$gene == 0 & SNP_gene_table$pseudogene > 0, ]) # 1890 in gene but not exon, so intronic?
intergenic_snps <- setdiff(snps_of_interest, genes$location) # 106
intergenic_close_snps <- setdiff(genes$location, genes2$location) # 901

df <- data.frame(category = c("exonic_snps", "intronic_snps", "pseudo_snps", "intergenic_snps", "intergenic_close_snps"), values = unlist(lapply(list(exonic_snps, intronic_snps, pseudo_snps, intergenic_snps, intergenic_close_snps), function(x) length(x))))
df <- df %>% arrange(desc(category)) %>% mutate(prop = values / sum(df$values) *100) %>% mutate(ypos = cumsum(prop)- 0.5*prop)
df <- df %>% mutate(values2 = values / sum(df$values) * 100)
pie.chart <- ggplot(df, aes(x = "", y = values2, fill = category)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) + theme_void() + theme(legend.position = "none") + geom_text(aes(y = ypos, label = category)) #+ theme_classic() + theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "PieChart_plot.png", sep = "_"), width = 5, height = 5, units = "in", res = 750)
pie.chart
dev.off()

################################################################################################################################################################################################
#### Look at the best hits, extract exonic non-synonomous changes  #############################################################################################################################
################################################################################################################################################################################################


# # Find the top 1000 SNPs by lowest combined p-values
# snps_of_interest <- dplyr::slice_min(perchr, ps+PGLSps, n=1000) %>% dplyr::select(location) # This finds the SNPS with the lowest combined value across both tests
# snps_of_interest <- snps_of_interest$location
# 
# # Annotate SNPs with close genes and save out csv
# genes <- getGenesEns(results.table = perchr[perchr$location %in% snps_of_interest,], gene_table = gtf2, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE)
# genes$location <- paste(genes$CHROM, genes$POS, sep = ":")
# write.csv(genes, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_genes.csv", sep = ""))

# Find SNPs in exons only (also includes untranslated regions) and save out csv
# exon_genes <- getGenesEns(results.table = perchr[perchr$location %in% snps_of_interest,], gene_table = gtf4, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE)
exon_genes <- genes2[genes2$feature == "exon",]


# exon_genes <- exon_genes[exon_genes$distance_to_gene == 0, ]
exon_genes <- exon_genes[!(duplicated(exon_genes$location)),]

## Find out if they are in wobble or whichever position
## I could definitely annotate the SNPs in exons as being either 1st, 2nd or 3rd/wobble position, find ones in 'CDS', if frame == +, then add the frame value, if - then subtract it, if .33 then first, .66 then 2nd, if .0 then wobble. (frame is how many nucleotides ahead the frame is)
## Have to pull the fasta file to get the AA change, which could be done, but then, probably faster to download the partial vcf file (with just the SNP and anno data), and load/query it ::shrugemoji::

# Take each row of exon_genes, find the CDS from gtf that matches, and calculate the aa position

exon_genes_2 <- wobbleAA(results.table = exon_genes, gene_table = gtf4, window = 0)

write.csv(exon_genes_2, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "SNPSinExons.csv", sep = "_"))

## The rest, plotting the manhattan plots, and saving the allele frequencies, is done in a separate script.
