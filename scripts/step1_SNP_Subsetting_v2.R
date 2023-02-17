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
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas")

## Make a script to read in all of the allelefreq files, combine them, plot the relationships, then subset by p-values
## This can then be modified, and re-run if we want other plots produced!

source('/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R')

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]
gtf3 <- gtf[gtf$feature == "mRNA",]
gtf4 <- gtf[gtf$feature == "CDS",]
# rm(gtf)

# # This can be run from the command line with trailing arguments, or within IDE by setting them manually
# # Load arguments
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 3) {
#   stop("must have 3 arguments for Rscript command")
# }

# # Setting arguments for testing script
# # Need to make sure to specify the same phenotype file as was used in the gwas (if testing fewer species)
# args <- c("../../pheno_data/combined_cichlid_data_2022-08-04.csv", 0.00005, "day_night_dif_55-species", "pgls")
# 
# Phenofile <- args[1]
# cutoff <- as.numeric(args[2])
# pheno_to_test <- as.character(args[3])
# statistical_test <- args[4]


# ################################################################################################################################################################################################
# #### Load the data  ############################################################################################################################################################################
# ################################################################################################################################################################################################
# 
# # Retrieve dataset names
# # Should have files per chromosome, per trait, and per test (pi or PGLSpi)
# if (pheno_to_test == "peak") {
#   datasets <- list.files()[grep(pheno_to_test, list.files())]
#   datasets <- datasets[!(grepl("amplitude", datasets))]
#   datasets <- datasets[grep(".txt", datasets)]
# } 
# if (pheno_to_test == "day_night_dif" | pheno_to_test == "day_night_dif_spd") {
#   datasets <- list.files()[grep(pheno_to_test, list.files())]
#   datasets <- datasets[!(grepl("-species", datasets))]
#   datasets <- datasets[grep(".txt", datasets)]
# }  else {
#   datasets <- list.files()[grep(pheno_to_test, list.files())]
#   datasets <- datasets[grep(".txt", datasets)]
# }
# 
# if (statistical_test == "pgls") {
#   datasets <- datasets[grep("_PGLSpiVals", datasets)]
# }
# 
# if (statistical_test == "gwas") {
#   datasets <- datasets[grep("_piVals", datasets)]
# }
# 
# names <- substr(datasets, start = 52, stop = 60)
# 
# 
# # The above are all of the GWAS results per chromosome
# # Read in and combine into one dataframe for plotting/subsetting
# 
# perchr <- lapply(datasets, function(x) fread(x, showProgress = F, drop = "V1"))
# names(perchr) <- names
# perchr <- Reduce(rbind, perchr)
# colnames(perchr) <- c("chr", "coord", "pvalue")
# 
# 
# ## Run the above for all comparisons, and than use a pvalue cutoff
# ## Any multiple testing correction methods seem to fail because of the vast number of tests (44 million) that were run
# ## For example, bonferroni multiples the pvalue by the # of tests (44,000,000), which means nothing is significant (or even lower than 0.5) after correction
# 
# perchr.list <- list()
# perchr.list[[1]] <- perchr.dn55
# perchr.list[[2]] <- perchr.dn58
# perchr.list[[3]] <- perchr.dn60
# perchr.list[[4]] <- perchr.tr
# perchr.list[[5]] <- perchr.peak

# perchr.list <- readRDS("sra_reads_nobackup/per_chromosome_list.rds")

# perchr.list <- lapply(perchr.list, function(x) {
#   x$location <- paste(x$chr, x$coord, sep = ":")
#   return(x)
# })
# 
# perchr.list.10k <- lapply(perchr.list, function(x) x %>% group_by(chr, ints = cut_width(coord, width = 10000, boundary = 0)) %>% summarise(start = min(coord), end = max(coord), mean_pvalue = mean(pvalue), pvalue = min(pvalue)) )

# saveRDS(perchr.list.10k, file = "sra_reads_nobackup/per_chromosome_list_10k.rds")

# saveRDS(perchr.list, file = "../per_chromosome_list.rds")
# 
# # lapply(perchr.list, function(x) dim(x[x$pvalue < 0.05,]))
# 
# perchr.list.cut <- lapply(perchr.list, function(x) x[x$pvalue < 0.05,])
# 
# saveRDS(perchr.list.cut, file = "../per_chromosome_pvalue_0.05_list.rds")
# 
# 
# fdr_list <- lapply(perchr.list, function(x) fdrtool(x$pvalue, statistic = "pvalue"))
# fdr_list_pct0 <- lapply(perchr.list, function(x) fdrtool(x$pvalue, statistic = "pvalue", cutoff.method = "pct0"))
# 
# 
# ### There are def too many markers (44million), and this kills any multiple testing correction (basically none are significant, because I'm testing too many things)
# ### What if I do it for 1 chromosome? NC_031970
# ### Still so many per chromosome (over a million)
# ### The only way might be to reduce the SNPs I'm testing, but that's not ideal
# ### I think the best is to just use a cutoff?
# 
# perchr_chr <- lapply(perchr.list, function(x) x[x$chr %in% "NC_031969",])
# fdr_list_chr <- lapply(perchr_chr, function(x) fdrtool(x$pvalue, statistic = "pvalue"))
# bonferroni_chr <- lapply(perchr_chr, function(x) p.adjust(x$pvalue, method = "BH"))
# lapply(bonferroni_chr, function(x) min(x))
# 
# ################################################################################################################################################################################################
# #### Find SNPS of interest and annotate genes and exons  #######################################################################################################################################
# ################################################################################################################################################################################################
# 
# ### Now for each list, pull the genes
# ### But only those that pass the stricter cutoff
# p_val_cutoff <- 0.0001
# perchr.list.cut <- lapply(perchr.list, function(x) x[x$pvalue < p_val_cutoff,])
# saveRDS(perchr.list.cut, file = "../per_chromosome_pvalue_0.0001_list.rds")
# 
# # Annotate all the snps that pass the cutoff for more Systemies biology
# snps_of_interest <- list()
# snps_of_interest <- lapply(perchr.list.cut, function(x) x$location)

perchr.list.10k <- readRDS(file = "sra_reads_nobackup/per_chromosome_list_10k.rds")

## I want to run the above per chromosome (which is faster)

## First make lists of gtf2 dataframes
gtf2_chrs <- unique(gtf2$seqname)[1:23]
gtf_df <- data.frame(gtf_chrs = gtf2_chrs, vcf_chrs = as.character(substr(gtf2_chrs, 1, 9)))
list_gtf2_sub <- lapply(gtf2_chrs, function(x) gtf2[gtf2$seqname %in% x,])
names(list_gtf2_sub) <- gtf2_chrs

list_gtf3_sub <- lapply(gtf2_chrs, function(x) gtf3[gtf3$seqname %in% x,])
names(list_gtf3_sub) <- gtf2_chrs


library(DescTools)

# To id if a genes' range overlaps with the window range

## This seems to work!

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

saveRDS(output_mean, file = "")

saveRDS(output_mean, file = "/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories.rds")































## Test some things for speed
tic()
gtf_sub <- gtf[gtf$seqname %in% "NC_031965.1",]
gtf2_sub <- gtf2[gtf2$seqname %in% "NC_031965.1",]
toc()
gtf3_sub  <- gtf3[gtf3$seqname %in% "NC_031965.1",]

tic()
test.cut <- perchr.list.cut[[1]][perchr.list.cut[[1]]$chr %in% "NC_031965",]
toc()

tic()
test <- getGenesEns(results.table = test.cut, gene_table = gtf2_sub, gene_table_2 = gtf3_sub, window = 25000, windowed = TRUE, combined = FALSE)
toc() # 1.781

tic()
test2 <- getGenesEns(results.table = test.cut, gene_table = gtf2, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE)
toc() # 1.622, this one usually takes a bit longer then the above

tic()
test3 <- getGenesEns(results.table = perchr.list.cut[[1]], gene_table = gtf2_sub, gene_table_2 = gtf3_sub, window = 25000, windowed = TRUE, combined = FALSE)
toc() # 6.38, this wastes some time searching for things that don't exist


## Seems fastest to use gtf2 (if I don't care about intronic vs exonic etc), and subset both lists by chr (least headache vs fastest tradeoff - if I go smaller than chr, I will need to rerun some (to catch everything))
## But let's see
tic()
per_chrs <- unique(perchr.list.cut[[1]]$chr)
gtf2_chrs <- unique(gtf2$seqname)[1:23]
gtf2_other_chrs <- unique(gtf2$seqname)[24:length(unique(gtf2$seqname))]

list_test.cut <- lapply(per_chrs, function(x) perchr.list.cut[[1]][perchr.list.cut[[1]]$chr %in% x,])

list_gtf2_sub <- lapply(gtf2_chrs, function(x) gtf2[gtf2$seqname %in% x,])
list_gtf2_sub[[24]] <- gtf2[gtf2$seqname %in% gtf2_other_chrs,]

list_gtf3_sub <- lapply(gtf2_chrs, function(x) gtf3[gtf3$seqname %in% x,])
list_gtf3_sub[[24]] <- gtf3[gtf3$seqname %in% gtf2_other_chrs,]
toc()

tic()
test4 <- getGenesEns(results.table = perchr.list.cut[[1]], gene_table = gtf2, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE)
toc()

tic()
test5 <- lapply(seq_along(list_test.cut), function(x) {
  output <- getGenesEns(results.table = list_test.cut[[x]], gene_table = list_gtf2_sub[[x]], gene_table_2 = list_gtf3_sub[[x]], window = 25000, windowed = TRUE, combined = FALSE)
})
test5 <- Reduce(rbind, test5)
toc()

## 16.846 sec elapsed and 61.177 sec elapsed for 1:100 and 25.039 sec elapsed and 46.227 sec elapsed for the subset which includes just one chromosome
## Remains about 50% of the time, to use the subsetted gtf

## Definitely speed improvements by subsetting the gtf (even if just by chromosome)
## If I want to go smaller (which also speeds up), then I would need to double up some of the comps
## Seems really just a 50% improvement on speed, maybe I should just look into snpsnff...



# Annotate SNPs with close genes and save out csv
genes <- lapply(seq_along(snps_of_interest), function(x) getGenesEns(results.table = perchr.list.cut[[x]], gene_table = gtf, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE))

genes1 <- lapply(genes, function(x) {
  x$location <- paste(x$CHROM, x$POS, sep = ":")
  x <- x[x$feature %in% c("gene", "mRNA", "exon", "lnc_RNA", "mRNA", "pseudogene", "transcript", "tRNA", "V_gene_segment"),]
  return(x)
})

lapply(genes1, function(x) length(unique(x[x$distance_to_gene == 0, "location"]))) # the number of unique SNPs that have a 'hit'

saveRDS(genes1, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "genes_ALL_SNPS.rds", sep = "_"))

genes1 <- readRDS(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/PGLS_genes_AllCategories", "pvalue", "0.0001", "genes_ALL_SNPS.rds", sep = "_"))
################################################################################################################################################################################################
#### Calculate % of snps in various categories  ################################################################################################################################################
################################################################################################################################################################################################

genes2 <- lapply(genes1, function(x) x[x$distance_to_gene == 0,])
lapply(genes2, function(x) length(unique(x[x$distance_to_gene == 0, "location"])))


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
