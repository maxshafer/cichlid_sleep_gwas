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

# # This can be run from the command line with trailing arguments, or within IDE by setting them manually
# # Load arguments
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 3) {
#   stop("must have 3 arguments for Rscript command")
# }

# # Setting arguments for testing script
# # Need to make sure to specify the same phenotype file as was used in the gwas (if testing fewer species)
args <- c("../../pheno_data/combined_cichlid_data_2022-11-17_new_peaks.csv", 0.005, "peak_dawn")

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


################################################################################################################################################################################################
#### Find SNPS of interest and annotate genes and exons  #######################################################################################################################################
################################################################################################################################################################################################

# We only care about SNPs that had low p-values in both experiments
snps <- intersect(gwas.perchr$location, pgls.perchr$location)
# snps <- intersect( slice_min(gwas.perschr, ps, n=100000) %>% select(location),  slice_min(pgls.perchr, PGLSps, n=100000) %>% select(location) ) # this would also work if I want a certain number of snps

# Make one dataframe and remove the extras to free up memory
perchr <- merge(gwas.perchr[gwas.perchr$location %in% snps,], pgls.perchr[pgls.perchr$location %in% snps,])
# rm(gwas.perchr, pgls.perchr)

# Preload all allelefreqs so I can use it also to make plots for the exonic genes (or at least fetch the REF/ALT nucleotides)

if (!(exists("allele_freq"))) {
  allele_freq <- lapply(unique(perchr$chr), function(x) fread(paste("speciesCodesAndGenomeIDs_for_GWAS_2_cohort_db_geno_", x, "_fromProbabilities_AF.txt", sep = "")))
  names(allele_freq) <- unique(perchr$chr)
}

################################################################################################################################################################################################
#### Load SNP files from previous script  ######################################################################################################################################################
################################################################################################################################################################################################

genes1 <- read.csv(file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "genes_ALL_SNPS.csv", sep = "_"), row.names = "X", header = T)

exon_genes_2 <- read.csv(file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "SNPSinExons.csv", sep = "_"), row.names = "X", header = T)

################################################################################################################################################################################################
#### Find tippy top SNPs and make plots  #######################################################################################################################################################
################################################################################################################################################################################################

# snps_of_interest <- dplyr::slice_min(perchr, ps+PGLSps, n=50) %>% dplyr::select(location)

snps_of_interest <- union(union(dplyr::slice_min(perchr, ps, n=100) %>% dplyr::select(location), dplyr::slice_min(perchr, PGLSps, n=100) %>% dplyr::select(location)), dplyr::slice_max(perchr, (-1*log(ps))+(-1*log(PGLSps)), n=100) %>% dplyr::select(location))

snps_of_interest <- snps_of_interest$location

# genes <- getGenesEns(results.table = perchr[perchr$location %in% snps_of_interest,], gene_table = gtf2, window = 25000, windowed = TRUE, combined = FALSE)
# genes$location <- paste(genes$CHROM, genes$POS, sep = ":")
genes_for_plotting <- genes1[genes1$location %in% snps_of_interest & genes1$closest == "yes" & genes1$feature == "gene",]

genes_for_plotting$label <- paste(genes_for_plotting$gene, genes_for_plotting$product, sep = " - ")

perchr$label <- genes_for_plotting$label[match(perchr$location, genes_for_plotting$location)]

gwas.pgls.plot <- ggplot(perchr, aes(x = log(ps), y = log(PGLSps), label = label)) + geom_point() + geom_text_repel() + theme_classic() + scale_x_reverse() + scale_y_reverse()

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "psVsPGLSps_plot.png", sep = "_"), width = 10, height = 10, units = "in", res = 750)
gwas.pgls.plot
dev.off()

################################################################################################################################################################################################
#### Make Manhattan Plots  #####################################################################################################################################################################
################################################################################################################################################################################################

# Add labels to data.frames

gwas.perchr$label <- genes_for_plotting$gene[match(gwas.perchr$location, genes_for_plotting$location)]
pgls.perchr$label <- genes_for_plotting$gene[match(pgls.perchr$location, genes_for_plotting$location)]
perchr$label2 <- genes_for_plotting$gene[match(perchr$location, genes_for_plotting$location)]

gwas.manhattan <- ggplot(gwas.perchr, aes(x = coord, y = -1*log(ps), colour = chr, label = label)) + geom_point() + geom_text_repel(size = 3, colour = "black", point.padding = 1) + facet_grid(~chr, scales = "free", space = "free", switch = "x") + theme_classic() + theme(panel.spacing = unit(0.5, "lines"), axis.text.x = element_text(angle = 90), legend.position = "none") + coord_cartesian(clip = "off")
pgls.manhattan <- ggplot(pgls.perchr, aes(x = coord, y = log(PGLSps), colour = chr, label = label)) + geom_point() + geom_text_repel(size = 3, colour = "black", point.padding = 1) + facet_grid(~chr, scales = "free", space = "free") + theme_classic() + theme(panel.spacing = unit(0.5, "lines"), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + coord_cartesian(clip = "off")

perchr.manhattan <- ggplot(perchr, aes(x = coord, y = (-1*log(ps))+(-1*log(PGLSps)), colour = chr, label = label2)) + geom_point() + geom_text_repel(size = 3, colour = "black", point.padding = 1) + facet_grid(~chr, scales = "free", space = "free", switch = "x") + theme_classic() + theme(panel.spacing = unit(0.5, "lines"), axis.text.x = element_text(angle = 90), legend.position = "none") + coord_cartesian(clip = "off")


png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "manhattan_plots.png", sep = "_"), width = 20, height = 10, units = "in", res = 750)
gwas.manhattan + pgls.manhattan + plot_layout(nrow = 2)
dev.off()

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "manhattan_plot_perchr.png", sep = "_"), width = 20, height = 5, units = "in", res = 750)
perchr.manhattan
dev.off()

# ## Make a plot for presentations
# genes_for_plotting <- genes_for_plotting[genes_for_plotting$gene == "LOC100711869",]
# gwas.perchr$label <- genes_for_plotting$gene[match(gwas.perchr$location, genes_for_plotting$location)]
# 
# gwas.manhattan <- ggplot(gwas.perchr2, aes(x = coord, y = -1*log(ps), colour = chr, label = label)) + geom_point() + geom_text_repel(size = 3, colour = "black", point.padding = 1) + facet_grid(~chr, scales = "free", space = "free", switch = "x") + theme_classic() + theme(panel.spacing = unit(0.5, "lines"), axis.text.x = element_text(angle = 90), legend.position = "none") + coord_cartesian(clip = "off")
# 
# perchr$label <- genes_for_plotting$gene[match(perchr$location, genes_for_plotting$location)]
# perchr2 <- perchr[perchr$chr != "UNPLACED",]
# 
# perchr.manhattan <- ggplot(perchr2, aes(x = coord, y = (-1*log(ps))+(-1*log(PGLSps)), colour = chr, label = label)) + geom_point() + facet_grid(~chr, scales = "free", space = "free", switch = "x") + theme_classic() + theme(panel.spacing = unit(0.1, "lines"), axis.text.x = element_text(angle = 90), legend.position = "none") + coord_cartesian(clip = "off")
# 
# perchr.manhattan <- perchr.manhattan + theme(axis.title = element_text(size = 20)) + xlab("Chromosomal position") + ylab("GWAS plus PGLS p-value") #+ geom_text(size = 5, colour = "black", nudge_y = 1)
# 
# png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "manhattan_plot_perchr_presentation.png", sep = "_"), width = 20, height = 5, units = "in", res = 750)
# perchr.manhattan
# dev.off()

################################################################################################################################################################################################
#### Find the allele frequencies for the top SNPs  #############################################################################################################################################
################################################################################################################################################################################################

# Load phenotype file for plotting and annotating with frequencies
phenotypes <- read.csv(Phenofile, row.names = "X")

if (grepl("day_night_dif", pheno_to_test)) {
  pheno_to_test2 <- "day_night_dif"
} else {
  pheno_to_test2 <- pheno_to_test
}

plots <- list()
for (i in 1:nrow(genes_for_plotting)) {
  gene <- genes_for_plotting$gene[i]
  chromosome <- genes_for_plotting$CHROM[i]
  coordinate <- genes_for_plotting$POS[i]
  
  # if (!(exists("allele_freq"))) {
  #   allele_freq <- fread(paste("speciesCodesAndGenomeIDs_for_GWAS_2_cohort_db_geno_",chromosome, "_fromProbabilities_AF.txt", sep = ""))
  # } else {
  #   if (chromosome %in% allele_freq$chr[1]) {
  #     allele_freq <- allele_freq
  #   } else {
  #     allele_freq <- fread(paste("speciesCodesAndGenomeIDs_for_GWAS_2_cohort_db_geno_",chromosome, "_fromProbabilities_AF.txt", sep = ""))
  #   }
  # }
  
  af <- allele_freq[[chromosome]]
  
  phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] <- t(af[af$coord == coordinate,5:ncol(af)])[match(phenotypes$six_letter_name_Ronco, row.names(t(af[af$coord == coordinate,5:ncol(af)])))]
  
  phenotypes$Genotype <- as.character(ifelse(phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] < 0.15, 0, ifelse(phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] > 0.15 & phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] < 0.35, 0.25, ifelse(phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] > 0.35 & phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] < 0.65, 0.5, ifelse(phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] > 0.65 & phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])] < 0.85, 0.75, 1.0)))))
  phenotypes <- phenotypes[!(is.na(phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])])),]
  
  lm <- lm(phenotypes[,pheno_to_test2] ~ phenotypes[,paste0("genotype_", genes_for_plotting$gene[i])])
  direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
  
  ref_allele <- as.character(af[af$coord == coordinate,3])
  alt_allele <- as.character(af[af$coord == coordinate,4])
  
  genes_for_plotting$ref_allele[i] <- ref_allele
  genes_for_plotting$alt_allele[i] <- alt_allele
  
  plots[[i]] <- ggplot(phenotypes, aes_string(x = "Genotype", y = pheno_to_test2, fill = "Genotype")) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic()  + scale_fill_brewer(palette = "RdBu", direction = direction)
  plots[[i]] <- plots[[i]] + ggtitle(paste("GWAS for ", gene, "\nPosition: ", chromosome, ":", coordinate, "   Ref: ", ref_allele, ", Alt: ", alt_allele, sep = "")) + theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20)) + ylab(pheno_to_test)
}

names(plots) <- genes_for_plotting$gene

# Just plot the top 50 by combined p-values
snps_of_interest_final <- dplyr::slice_max(perchr, (-1*log(ps))+(-1*log(PGLSps)), n=50) %>% dplyr::select(location)
snps_of_interest_final <- snps_of_interest_final$location
genes_for_plotting_2 <- genes_for_plotting[genes_for_plotting$location %in% snps_of_interest_final,]

write.csv(genes_for_plotting, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "genes_for_plotting.csv", sep = "_"))

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "ExampleGenes.png", sep = "_"), width = 25, height = 25, units = "in", res = 500)
wrap_plots(plots[genes_for_plotting_2$gene])
dev.off()

pdf(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "ExampleGenes.pdf", sep = "_"), width = 25, height = 25)
wrap_plots(plots[genes_for_plotting_2$gene])
dev.off()

exon_genes_3 <- exon_genes_2[!(duplicated(exon_genes_2$location)) & exon_genes_2$wobble == F,]
exon_genes_3 <- dplyr::slice_min(exon_genes_3, ps+PGLSps, n=100)

plots <- list()
for (i in 1:nrow(exon_genes_3)) {
  gene <- exon_genes_3$gene[i]
  chromosome <- exon_genes_3$CHROM[i]
  coordinate <- exon_genes_3$POS[i]
  
  af <- allele_freq[[chromosome]]
  
  phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] <- t(af[af$coord == coordinate,5:ncol(af)])[match(phenotypes$six_letter_name_Ronco, row.names(t(af[af$coord == coordinate,5:ncol(af)])))]
  
  phenotypes$Genotype <- as.character(ifelse(phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] < 0.15, 0, ifelse(phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] > 0.15 & phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] < 0.35, 0.25, ifelse(phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] > 0.35 & phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] < 0.65, 0.5, ifelse(phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] > 0.65 & phenotypes[,paste0("genotype_", exon_genes_3$gene[i])] < 0.85, 0.75, 1.0)))))
  phenotypes <- phenotypes[!(is.na(phenotypes[,paste0("genotype_", exon_genes_3$gene[i])])),]
  
  lm <- lm(phenotypes[,pheno_to_test2] ~ phenotypes[,paste0("genotype_", exon_genes_3$gene[i])])
  direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
  
  ref_allele <- as.character(af[af$coord == coordinate,3])
  alt_allele <- as.character(af[af$coord == coordinate,4])
  
  exon_genes_3$ref_allele[i] <- ref_allele
  exon_genes_3$alt_allele[i] <- alt_allele
  
  plots[[i]] <- ggplot(phenotypes, aes_string(x = "Genotype", y = pheno_to_test2, fill = "Genotype")) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "RdBu", direction = direction)
  plots[[i]] <- plots[[i]] + ggtitle(paste("GWAS for ", gene, "\nPosition: ", chromosome, ":", coordinate, "   Ref: ", ref_allele, ", Alt: ", alt_allele, sep = "")) + theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20)) + ylab(pheno_to_test)
}

names(plots) <- exon_genes_3$gene

write.csv(exon_genes_3, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "SNPSinExons_nonwobble.csv", sep = "_"))

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "ExampleExonGenes.png", sep = "_"), width = 75, height = 75, units = "in", res = 500)
wrap_plots(plots)
dev.off()

pdf(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "ExampleExonGenes.pdf", sep = "_"), width = 75, height = 75)
wrap_plots(plots)
dev.off()


# phenotypes is also appended with the genotype for each of the genes so we can save that
write.csv(phenotypes, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS", pheno_to_test, cutoff, "allele_freqs_topSNPs.csv", sep = "_"))
