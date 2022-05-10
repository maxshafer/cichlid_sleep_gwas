library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)

## Make a script to read in all of the allelefreq files, combine them, plot the relationships, then subset by p-values
## This can then be modified, and re-run if we want other plots produced!

source('/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R')

# This can be run from the command line with trailing arguments, or within IDE by setting them manually
# # Load arguments
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 2) {
#   stop("must have 3 arguments for Rscript command")
# }

# Setting arguments for testing script
args <- c("../cichlid_sleep_gwas/pheno_data/combined_cichlid_data_2022-05-05.csv", 0.001, "peak")

Phenofile <- args[1]
cutoff <- as.numeric(args[2])
pheno_to_test <- as.character(args[3])

# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/alleleFreqs")

################################################################################################################################################################################################
#### Load the data  ############################################################################################################################################################################
################################################################################################################################################################################################

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene",]
gtf3 <- gtf[gtf$feature == "mRNA",]
gtf4 <- gtf[gtf$feature == "exon",]
rm(gtf)

# Retrieve dataset names
# Should have files per chromosome, per trait, and per test (pi or PGLSpi)
if (pheno_to_test == "peak") {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
  datasets <- datasets[!(grepl("amplitude", datasets))]
} else {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
}

datasets.gwas <- datasets[grep("_piVals", datasets)]
datasets.pgls <- datasets[grep("_PGLSpiVals", datasets)]

# The above are all of the GWAS results per chromosome
# Read in and combine into one dataframe for plotting/subsetting

gwas.perchr <- lapply(datasets.gwas, function(x) fread(x, showProgress = F, drop = "V1"))
gwas.perchr <- lapply(gwas.perchr, function(x) x[x$ps < cutoff,]) # this keeps only the top ~10% of snps
gwas.perchr <- Reduce(rbind, gwas.perchr)
colnames(gwas.perchr) <- c("chr", "coord", "ps")
gwas.perchr$location <- paste(gwas.perchr$chr, gwas.perchr$coord, sep = ":")

pgls.perchr <- lapply(datasets.pgls, function(x) fread(x, showProgress = F, drop = "V1"))
pgls.perchr <- lapply(pgls.perchr, function(x) {
  colnames(x) <- c("chr", "coord", "PGLSps")
  return(x)
})
pgls.perchr <- lapply(pgls.perchr, function(x) x[x$PGLSps < cutoff,]) # this keeps only the top ~10% of snps
pgls.perchr <- Reduce(rbind, pgls.perchr)
pgls.perchr$location <- paste(pgls.perchr$chr, pgls.perchr$coord, sep = ":")

# # Make sure to use the intersect of what is left
# locations <- intersect(gwas.perchr$location, pgls.perchr$location)
# gwas.perchr <- gwas.perchr[gwas.perchr$location %in% locations,]
# pgls.perchr <- pgls.perchr[pgls.perchr$location %in% locations,]

################################################################################################################################################################################################
#### Find SNPS of interest and annotate genes and exons  #######################################################################################################################################
################################################################################################################################################################################################

# We only care about SNPs that were in both experiments (some were removed by QC before statistical analysis)
snps <- intersect(gwas.perchr$location, pgls.perchr$location) # This keeps about 300k snps
# snps <- intersect( slice_min(gwas.perchr, ps, n=100000) %>% select(location),  slice_min(pgls.perchr, PGLSps, n=n=100000) %>% select(location) ) # this would also work if I want a certain number of snps

# Make one dataframe and remove the extras to free up memory
perchr <- merge(gwas.perchr[gwas.perchr$location %in% snps,], pgls.perchr[pgls.perchr$location %in% snps,])
# rm(gwas.perchr, pgls.perchr)

# Find the top 1000 SNPs by lowest combined p-values
snps_of_interest <- dplyr::slice_min(perchr, ps+PGLSps, n=1000) %>% dplyr::select(location) # This finds the SNPS with the lowest combined value across both tests
snps_of_interest <- snps_of_interest$location

# Annotate SNPs with close genes and save out csv
genes <- getGenesEns(results.table = perchr[perchr$location %in% snps_of_interest,], gene_table = gtf2, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE)
genes$location <- paste(genes$CHROM, genes$POS, sep = ":")
write.csv(genes, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_genes.csv", sep = ""))

# Find SNPs in exons only (also includes untranslated regions) and save out csv
exon_genes <- getGenesEns(results.table = perchr[perchr$location %in% snps_of_interest,], gene_table = gtf4, gene_table_2 = gtf3, window = 25000, windowed = TRUE, combined = FALSE)
exon_genes <- exon_genes[exon_genes$distance_to_gene == 0, ]
exon_genes <- exon_genes[!(duplicated(exon_genes$gene)),]
exon_genes$location <- paste(exon_genes$CHROM, exon_genes$POS, sep = ":")
write.csv(exon_genes, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_SNPSinExons.csv", sep = ""))

################################################################################################################################################################################################
#### Find tippy top SNPs and make plots  #######################################################################################################################################################
################################################################################################################################################################################################

snps_of_interest <- dplyr::slice_min(perchr, ps+PGLSps, n=20) %>% dplyr::select(location)
snps_of_interest <- snps_of_interest$location

genes <- getGenesEns(results.table = perchr[perchr$location %in% snps_of_interest,], gene_table = gtf2, window = 25000, windowed = TRUE, combined = FALSE)
genes$location <- paste(genes$CHROM, genes$POS, sep = ":")
genes <- genes[genes$closest == "yes",]
genes$label <- paste(genes$gene, genes$product, sep = " - ")

perchr$label <- genes$label[match(perchr$location, genes$location)]

gwas.pgls.plot <- ggplot(perchr, aes(x = log(ps), y = log(PGLSps), label = label)) + geom_point() + geom_text_repel() + theme_classic() + scale_x_reverse() + scale_y_reverse()

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_psVsPGLSps_plot.png", sep = ""), width = 10, height = 10, units = "in", res = 750)
gwas.pgls.plot
dev.off()

################################################################################################################################################################################################
#### Make Manhattan Plots  #####################################################################################################################################################################
################################################################################################################################################################################################

# Add labels to data.frames

gwas.perchr$label <- genes$gene[match(gwas.perchr$location, genes$location)]
pgls.perchr$label <- genes$gene[match(pgls.perchr$location, genes$location)]

gwas.manhattan <- ggplot(gwas.perchr, aes(x = coord, y = -1*log(ps), colour = chr, label = label)) + geom_point() + geom_text_repel(size = 5, colour = "black", point.padding = 1) + facet_grid(~chr, scales = "free", space = "free", switch = "x") + theme_classic() + theme(panel.spacing = unit(0.5, "lines"), axis.text.x = element_text(angle = 90), legend.position = "none") + coord_cartesian(clip = "off")
pgls.manhattan <- ggplot(pgls.perchr, aes(x = coord, y = log(PGLSps), colour = chr, label = label)) + geom_point() + geom_text_repel(size = 5, colour = "black", point.padding = 1) + facet_grid(~chr, scales = "free", space = "free") + theme_classic() + theme(panel.spacing = unit(0.5, "lines"), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + coord_cartesian(clip = "off")

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_manhattan_plots.png", sep = ""), width = 20, height = 10, units = "in", res = 750)
gwas.manhattan + pgls.manhattan + plot_layout(nrow = 2)
dev.off()

################################################################################################################################################################################################
#### Find the allele frequencies for the top SNPs  #############################################################################################################################################
################################################################################################################################################################################################

# Load phenotype file for plotting and annotating with frequencies
phenotypes <- read.csv("../cichlid_sleep_gwas/pheno_data/combined_cichlid_data_2022-05-05.csv", row.names = "X")

# Sort by Chromosome for faster finding (loads fewer allele_freq tables)
genes <- genes[order(genes$CHROM),]

plots <- list()
for (i in 1:nrow(genes)) {
  gene <- genes$gene[i]
  chromosome <- genes$CHROM[i]
  coordinate <- genes$POS[i]
  
  if (!(exists("allele_freq"))) {
    allele_freq <- fread(paste("speciesCodesAndGenomeIDs_for_GWAS_2_cohort_db_geno_",chromosome, "_fromProbabilities_AF.txt", sep = ""))
  } else {
    if (chromosome %in% allele_freq$chr[1]) {
      allele_freq <- allele_freq
    } else {
      allele_freq <- fread(paste("speciesCodesAndGenomeIDs_for_GWAS_2_cohort_db_geno_",chromosome, "_fromProbabilities_AF.txt", sep = ""))
    }
  }

  phenotypes[,paste0("genotype_", genes$gene[i])] <- t(allele_freq[allele_freq$coord == coordinate,5:ncol(allele_freq)])[match(phenotypes$six_letter_name_Ronco, row.names(t(allele_freq[allele_freq$coord == coordinate,5:ncol(allele_freq)])))]
  
  phenotypes$Genotype <- as.character(ifelse(phenotypes[,paste0("genotype_", genes$gene[i])] < 0.15, 0, ifelse(phenotypes[,paste0("genotype_", genes$gene[i])] > 0.15 & phenotypes[,paste0("genotype_", genes$gene[i])] < 0.35, 0.25, ifelse(phenotypes[,paste0("genotype_", genes$gene[i])] > 0.35 & phenotypes[,paste0("genotype_", genes$gene[i])] < 0.65, 0.5, ifelse(phenotypes[,paste0("genotype_", genes$gene[i])] > 0.65 & phenotypes[,paste0("genotype_", genes$gene[i])] < 0.85, 0.75, 1.0)))))
  phenotypes <- phenotypes[!(is.na(phenotypes[,paste0("genotype_", genes$gene[i])])),]
  
  lm <- lm(phenotypes$day_night_dif ~ phenotypes[,paste0("genotype_", genes$gene[i])])
  direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
  
  ref_allele <- as.character(allele_freq[allele_freq$coord == coordinate,3])
  alt_allele <- as.character(allele_freq[allele_freq$coord == coordinate,4])
  
  plots[[i]] <- ggplot(phenotypes, aes_string(x = "Genotype", y = pheno_to_test, fill = "Genotype")) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic()  + scale_fill_brewer(palette = "RdBu", direction = direction)
  plots[[i]] <- plots[[i]] + ggtitle(paste("GWAS for ", gene, "\nPosition: ", chromosome, ":", coordinate, "   Ref: ", ref_allele, ", Alt: ", alt_allele, sep = "")) + theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20)) + ylab(pheno_to_test)
}

names(plots) <- genes$gene

png(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_ExampleGenes.png", sep = ""), width = 25, height = 25, units = "in", res = 500)
wrap_plots(plots)
dev.off()

pdf(paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_ExampleGenes.pdf", sep = ""), width = 25, height = 25)
wrap_plots(plots)
dev.off()

# phenotypes is also appended with the genotype for each of the genes so we can save that
write.csv(phenotypes, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/GWAS_vs_PGLS_", pheno_to_test, "_allele_freqs_topSNPs.csv", sep = ""))



