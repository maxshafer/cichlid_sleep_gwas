library(ape)
library(ggtree)
library(ggplot2)
library(viridis)
library(stringr)
library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)

# Source this gist (From "https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32"), which allows using multiple of the same aes for the same plot
# Here I use it to make multiple heatmaps for Annika's figure (but should work better than the super hack I used for Figure 5 of the Cavefish single-cell paper)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/sra_reads_nobackup/alleleFreqs")
source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")
source('/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R')

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene",]
gtf3 <- gtf[gtf$feature == "mRNA",]
gtf4 <- gtf[gtf$feature == "CDS",]

# # Setting arguments for testing script
# # Need to make sure to specify the same phenotype file as was used in the gwas (if testing fewer species)

Phenofile <- "../../pheno_data/combined_cichlid_data_2022-11-17.csv"
# pheno_to_test <- "total_rest"
# pheno_to_test <- "day_night_dif_spd_60-species"
# pheno_to_test <- "day_night_dif_spd_58-species"
# pheno_to_test <- "day_night_dif_spd_55-species"
# pheno_to_test <- "day_night_dif_60-species"
# pheno_to_test <- "day_night_dif_58-species"
# pheno_to_test <- "day_night_dif_55-species"
pheno_to_test <- "peak_60-species"
# pheno_to_test <- "peak_amplitude"

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
# gwas.perchr <- lapply(gwas.perchr, function(x) x[x$ps < cutoff,]) # this keeps only the top ~10% of snps
names(gwas.perchr) <- gwas.names
gwas.perchr <- Reduce(rbind, gwas.perchr)
colnames(gwas.perchr) <- c("chr", "coord", "ps")
gwas.perchr$location <- paste(gwas.perchr$chr, gwas.perchr$coord, sep = ":")

pgls.perchr <- lapply(datasets.pgls, function(x) fread(x, showProgress = F, drop = "V1"))
pgls.perchr <- lapply(pgls.perchr, function(x) {
  colnames(x) <- c("chr", "coord", "PGLSps")
  return(x)
})
# pgls.perchr <- lapply(pgls.perchr, function(x) x[x$PGLSps < cutoff,]) # this keeps only the top ~10% of snps
names(pgls.perchr) <- pgls.names
pgls.perchr <- Reduce(rbind, pgls.perchr)
pgls.perchr$location <- paste(pgls.perchr$chr, pgls.perchr$coord, sep = ":")


################################################################################################################################################################################################
#### Find SNPS of interest and annotate genes and exons  #######################################################################################################################################
################################################################################################################################################################################################

# Make one dataframe and remove the extras to free up memory
perchr <- merge(gwas.perchr, pgls.perchr)

rm(gwas.perchr, pgls.perchr)

perchr$comb_ps <- (-1*log(perchr$ps))+(-1*log(perchr$PGLSps))

chromosomes <- unique(perchr$chr)

################################################################################################################################################################################################
#### Find tippy top SNPs and make plots  #######################################################################################################################################################
################################################################################################################################################################################################

# snps_of_interest <- dplyr::slice_min(perchr, ps+PGLSps, n=50) %>% dplyr::select(location)

cutoff <- 50

snps_of_interest <- union(union(dplyr::slice_min(perchr, ps, n=cutoff) %>% dplyr::select(location), dplyr::slice_min(perchr, PGLSps, n=cutoff) %>% dplyr::select(location)), dplyr::slice_max(perchr, comb_ps, n=cutoff) %>% dplyr::select(location))

snps_of_interest <- snps_of_interest$location

perchr_soi <- perchr[perchr$location %in% snps_of_interest,]

rm(perchr)

genes_for_plotting <- getGenesEns(results.table = perchr_soi, gene_table = gtf2, window = 25000, windowed = TRUE, combined = FALSE)
genes_for_plotting$location <- paste(genes_for_plotting$CHROM, genes_for_plotting$POS, sep = ":")
genes_for_plotting$comb_ps <- (-1*log(genes_for_plotting$ps))+(-1*log(genes_for_plotting$PGLSps))

genes_for_plotting <- genes_for_plotting[genes_for_plotting$closest == "yes" & genes_for_plotting$feature == "gene",]

genes_for_plotting <- genes_for_plotting[genes_for_plotting$gene_biotype != "tRNA",]

write.csv(genes_for_plotting, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/", pheno_to_test, "_PhyloPlots", "/GWAS_vs_PGLS_", pheno_to_test, "_top", cutoff, "_genes_for_plotting.csv", sep = ""))

# Preload all allelefreqs so I can use it also to make plots for the exonic genes (or at least fetch the REF/ALT nucleotides)

if (!(exists("allele_freq"))) {
  allele_freq <- lapply(chromosomes, function(x) fread(paste("speciesCodesAndGenomeIDs_for_GWAS_2_cohort_db_geno_", x, "_fromProbabilities_AF.txt", sep = "")))
  names(allele_freq) <- chromosomes
}

################################################################################################################################################################################################
#### Find the allele frequencies for the top SNPs  #############################################################################################################################################
################################################################################################################################################################################################

# Load phenotype file for plotting and annotating with frequencies
phenotypes <- read.csv(Phenofile, row.names = "X")

lt_phylo <- read.nexus("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/05_BEAST_RAxML.tre")

if (grepl("day_night_dif", pheno_to_test)) {
  pheno_to_test2 <- "day_night_dif"
} 
if (grepl("peak", pheno_to_test)) {
  pheno_to_test2 <- "peak"
} else {
  pheno_to_test2 <- pheno_to_test
}

# i <- grep("LOC100711869", genes_for_plotting$gene)

plots <- list()
for (i in 1:nrow(genes_for_plotting)) {
  
  # Subset the phylo tree
  subset <- keep.tip(lt_phylo, tip = phenotypes$six_letter_name_Ronco[phenotypes$six_letter_name_Ronco %in% lt_phylo$tip.label])
  d <- fortify(subset)
  
  # modify the SNP data
  phenotypes$tips <- phenotypes$six_letter_name_Ronco
  phenotypes$x <- 11
  phenotypes <- phenotypes[phenotypes$tips %in% subset$tip.label,]
  phenotypes$tips <- factor(phenotypes$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
  phenotypes$y <- as.numeric(phenotypes$tips)
  
  # extract data from genes_for_plotting
  gene <- genes_for_plotting$gene[i]
  chromosome <- genes_for_plotting$CHROM[i]
  coordinate <- genes_for_plotting$POS[i]
  
  # pick correct allefreqs file
  af <- allele_freq[[chromosome]]
  
  # generat col_name
  col_name <- paste("genotype", genes_for_plotting$gene[i], genes_for_plotting$location[i], sep = "_")
  col_name <- str_replace(col_name, ":", "_")
  col_name <- str_replace(col_name, "-", "_")
  
  # add genotype to phenotype file
  phenotypes[,col_name] <- t(af[af$coord == coordinate,5:ncol(af)])[match(phenotypes$six_letter_name_Ronco, row.names(t(af[af$coord == coordinate,5:ncol(af)])))]
  
  # round genotype frequencies for plotting
  phenotypes$Genotype <- as.character(ifelse(phenotypes[,col_name] < 0.15, 0, ifelse(phenotypes[,col_name] > 0.15 & phenotypes[,col_name] < 0.35, 0.25, ifelse(phenotypes[,col_name] > 0.35 & phenotypes[,col_name] < 0.65, 0.5, ifelse(phenotypes[,col_name] > 0.65 & phenotypes[,col_name] < 0.85, 0.75, 1.0)))))
  phenotypes <- phenotypes[!(is.na(phenotypes[,col_name])),]
  
  lm <- lm(phenotypes[,pheno_to_test2] ~ phenotypes[,col_name])
  direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
  
  ref_allele <- as.character(af[af$coord == coordinate,3])
  alt_allele <- as.character(af[af$coord == coordinate,4])
  
  genes_for_plotting$ref_allele[i] <- ref_allele
  genes_for_plotting$alt_allele[i] <- alt_allele
  
  # Make phenotype plots
  phylo.plot <- ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 3)
  
  if(grepl("day_night_dif", pheno_to_test)) {
    p3 <- phylo.plot + geom_tile(data = phenotypes, aes(y=y, x=x, fill = day_night_dif), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu", limits = c(-0.75,0.75), name = "Day/Night Diff. (rest)")
  }
  if(grepl("day_night_dif_spd", pheno_to_test)) {
    p3 <- phylo.plot + geom_tile(data = phenotypes, aes(y=y, x=x, fill = day_night_dif_spd), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu", limits = c(-50,50), name = "Day/Night Diff. (speed)")
  }
  if(grepl("peak", pheno_to_test)) {
    p3 <- phylo.plot + geom_tile(data = phenotypes, aes_string(y="y", x="x", fill = "peak"), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "BrBG", direction = 1, name = "Peak %")
  }
  if(grepl("total_rest", pheno_to_test)) {
    p3 <- phylo.plot + geom_tile(data = phenotypes, aes_string(y="y", x="x", fill = total_rest), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "PRGn", direction = -1, name = "Total rest (hrs)")
  }
  
  p3 <- p3 + new_scale("size") + new_scale("fill") + geom_tile(data = phenotypes, aes_string(y="y", x="x+1", fill = col_name)) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
  # p3 <- p3 + ggtitle(paste(pheno_to_test, " vs ", col_name, sep = ""))
  # p3 <- p3 + xlim(c(0,14))
  
  # Make boxplots
  
  if(grepl("day_night_dif", pheno_to_test)) {
    lm <- lm(phenotypes[,"day_night_dif"] ~ phenotypes[,col_name])
    direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
    p4 <- ggplot(phenotypes, aes(y=day_night_dif, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "RdBu", direction = direction)
  }
  if(grepl("day_night_dif_spd", pheno_to_test)) {
    lm <- lm(phenotypes[,"day_night_dif"] ~ phenotypes[,col_name])
    direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
    p4 <- ggplot(phenotypes, aes(y=day_night_dif_spd, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "RdBu", direction = direction)
  }
  if(grepl("peak", pheno_to_test)) {
    lm <- lm(phenotypes[,"peak"] ~ phenotypes[,col_name])
    direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
    p4 <- ggplot(phenotypes, aes(y=peak, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "BrBG", direction = direction)
  }
  if(grepl("total_rest", pheno_to_test)) {
    lm <- lm(phenotypes[,"total_rest"] ~ phenotypes[,col_name])
    direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
    p4 <- ggplot(phenotypes, aes(y=total_rest, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "PRGn", direction = direction)
  }
  
  # p4 <- p4 + ggtitle(paste(pheno_to_test, " vs ", col_name, sep = ""))
  p5 <- p3 + p4 + plot_layout(nrow = 2, guides = 'collect') + plot_annotation(title = paste(pheno_to_test, " vs ", col_name, sep = ""), subtitle = paste("GWAS p-value", genes_for_plotting$ps[i], "PGLS p-value", genes_for_plotting$PGLSps[i], sep = ", "))
  ggsave(p5, file = paste("../../outs/", pheno_to_test, "_PhyloPlots", "/", pheno_to_test, "_vs_", col_name, ".pdf", sep = ""))
  
}


# phenotypes is also appended with the genotype for each of the genes so we can save that
write.csv(phenotypes, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/", pheno_to_test, "_PhyloPlots", "/Pheno_vs_Geno_", pheno_to_test, "_top", cutoff, "allele_freqs_topSNPs.csv", sep = ""))

# Write out again the genes, this time with the appended alleles
write.csv(genes_for_plotting, file = paste("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs/", pheno_to_test, "_PhyloPlots", "/GWAS_vs_PGLS_", pheno_to_test, "_top", cutoff, "_genes_for_plotting.csv", sep = ""))


## If making for presentation

p3
p4 <- p4 + ylab("Total rest (hrs)") + xlab("Genotype allele frequency")

p5 <- p3 + p4 + plot_layout(nrow = 2, guides = 'collect') + plot_annotation(title = paste(pheno_to_test, " vs ", col_name, sep = ""), subtitle = paste("GWAS p-value", genes_for_plotting$ps[i], "PGLS p-value", genes_for_plotting$PGLSps[i], sep = ", "))


# # For each of the comparisons, load the topSNP allele freqs, and plot them against the corresponding phenotype
# # Then save pdfs for each in a subfolder
# 
# comparisons <- list.files()
# comparisons <- comparisons[grepl("allele_freqs_topSNPs.csv", comparisons)]
# 
# comparison.names <- c("day_night_dif", "day_night_dif_55-species", "day_night_dif_58-species", "day_night_dif_spd_55-species", "peak", "peak_amplitude", "total_rest")
# 
# for (i in 1:length(comparisons)) {
#   # load in the snp data
#   snp.data <- read.csv(comparisons[i])
#   snp.data <- snp.data[,c(1:8,10:ncol(snp.data))]
#   last.col <- ncol(snp.data)
#   
#   # Subset the phylo tree
#   subset <- keep.tip(lt_phylo, tip = snp.data$six_letter_name_Ronco[snp.data$six_letter_name_Ronco %in% lt_phylo$tip.label])
#   d <- fortify(subset)
#   
#   # modify the SNP data
#   snp.data$tips <- snp.data$six_letter_name_Ronco
#   snp.data$x <- 11
#   snp.data <- snp.data[snp.data$tips %in% subset$tip.label,]
#   snp.data$tips <- factor(snp.data$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
#   snp.data$y <- as.numeric(snp.data$tips)
#   
#   #plots <- list()
#   for (j in 8:last.col) {
#     phylo.plot <- ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 3)
#     
#     if(grepl("day_night_dif", comparison.names[i])) {
#       p3 <- phylo.plot + geom_tile(data = snp.data, aes(y=y, x=x, fill = day_night_dif), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu")
#     }
#     if(grepl("peak", comparison.names[i])) {
#       p3 <- phylo.plot + geom_tile(data = snp.data, aes_string(y="y", x="x", fill = comparison.names[i]), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "BrBG", direction = 1)
#     }
#     if(grepl("total_rest", comparison.names[i])) {
#       p3 <- phylo.plot + geom_tile(data = snp.data, aes_string(y="y", x="x", fill = comparison.names[i]), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "PRGn", direction = -1)
#     }
#     
#     p3 <- p3 + new_scale("size") + new_scale("fill") + geom_tile(data = snp.data, aes_string(y="y", x="x+1", fill = colnames(snp.data)[j])) + scale_fill_distiller(palette = "Greys", direction = 1) 
#     setwd
#     p3 <- p3 + ggtitle(paste(comparison.names[i], " vs ", colnames(snp.data)[j], sep = ""))
#     p3 <- p3 + xlim(c(0,14))
#     
#     # Make boxplots
#     
#     snp.data$Genotype <- as.character(ifelse(snp.data[,colnames(snp.data)[j]] < 0.15, 0, ifelse(snp.data[,colnames(snp.data)[j]] > 0.15 & snp.data[,colnames(snp.data)[j]] < 0.35, 0.25, ifelse(snp.data[,colnames(snp.data)[j]] > 0.35 & snp.data[,colnames(snp.data)[j]] < 0.65, 0.5, ifelse(snp.data[,colnames(snp.data)[j]] > 0.65 & snp.data[,colnames(snp.data)[j]] < 0.85, 0.75, 1.0)))))
#     snp.data <- snp.data[!(is.na(snp.data[,colnames(snp.data)[j]])),]
#     
#     if(grepl("day_night_dif", comparison.names[i])) {
#       lm <- lm(snp.data[,"day_night_dif"] ~ snp.data[,colnames(snp.data)[j]])
#       direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
#       p4 <- ggplot(snp.data, aes(y=day_night_dif, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "RdBu", direction = direction)
#     }
#     if(grepl("peak", comparison.names[i])) {
#       lm <- lm(snp.data[,"peak"] ~ snp.data[,colnames(snp.data)[j]])
#       direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
#       p4 <- ggplot(snp.data, aes(y=peak, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "BrBG", direction = direction)
#     }
#     if(grepl("total_rest", comparison.names[i])) {
#       lm <- lm(snp.data[,"total_rest"] ~ snp.data[,colnames(snp.data)[j]])
#       direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
#       p4 <- ggplot(snp.data, aes(y=total_rest, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "PRGn", direction = direction)
#     }
#     
#     p4 <- p4 + ggtitle(paste(comparison.names[i], " vs ", colnames(snp.data)[j], sep = ""))
#     p5 <- p3 + p4 + plot_layout(nrow = 1)
#     ggsave(p5, file = paste(comparison.names[i], "_PhyloPlots", "/", comparison.names[i], "_vs_", colnames(snp.data)[j], ".pdf", sep = ""))
#   }
# }