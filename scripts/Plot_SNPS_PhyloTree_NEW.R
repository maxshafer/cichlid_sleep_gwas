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
library(here)
library(ggnewscale)

# Source this gist (From "https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32"), which allows using multiple of the same aes for the same plot
# This is now a library (ggnewscale), which works as intended.

# source('/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R')

## Load in the snps, and subset to unique gene/location matches
filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))
snps <- filter_snps[2:4]

snps <- lapply(snps, function(x) {
  x <- x %>% distinct(ANN_GENE, location, .keep_all = T)
  return(x)
})
names(snps) <- c("pc1", "pc2", "tr")

# Preload all allelefreqs so I can use it also to make plots for the exonic genes (or at least fetch the REF/ALT nucleotides)
if (!(exists("allele_freq"))) {
  af_files <- list.files("sra_reads_nobackup/alleleFreqs/")
  chromosomes <- paste(str_split_i(af_files, "_", 9), str_split_i(af_files, "_", 10), sep = "_")
  allele_freq <- lapply(af_files, function(x) fread(paste("sra_reads_nobackup/alleleFreqs/",x, sep = "")))
  names(allele_freq) <- chromosomes
}

# Load phenotype file for plotting and annotating with frequencies
phenotypes <- read.csv("pheno_data/GWAS_data_PCs_TR_2024-01-28.csv", row.names = "X")

# Load the tree
lt_phylo <- read.nexus("scripts/05_BEAST_RAxML.tre")

################################################################################################################################################################################################
#### Find the allele frequencies for the top SNPs  #############################################################################################################################################
################################################################################################################################################################################################

genes_for_plotting <- snps[[1]]
pheno_to_test2 <- "pc1"

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
  gene <- genes_for_plotting$ANN_GENE[i]
  chromosome <- genes_for_plotting$CHROM[i]
  coordinate <- genes_for_plotting$POS[i]
  
  # pick correct allefreqs file
  af <- allele_freq[[chromosome]]
  
  # generat col_name
  col_name <- paste("genotype", genes_for_plotting$ANN_GENE[i], genes_for_plotting$location[i], sep = "_")
  col_name <- str_replace(col_name, ":", "_")
  col_name <- str_replace(col_name, "-", "_")
  col_name <- str_replace(col_name, "-", "_")
  
  # add genotype to phenotype file
  phenotypes[,col_name] <- t(af[af$coord == coordinate,5:ncol(af)])[match(phenotypes$six_letter_name_Ronco, row.names(t(af[af$coord == coordinate,5:ncol(af)])))]
  
  # round genotype frequencies for plotting
  phenotypes$Genotype <- as.character(ifelse(phenotypes[,col_name] < 0.15, 0, ifelse(phenotypes[,col_name] > 0.15 & phenotypes[,col_name] < 0.35, 0.25, ifelse(phenotypes[,col_name] > 0.35 & phenotypes[,col_name] < 0.65, 0.5, ifelse(phenotypes[,col_name] > 0.65 & phenotypes[,col_name] < 0.85, 0.75, 1.0)))))
  phenotypes <- phenotypes[!(is.na(phenotypes[,col_name])),]
  
  # lm <- lm(phenotypes[,pheno_to_test2] ~ phenotypes[,col_name])
  # direction <- ifelse(lm$coefficients[2] > 0, 1, -1)
  
  ref_allele <- as.character(af[af$coord == coordinate,3])
  alt_allele <- as.character(af[af$coord == coordinate,4])
  
  # genes_for_plotting$ref_allele[i] <- ref_allele
  # genes_for_plotting$alt_allele[i] <- alt_allele
  
  # Make phenotype plots
  phylo.plot <- ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 3)
  
  if(grepl("pc1", pheno_to_test2)) {
    p3 <- phylo.plot + geom_tile(data = phenotypes, aes(y=y, x=x, fill = pc1), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu", name = "PC1 Value", direction = 1)
  }
  if(grepl("pc2", pheno_to_test2)) {
    p3 <- phylo.plot + geom_tile(data = phenotypes, aes_string(y="y", x="x", fill = "pc2"), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "BrBG", direction = -1, name = "PC2 Value")
  }
  if(grepl("tr", pheno_to_test2)) {
    p3 <- phylo.plot + geom_tile(data = phenotypes, aes_string(y="y", x="x", fill = "total_rest"), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "PRGn", direction = -1, name = "Total rest (hrs)")
  }
  
  p3 <- p3 + new_scale("size") + new_scale("fill") + geom_tile(data = phenotypes, aes_string(y="y", x="x+1", fill = as.character(col_name))) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
  
  # Make boxplots
  
  if(grepl("pc1", pheno_to_test2)) {
    lm <- lm(phenotypes[,"pc1"] ~ phenotypes[,col_name])
    direction <- ifelse(lm$coefficients[2] > 0, 1, -1)
    p4 <- ggplot(phenotypes, aes(y=pc1, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "RdBu", direction = direction)
  }
  if(grepl("pc2", pheno_to_test2)) {
    lm <- lm(phenotypes[,"pc2"] ~ phenotypes[,col_name])
    direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
    p4 <- ggplot(phenotypes, aes(y=pc2, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "BrBG", direction = direction)
  }
  if(grepl("tr", pheno_to_test2)) {
    lm <- lm(phenotypes[,"total_rest"] ~ phenotypes[,col_name])
    direction <- ifelse(lm$coefficients[2] > 0, -1, 1)
    p4 <- ggplot(phenotypes, aes(y=total_rest, x=Genotype, group=Genotype, fill=Genotype)) + geom_boxplot(outlier.colour = "transparent") + geom_point(position = "jitter") + theme_classic() + scale_fill_brewer(palette = "PRGn", direction = direction)
  }
  
  # p4 <- p4 + ggtitle(paste(pheno_to_test, " vs ", col_name, sep = ""))
  p5 <- p3 + p4 + plot_layout(nrow = 2, guides = 'collect') 
  p5 <- p5 + plot_annotation(title = paste(pheno_to_test2, " vs ", as.character(col_name), sep = ""))
  p5 <- p5 + plot_annotation(subtitle = paste("GWAS p-value ", 
                                              format(genes_for_plotting$piVals[i], scientific = T, digits = 2), 
                                              ", PGLS p-value ", 
                                              format(genes_for_plotting$PGLSpiVals[i], scientific = T, digits = 2), 
                                              "\nREF: ", ref_allele, ", ALT: ", alt_allele, 
                                              ", ", genes_for_plotting$ANN_IMPACT[i], " / ", genes_for_plotting$ANN_EFFECT[i], " / ", genes_for_plotting$ANN_DISTANCE[i], " bp", sep = ""))
  
  ggsave(p5, file = paste("outs/", pheno_to_test2, "_PhyloPlots", "/", pheno_to_test2, "_vs_", col_name, ".pdf", sep = ""))
  
}




# phenotypes is also appended with the genotype for each of the genes so we can save that
write.csv(phenotypes, file = paste("outs/", pheno_to_test2, "_PhyloPlots", "/Pheno_vs_Geno_", pheno_to_test2, "_perchr_1e-05_allele_freqs.csv", sep = ""))

