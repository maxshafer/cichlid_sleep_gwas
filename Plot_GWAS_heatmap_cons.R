library(data.table)
library(dplyr)
library(readr)
library(patchwork)
library(here)
library(ggnewscale)
library(stringr)
library(ggtree)
library(tidyr)
library(ape)


# Load the phylogeny

lt_phylo <- read.nexus("scripts/05_BEAST_RAxML.tre")

# Load phenotype file for plotting and annotating with frequencies
phenotypes <- read.csv("pheno_data/cichlids_consolidation_values_full.csv")

## Load in the snps, and subset to unique gene/location matches
filter_snps <- readRDS(here("filter_SNPs_perchr_1e-05_percentile"))
snps <- filter_snps[2:3]

snps <- lapply(snps, function(x) {
  x <- x %>% distinct(ANN_GENE, location, .keep_all = T)
  return(x)
})
names(snps) <- c("sfi", "cons")

# Preload all allelefreqs so I can use it also to make plots for the exonic genes (or at least fetch the REF/ALT nucleotides)
if (!(exists("allele_freq"))) {
  af_files <- list.files("../alleleFreqs/")
  chromosomes <- paste(str_split_i(af_files, "_", 9), str_split_i(af_files, "_", 10), sep = "_")
  allele_freq <- lapply(af_files, function(x) fread(paste("../alleleFreqs/",x, sep = "")))
  names(allele_freq) <- chromosomes
}

## Extract allelefreqs for all significant snps

extracted_af <- lapply(allele_freq, function(x) {
  x$location <- paste(x$chr, x$coord, sep = ":")
  x.sfi <- x[x$location %in% snps$sfi$location]
  x.cons <- x[x$location %in% snps$cons$location]
  x <- list(x.sfi, x.cons)
  names(x) <- c("sfi", "cons")
  return(x)
})

extracted_af <- extracted_af[2:23]
sfi_alleles <- Reduce(rbind, lapply(extracted_af, function(x) x[[1]]))
cons_alleles <- Reduce(rbind, lapply(extracted_af, function(x) x[[2]]))

## Convert allele frequencies into Nocturnal (1) - Diurnal (0) based on lm result

transformAlleleFreqs <- function(alleles = sfi_alleles, pheno_file = phenotypes, pheno = "pref_mean_sfi") {
  
  alleles_trans <- t(alleles[,5:65])
  
  test_convert <- apply(alleles_trans, 2, function(x) {
    x <- x[names(x) %in% phenotypes$six_letter_name_Ronco]
    lm <- lm(phenotypes[,pheno] ~ x)
    direction <- ifelse(lm$coefficients[2] > 0, 1, -1)
    
    if(direction == 1) {
      x <- x
    }
    if(direction == -1) {
      x <- 1-x
    }
    return(x)
  })
  
  alleles_convert <- t(test_convert)
  alleles_convert <- cbind(alleles[,c(1:4,66)], alleles_convert)
  return(alleles_convert)
}

sfi_convert <- transformAlleleFreqs(alleles = sfi_alleles, pheno_file = phenotypes, pheno = "pref_sfi")
cons_convert <- transformAlleleFreqs(alleles = cons_alleles, pheno_file = phenotypes, pheno = "pref_max_cons")


###### ###### ###### ###### ###### ###### ###### ###### ###### 
##### Make dendrogram for species
###### ###### ###### ###### ###### ###### ###### ###### ###### 

## Need to modify so that it colours things correctly
## And double check that it works properly for other phenos

returnHeatMap <- function(af_data_convert = pc1_convert, order_by = c("snps", "phylo"), phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pc1") {
  require(ape)
  require(ggtree)
  require(ggplot2)
  #########  #########  #########  #########
  ## Make dendrogram based on phylogeny   ##
  #########  #########  #########  #########
  
  subset <- keep.tip(phylo_tree, tip = phenotype_data$six_letter_name_Ronco[phenotype_data$six_letter_name_Ronco %in% phylo_tree$tip.label])
  dendro_phylo <- ggtree(subset) #+ geom_tiplab(offset = 1.25)
  phenotype_data$y_phy <- dendro_phylo$data$y[match(phenotype_data$six_letter_name_Ronco, dendro_phylo$data$label)]
  phenotype_data$x_phy <- dendro_phylo$data$x[match(phenotype_data$six_letter_name_Ronco, dendro_phylo$data$label)]
  dendro_phylo_pheno <- dendro_phylo + geom_tile(data = phenotype_data, aes_string(y="y_phy", x="x_phy+2.5", fill = pheno_plot), width = 4, height = 0.9, size = 1.5, inherit.aes = FALSE)
  
  
  #########  #########  #########  #########
  ## Run clustering by SNPS               ##
  #########  #########  #########  #########
  
  phenotypes_matrix <- as.matrix(t(af_data_convert[, -c(1:5)]))
  dendro_snps <- ggtree(as.phylo(hclust(d = dist(x = phenotypes_matrix)))) #+ geom_tiplab(offset = 1.25) 
  phenotype_data$y <- dendro_snps$data$y[match(phenotype_data$six_letter_name_Ronco, dendro_snps$data$label)]
  phenotype_data$x <- dendro_snps$data$x[match(phenotype_data$six_letter_name_Ronco, dendro_snps$data$label)]
  dendro_snps_pheno <- dendro_snps + geom_tile(data = phenotype_data, aes_string(y="y", x="x_phy+4.25", fill = pheno_plot), width = 4, height = 0.9, size = 1.5, inherit.aes = FALSE)
  
  ## Add phenotypes to plots
  
  if (pheno_plot == "mean_sfi") {
    dendro_phylo_pheno <- dendro_phylo_pheno + scale_fill_distiller(palette = "RdBu", direction = 1)
    dendro_snps_pheno <- dendro_snps_pheno + scale_fill_distiller(palette = "RdBu", direction = 1)
  }
  if (pheno_plot == "mean_max_cons") {
    dendro_phylo_pheno <- dendro_phylo_pheno + scale_fill_distiller(palette = "BrBG", direction = -1)
    dendro_snps_pheno <- dendro_snps_pheno + scale_fill_distiller(palette = "BrBG", direction = -1)
  }
  if (pheno_plot == "total_rest") {
    dendro_phylo_pheno <- dendro_phylo_pheno + scale_fill_distiller(palette = "PRGn", direction = -1)
    dendro_snps_pheno <- dendro_snps_pheno + scale_fill_distiller(palette = "PRGn", direction = -1)
  }
  
  
  ###### ###### ###### ###### ###### ###### ###### 
  ## Run clustering of species using SNPs       ##
  ###### ###### ###### ###### ###### ###### ######  
  
  # Run clustering
  phenotypes_matrix_r <- as.matrix(af_data_convert[,-c(1:5)])
  rownames(phenotypes_matrix_r) <- af_data_convert$location
  
  # Create dendro
  dendro_plot_r <- ggtree(as.phylo(hclust(d = dist(x = phenotypes_matrix_r)))) + coord_flip()
  phenotypes_scaled_r_x_order <- dendro_plot_r$data$label[order(dendro_plot_r$data$y)]
  phenotypes_scaled_r_x_order <- phenotypes_scaled_r_x_order[!(is.na(phenotypes_scaled_r_x_order))]
  
  
  ###### ###### ###### ###### ###### ######
  ## Wrangle data                        ##
  ###### ###### ###### ###### ###### ###### 
  
  rownames(af_data_convert) <- af_data_convert$location
  
  phenotypes_dnd_long <- pivot_longer(data = af_data_convert,
                                      cols = -c(chr, coord, ref, alt, location),
                                      names_to = "measurement",
                                      values_to = "value")
  
  ###### ###### ###### ###### ###### ######
  ## Plot Heatmaps                       ##
  ###### ###### ###### ###### ###### ######
  
  design = "
  AB
  CD"
  
  if (order_by == "snps") {
    
    ## Change factor levels to modify ordering of X and Y
    phenotypes_scaled_y_order <- dendro_snps$data$label[order(dendro_snps$data$y)]
    phenotypes_scaled_y_order <- phenotypes_scaled_y_order[!(is.na(phenotypes_scaled_y_order))]
    phenotypes_dnd_long$measurement <- factor(phenotypes_dnd_long$measurement, levels = phenotypes_scaled_y_order)
    phenotypes_dnd_long$location <- factor(phenotypes_dnd_long$location, levels = phenotypes_scaled_r_x_order)
    
    ## Plot heatmap
    plot <- ggplot(phenotypes_dnd_long, aes(x=location, y = measurement, fill = value)) + geom_tile() + theme(axis.title = element_blank(), axis.text.x = element_blank()) + scale_fill_distiller(palette = "Greys", direction = 1) + guides(fill=guide_legend(title= paste(pheno_plot, " allele freq", sep = "") ))
    
    plot_combined <- dendro_snps_pheno + plot + plot_spacer() + dendro_plot_r + plot_layout(design = design, guides = "collect", width = unit(c(50,200), "mm"), height = unit(c(150,25), "mm"))
  }
  
  if (order_by == "phylo") {
    ## Just order by phylogeny and dendrogram of snps
    phenotypes_dnd_long$measurement <- factor(phenotypes_dnd_long$measurement, levels = dendro_phylo_pheno$data$label[order(dendro_phylo_pheno$data$y)])
    phenotypes_dnd_long$location <- factor(phenotypes_dnd_long$location, levels = phenotypes_scaled_r_x_order)
    
    ## Plot heatmap
    plot <- ggplot(phenotypes_dnd_long, aes(x=location, y = measurement, fill = value)) + geom_tile() + theme(axis.title = element_blank(), axis.text.x = element_blank()) + scale_fill_distiller(palette = "Greys", direction = 1) + guides(fill=guide_legend(title= paste(pheno_plot, " allele freq", sep = "") ))
    
    plot_combined <- dendro_phylo_pheno + plot + plot_spacer() + dendro_plot_r + plot_layout(design = design, guides = "collect", width = unit(c(50,200), "mm"), height = unit(c(150,25), "mm"))
  }
  
  return(plot_combined)
  
}


heatmap_sfi_phylo <- returnHeatMap(af_data_convert = sfi_convert, order_by = "phylo", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pref_sfi")
heatmap_sfi_snps <- returnHeatMap(af_data_convert = sfi_convert, order_by = "snps", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pref_sfi")
## OK for above, it really separates by lamprologs, ectos (and allies) and trophs, where the lamprologs use a much larger # of SNPS (more heterogeniety)
## Neobue has 0 nocturnal alleles (there is a group of alleles that are nocturnal, that all other species have)
## Neomul and Neobre have all nocturnal alleles

heatmap_cons_phylo <- returnHeatMap(af_data_convert = cons_convert, order_by = "phylo", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pref_max_cons")
heatmap_cons_snps <- returnHeatMap(af_data_convert = cons_convert, order_by = "snps", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pref_max_cons")
## This is the least clear
## One large group of SNPs is mostly found in high pc2, and not low
## Another group is sort of lost in some that have low pc2 values
## My guess is that this is because dawn/dusk peaks are different genes (at least in lamprologs)


heatmap_tr_phylo <- returnHeatMap(af_data_convert = tr_convert, order_by = "phylo", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "total_rest")
heatmap_tr_snps <- returnHeatMap(af_data_convert = tr_convert, order_by = "snps", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "total_rest")
## For total_rest, it also makes some sense
## There are three groups:
## Lamprologs are split in two, with one group using a large swatch of linked SNPS, and the other using some of those, but also some novel ones
## The third group is the ectos, which use different subsets of teh large linked swath


### Save out the plots

pdf("plots/heatmap_sfi_phylo_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_pc1_phylo
dev.off()

png("plots/heatmap_sfi_phylo_perchr_1e-05_percentile.png", width = 15, height = 10, units = "in", res = 500)
heatmap_sfi_phylo
dev.off()

pdf("plots/heatmap_sfi_snps_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_pc1_snps
dev.off()

pdf("plots/heatmap_cons_phylo_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_pc2_phylo
dev.off()

pdf("outs/heatmap_cons_snps_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_pc2_snps
dev.off()

pdf("outs/heatmap_tr_phylo_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_tr_phylo
dev.off()

pdf("outs/heatmap_tr_snps_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_tr_snps
dev.off()



### Below is for finding example SNPs for phyloplots
### First export a large version of the heatmap with snp labels
### Use the groupings to find a label, and then make a PhyloPlot with it

plot <- ggplot(phenotypes_dnd_long, aes(x=location, y = measurement, fill = value)) + geom_tile() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_distiller(palette = "Greys", direction = 1) + guides(fill=guide_legend(title= paste(pheno_plot, " allele freq", sep = "") ))

plot_combined <- dendro_phylo_pheno + plot + plot_layout(design = design, guides = "collect", width = unit(c(25,75), "in"), height = unit(c(10), "in"))


pdf("outs/heatmap_pc1_phylo_labelled_perchr_1e-05_percentile.pdf", width = 100, height = 10)
plot
dev.off()




##### Can I modify the function to return a plot with the SNPs of interest also shown?

subset <- keep.tip(lt_phylo, tip = phenotypes$six_letter_name_Ronco[phenotypes$six_letter_name_Ronco %in% lt_phylo$tip.label])
dendro_phylo <- ggtree(subset) #+ geom_tiplab(offset = 13)

new_pheno <- as.matrix(pc1_convert[,-c(1:5)])
rownames(new_pheno) <- pc1_convert$location
new_pheno <- t(new_pheno)

pheno_data <- phenotypes[match(rownames(new_pheno), phenotypes$six_letter_name_Ronco),]
pheno_data$y_phy <- dendro_phylo$data$y[match(pheno_data$six_letter_name_Ronco, dendro_phylo$data$label)]
pheno_data$x_phy <- dendro_phylo$data$x[match(pheno_data$six_letter_name_Ronco, dendro_phylo$data$label)]

pheno_data$SNP1 <- new_pheno[,grep("NC_031982.2:23856121", colnames(new_pheno))] # Group i) Ectos/limnos/cyphos + outgroup lamps
pheno_data$SNP2 <- new_pheno[,grep("NC_031971.2:24810009", colnames(new_pheno))] # Group ii) + iv) Just lamprologini NC_031974.2:29083768
pheno_data$SNP3 <- new_pheno[,grep("NC_031978.2:6846904", colnames(new_pheno))] # Group iii) Opposite pattern for strong dirunal lampros NC_031972.2:60671126
pheno_data$SNP4 <- new_pheno[,grep("NC_031975.2:18513666", colnames(new_pheno))] # Group v) all but the mostly diurnal tropheini NC_031975.2:18513666 (crcp)

dendro_phylo_pheno <- dendro_phylo + geom_tile(data = pheno_data, aes_string(y="y_phy", x="x_phy+1.5", fill = "pc1"), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", direction = 1)

dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+4, fill = SNP4), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+6.5, fill = SNP1), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+9, fill = SNP2), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+11.5, fill = SNP3), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 

#dendro_phylo_pheno

plot <- ggplot(phenotypes_dnd_long, aes(x=location, y = measurement, fill = value)) + geom_tile() + theme(axis.title = element_blank(), axis.text.x = element_blank()) + scale_fill_distiller(palette = "Greys", direction = 1) + guides(fill=guide_legend(title= paste(pheno_plot, " allele freq", sep = "") ))

figure_first <- dendro_phylo_pheno + plot + plot_spacer() + dendro_plot_r + plot_layout(design = design, guides = "collect", width = unit(c(35,150), "mm"), height = unit(c(50,10), "mm"))


pdf("outs/first_figure_test.pdf", width = 10, height = 5)
figure_first
dev.off()




