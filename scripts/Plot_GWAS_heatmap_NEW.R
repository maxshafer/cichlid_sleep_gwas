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
phenotypes <- read.csv("pheno_data/GWAS_data_PCs_TR_2024-01-28.csv", row.names = "X")

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

## Extract allelefreqs for all significant snps

extracted_af <- lapply(allele_freq, function(x) {
  x$location <- paste(x$chr, x$coord, sep = ":")
  x.pc1 <- x[x$location %in% snps$pc1$location]
  x.pc2 <- x[x$location %in% snps$pc2$location]
  x.tr <- x[x$location %in% snps$tr$location]
  x <- list(x.pc1, x.pc2, x.tr)
  names(x) <- c("pc1", "pc2", "tr")
  return(x)
})

extracted_af <- extracted_af[2:23]
pc1_alleles <- Reduce(rbind, lapply(extracted_af, function(x) x[[1]]))
pc2_alleles <- Reduce(rbind, lapply(extracted_af, function(x) x[[2]]))
tr_alleles <- Reduce(rbind, lapply(extracted_af, function(x) x[[3]]))

## Convert allele frequencies into Nocturnal (1) - Diurnal (0) based on lm result

transformAlleleFreqs <- function(alleles = pc1_alleles, pheno_file = phenotypes) {
  
  alleles_trans <- t(alleles[,5:65])
  
  test_convert <- apply(alleles_trans, 2, function(x) {
    x <- x[names(x) %in% phenotypes$six_letter_name_Ronco]
    lm <- lm(phenotypes[,"pc1"] ~ x)
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

pc1_convert <- transformAlleleFreqs(alleles = pc1_alleles, pheno_file = phenotypes)
pc2_convert <- transformAlleleFreqs(alleles = pc2_alleles, pheno_file = phenotypes)
tr_convert <- transformAlleleFreqs(alleles = tr_alleles, pheno_file = phenotypes)


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
  
  if (pheno_plot == "pc1") {
    dendro_phylo_pheno <- dendro_phylo_pheno + scale_fill_distiller(palette = "RdBu", direction = 1)
    dendro_snps_pheno <- dendro_snps_pheno + scale_fill_distiller(palette = "RdBu", direction = 1)
  }
  if (pheno_plot == "pc2") {
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


heatmap_pc1_phylo <- returnHeatMap(af_data_convert = pc1_convert, order_by = "phylo", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pc1")
heatmap_pc1_snps <- returnHeatMap(af_data_convert = pc1_convert, order_by = "snps", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pc1")
## OK for above, it really separates by lamprologs, ectos (and allies) and trophs, where the lamprologs use a much larger # of SNPS (more heterogeniety)
## Neobue has 0 nocturnal alleles (there is a group of alleles that are nocturnal, that all other species have)
## Neomul and Neobre have all nocturnal alleles

heatmap_pc2_phylo <- returnHeatMap(af_data_convert = pc2_convert, order_by = "phylo", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pc2")
heatmap_pc2_snps <- returnHeatMap(af_data_convert = pc2_convert, order_by = "snps", phylo_tree = lt_phylo, phenotype_data = phenotypes, pheno_plot = "pc2")
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

pdf("outs/heatmap_pc1_phylo_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_pc1_phylo
dev.off()

png("outs/heatmap_pc1_phylo_perchr_1e-05_percentile.png", width = 15, height = 10, units = "in", res = 500)
heatmap_pc1_phylo
dev.off()

pdf("outs/heatmap_pc1_snps_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_pc1_snps
dev.off()

pdf("outs/heatmap_pc2_phylo_perchr_1e-05_percentile.pdf", width = 15, height = 10)
heatmap_pc2_phylo
dev.off()

pdf("outs/heatmap_pc2_snps_perchr_1e-05_percentile.pdf", width = 15, height = 10)
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

plot_combined <- dendro_phylo_pheno + plot + plot_spacer() + dendro_plot_r + plot_layout(design = design, guides = "collect", width = unit(c(200,800), "mm"), height = unit(c(600,100), "mm"))


pdf("outs/heatmap_pc1_phylo_labelled_perchr_1e-05_percentile.pdf", width = 100, height = 10)
plot
dev.off()




##### Can I modify the function to return a plot with the SNPs of interest also shown?

subset <- keep.tip(phylo_tree, tip = phenotypes$six_letter_name_Ronco[phenotypes$six_letter_name_Ronco %in% phylo_tree$tip.label])
dendro_phylo <- ggtree(subset) #+ geom_tiplab(offset = 13)

new_pheno <- as.matrix(pc1_convert[,-c(1:5)])
rownames(new_pheno) <- pc1_convert$location
new_pheno <- t(new_pheno)

pheno_data <- phenotypes[match(rownames(new_pheno), phenotypes$six_letter_name_Ronco),]
pheno_data$y_phy <- dendro_phylo$data$y[match(pheno_data$six_letter_name_Ronco, dendro_phylo$data$label)]
pheno_data$x_phy <- dendro_phylo$data$x[match(pheno_data$six_letter_name_Ronco, dendro_phylo$data$label)]

pheno_data$SNP1 <- new_pheno[,grep("NC_031969.2:31739582", colnames(new_pheno))] # NC_031969.2:31739582 Lamps and ectos
pheno_data$SNP2 <- new_pheno[,grep("NC_031976.2:30604722", colnames(new_pheno))] # NC_031978.2:19381306 NC_031976.2:30604722, NC_031979.2:1797468 (just lampros)
pheno_data$SNP3 <- new_pheno[,grep("NC_031972.2:22216041", colnames(new_pheno))] # NC_031965.2:37067101, NC_031978.2:10013506, NC_031972.2:22216040, NC_031972.2:22216041, NC_031966.2:13031943 (just ectos)
pheno_data$SNP4 <- new_pheno[,grep("NC_031969.2:28463972", colnames(new_pheno))] # NC_031987.2:1114119, NC_031973.2:22078731, NC_031973.2:22078731, NC_031969.2:28463972 (just buescheri)

dendro_phylo_pheno <- dendro_phylo + geom_tile(data = pheno_data, aes_string(y="y_phy", x="x_phy+1.5", fill = pheno_plot), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", direction = 1)

dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+4, fill = SNP3), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+6.5, fill = SNP2), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+9, fill = SNP1), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 
dendro_phylo_pheno <- dendro_phylo_pheno + new_scale("size") + new_scale("fill") + geom_tile(data = pheno_data, aes(y=y_phy, x=x_phy+11.5, fill = SNP4), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "Greys", direction = 1, name = "Genotype Freq.") 

#dendro_phylo_pheno

plot <- ggplot(phenotypes_dnd_long, aes(x=location, y = measurement, fill = value)) + geom_tile() + theme(axis.title = element_blank(), axis.text.x = element_blank()) + scale_fill_distiller(palette = "Greys", direction = 1) + guides(fill=guide_legend(title= paste(pheno_plot, " allele freq", sep = "") ))

figure_first <- dendro_phylo_pheno + plot + plot_spacer() + dendro_plot_r + plot_layout(design = design, guides = "collect", width = unit(c(35,150), "mm"), height = unit(c(50,10), "mm"))


pdf("outs/first_figure_test.pdf", width = 10, height = 5)
figure_first
dev.off()




