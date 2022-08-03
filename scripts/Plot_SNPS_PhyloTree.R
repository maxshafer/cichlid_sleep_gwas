library(ape)
library(ggtree)
library(ggplot2)
library(viridis)
library(stringr)

# Source this gist (From "https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32"), which allows using multiple of the same aes for the same plot
# Here I use it to make multiple heatmaps for Annika's figure (but should work better than the super hack I used for Figure 5 of the Cavefish single-cell paper)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/outs")
source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")

lt_phylo <- read.nexus("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/05_BEAST_RAxML.tre")


# For each of the comparisons, load the topSNP allele freqs, and plot them against the corresponding phenotype
# Then save pdfs for each in a subfolder

comparisons <- list.files()
comparisons <- comparisons[grepl("allele_freqs_topSNPs.csv", comparisons)]

comparison.names <- str_sub(comparisons, start = 14, end = -26)

for (i in 1:length(comparisons)) {
  # load in the snp data
  snp.data <- read.csv(comparisons[i])
  last.col <- ncol(snp.data)
  
  # Subset the phylo tree
  subset <- keep.tip(lt_phylo, tip = snp.data$six_letter_name_Ronco[snp.data$six_letter_name_Ronco %in% lt_phylo$tip.label])
  d <- fortify(subset)
  
  # modify the SNP data
  snp.data$tips <- snp.data$six_letter_name_Ronco
  snp.data$x <- 11
  snp.data <- snp.data[snp.data$tips %in% subset$tip.label,]
  snp.data$tips <- factor(snp.data$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
  snp.data$y <- as.numeric(snp.data$tips)
  
  #plots <- list()
  for (j in 8:last.col) {
    phylo.plot <- ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 3)
    
    if(grepl("day_night_dif", comparison.names[i])) {
      p3 <- phylo.plot + geom_tile(data = snp.data, aes(y=y, x=x, fill = day_night_dif), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu")
    }
    if(grepl("peak", comparison.names[i])) {
      p3 <- phylo.plot + geom_tile(data = snp.data, aes_string(y="y", x="x", fill = comparison.names[i]), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "BrBG", direction = 1)
    }
    if(grepl("total_rest", comparison.names[i])) {
      p3 <- phylo.plot + geom_tile(data = snp.data, aes_string(y="y", x="x", fill = comparison.names[i]), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "PRGn", direction = -1)
    }
    
    p3 <- p3 + new_scale("size") + new_scale("fill") + geom_tile(data = snp.data, aes_string(y="y", x="x+1", fill = colnames(snp.data)[j])) + scale_fill_distiller(palette = "Greys", direction = 1) 
    
    p3 <- p3 + ggtitle(paste(comparison.names[i], " vs ", colnames(snp.data)[j], sep = ""))
    p3 <- p3 + xlim(c(0,14))
    ggsave(p3, file = paste(comparison.names[i], "_PhyloPlots", "/", comparison.names[i], "_vs_", colnames(snp.data)[j], ".pdf", sep = ""))
    #plots[[j]] <- p3
  }
}


snp.data <- read.csv("outs/GWAS_vs_PGLS_day_night_dif_58-species_allele_freqs_topSNPs.csv")

for (j in 8:last.col) {
  pdf(file = paste(comparison.names[i], "_PhyloPlots", "/", comparison.names[i], "_vs_", colnames(snp.data)[j], ".pdf", sep = ""))
  plots[[j]]
  dev.off()
}


