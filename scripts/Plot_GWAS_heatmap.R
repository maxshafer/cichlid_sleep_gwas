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
library(gsheet)
library(ggdendro)
library(tidyr)

######## OK, this script is to take the lists of annotated genes (top50 for ps, PGLS, and comb_ps) for total_rest and day_night_dif and make plots based on ontology analysis (ND and SC genes)

### LOAD IN GOOGLE SHEET DATA ### 

# Load in the data from google sheets
# and do some clean up on the data

url <- 'https://docs.google.com/spreadsheets/d/1MIzkbeFkARwfA2YEAfV1v1_2eECi_A-tZ5PEwYMZAGY/edit#gid=435005342'
gwas <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

# The above gives me pvalues and some ontology annotation (should improve this at some point)
# I need a phylogeny, and allele frequencies for these snps, and the phenotypes as well.
# I want to make a big heatmap to compare across species

## Load allele frequncies

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/sra_reads_nobackup/alleleFreqs")

file_list <- list.files()
file_list <- file_list[grep("_fromProbabilities_AF.txt", file_list)]

if (!(exists("allele_freq"))) {
  allele_freq <- lapply(file_list, function(x) fread(x))
}

names(allele_freq) <- substr(file_list, start = 52, stop = 60)


# Load the phylogeny

lt_phylo <- read.nexus("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/05_BEAST_RAxML.tre")

# Load phenotype file for plotting and annotating with frequencies
phenotypes <- read.csv("../../pheno_data/combined_cichlid_data_2022-08-04.csv", row.names = "X")


# OK, need to pull the frequencies for each of the SNPs (in gwas)
phenotypes_dnd <- phenotypes
for (i in 1:nrow(gwas[gwas$comparison == "day_night_dif",])) {
  gene <- gwas$gene[i]
  location <- gwas$location[i]
  chromosome <- str_split(location, ":")[[1]][[1]]
  coord <- as.numeric(str_split(location, ":")[[1]][[2]])
  
  af <- as.data.frame(allele_freq[[chromosome]])
  af2 <- af[af$coord %in% coord,5:65]
  
  phenotypes_dnd <- phenotypes_dnd[phenotypes_dnd$six_letter_name_Ronco %in% names(af2),]
  
  phenotypes_dnd[,paste("genotype", gene, location, sep = "_")] <- as.numeric(af2)
  
}

phenotypes_tr <- phenotypes
for (i in 1:nrow(gwas[gwas$comparison == "total_rest",])) {
  gene <- gwas$gene[i]
  location <- gwas$location[i]
  chromosome <- str_split(location, ":")[[1]][[1]]
  coord <- as.numeric(str_split(location, ":")[[1]][[2]])
  
  af <- as.data.frame(allele_freq[[chromosome]])
  af2 <- af[af$coord %in% coord,5:65]
  
  phenotypes_tr <- phenotypes_tr[phenotypes_tr$six_letter_name_Ronco %in% names(af2),]
  
  phenotypes_tr[,paste("genotype", gene, location, sep = "_")] <- as.numeric(af2)
  
}

## There are some -1 values (don't know what they mean), but can replace with 0s or NAs?

phenotypes_tr[phenotypes_tr == -1] <- NA
phenotypes_dnd[phenotypes_dnd == -1] <- NA


###### ###### ###### ###### ###### ###### ###### ###### ###### 
##### Make dendrogram for species
###### ###### ###### ###### ###### ###### ###### ###### ###### 

phenotypes_dnd_scaled <- phenotypes_dnd

# Scale the values, or not
# phenotypes_dnd_scaled[,c(8:ncol(phenotypes_dnd_scaled))] <- scale(phenotypes_dnd_scaled[,c(8:ncol(phenotypes_dnd_scaled))])

# Run clustering
phenotypes_dnd_matrix <- as.matrix(phenotypes_dnd_scaled[, -c(1:7)])
rownames(phenotypes_dnd_matrix) <- phenotypes_dnd_scaled$six_letter_name_Ronco
phenotypes_dnd_dendro <- as.phylo(hclust(d = dist(x = phenotypes_dnd_matrix)))

# Create dendro
dendro_dnd_plot <- ggtree(phenotypes_dnd_dendro) # + geom_tiplab() 

phenotypes_dnd_scaled$y <- dendro_dnd_plot$data$y[match(phenotypes_dnd_scaled$six_letter_name_Ronco, dendro_dnd_plot$data$label)]

dendro_dnd_plot <- dendro_dnd_plot + geom_tile(data = phenotypes_dnd_scaled, aes(y=y, x=7, fill = day_night_dif), width = 1, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu")
# Preview the plot
# dendro_dnd_plot


###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Make dendrogram for SNPS
###### ###### ###### ###### ###### ###### ###### ###### ###### 

phenotypes_dnd_scaled_r <- t(phenotypes_dnd)
colnames(phenotypes_dnd_scaled_r) <- phenotypes_dnd$six_letter_name_Ronco

# Scale the values, or not
# phenotypes_dnd_scaled_r[c(8:nrow(phenotypes_dnd_scaled_r)),] <- scale(as.numeric(phenotypes_dnd_scaled_r[c(8:nrow(phenotypes_dnd_scaled_r)),]))

# Run clustering
phenotypes_dnd_matrix_r <- as.matrix(phenotypes_dnd_scaled_r[-c(1:7),])
phenotypes_dnd_dendro_r <- as.phylo(hclust(d = dist(x = phenotypes_dnd_matrix_r)))

# Create dendro
dendro_dnd_plot_r <- ggtree(phenotypes_dnd_dendro_r) + coord_flip()

# Preview the plot
# print(dendro_dnd_plot_r)

# Heatmap
# Data wrangling
phenotypes_dnd_long <- pivot_longer(data = phenotypes_dnd_scaled,
                           cols = -c(six_letter_name_Ronco, total_rest, peak_amplitude, peak, day_night_dif, cluster, day_night_dif_spd, y),
                           names_to = "measurement",
                           values_to = "value")

phenotypes_dnd_order <- dendro_dnd_plot$data$label[match(c(1:length(phenotypes_dnd_dendro$tip.label)), dendro_dnd_plot$data$y)]
phenotypes_dnd_order_r <- dendro_dnd_plot_r$data$label[match(c(1:length(phenotypes_dnd_dendro_r$tip.label)), dendro_dnd_plot_r$data$y)]

# Order the levels according to their position in the cluster
phenotypes_dnd_long$six_letter_name_Ronco <- factor(x = phenotypes_dnd_long$six_letter_name_Ronco,
                               levels = phenotypes_dnd_order, 
                               ordered = TRUE)

phenotypes_dnd_long$measurement <- factor(x = phenotypes_dnd_long$measurement,
                                                    levels = phenotypes_dnd_order_r, 
                                                    ordered = TRUE)


heatmap_plot_dnd <- ggplot(data = phenotypes_dnd_long, aes(x = measurement, y = six_letter_name_Ronco)) + geom_tile(aes(fill = value)) + theme(axis.text = element_text(size = 6), axis.text.x = element_text(angle = 90), axis.title = element_blank()) + scale_fill_distiller(palette = "Greys", direction = 1) #+ scale_fill_gradient2(midpoint = 0.5)

# Preview the heatmap
# print(heatmap_plot_dnd)


# dendro_dnd_plot_r + plot_spacer() + heatmap_plot_dnd + dendro_dnd_plot



###############################################################################################################################################################################################################
###############################################################################################################################################################################################################




###### ###### ###### ###### ###### ###### ###### ###### ###### 
##### Make dendrogram for species
###### ###### ###### ###### ###### ###### ###### ###### ###### 

phenotypes_tr_scaled <- phenotypes_tr

# Scale the values, or not
# phenotypes_tr_scaled[,c(8:ncol(phenotypes_tr_scaled))] <- scale(phenotypes_tr_scaled[,c(8:ncol(phenotypes_tr_scaled))])

# Run clustering
phenotypes_tr_matrix <- as.matrix(phenotypes_tr_scaled[, -c(1:7)])
rownames(phenotypes_tr_matrix) <- phenotypes_tr_scaled$six_letter_name_Ronco
phenotypes_tr_dendro <- as.phylo(hclust(d = dist(x = phenotypes_tr_matrix)))


# Create dendro
dendro_tr_plot <- ggtree(phenotypes_tr_dendro) #+ geom_tiplab() 
phenotypes_tr_scaled$y <- dendro_tr_plot$data$y[match(phenotypes_tr_scaled$six_letter_name_Ronco, dendro_tr_plot$data$label)]

dendro_tr_plot <- dendro_tr_plot + geom_tile(data = phenotypes_tr_scaled, aes(y=y, x=5, fill = total_rest), width = 1, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "PRGn", direction = 1) # + scale_fill_distiller(palette = "RdBu")
# Preview the plot
# dendro_tr_plot

###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Make dendrogram for SNPS
###### ###### ###### ###### ###### ###### ###### ###### ###### 

phenotypes_tr_scaled_r <- t(phenotypes_tr)
colnames(phenotypes_tr_scaled_r) <- phenotypes_tr$six_letter_name_Ronco

# Scale the values, or not
# phenotypes_tr_scaled_r[c(8:nrow(phenotypes_tr_scaled_r)),] <- scale(as.numeric(phenotypes_tr_scaled_r[c(8:nrow(phenotypes_tr_scaled_r)),]))

# Run clustering
phenotypes_tr_matrix_r <- as.matrix(phenotypes_tr_scaled_r[-c(1:7),])
phenotypes_tr_dendro_r <- as.phylo(hclust(d = dist(x = phenotypes_tr_matrix_r)))

# Create dendro
dendro_tr_plot_r <- ggtree(phenotypes_tr_dendro_r) + coord_flip()

# Preview the plot
# print(dendro_tr_plot_r)



# Heatmap
# Data wrangling
phenotypes_tr_long <- pivot_longer(data = phenotypes_tr_scaled, cols = -c(six_letter_name_Ronco, total_rest, peak_amplitude, peak, day_night_dif, cluster, day_night_dif_spd, y), names_to = "measurement", values_to = "value")


phenotypes_tr_order <- dendro_tr_plot$data$label[match(c(1:length(phenotypes_tr_dendro$tip.label)), dendro_tr_plot$data$y)]
phenotypes_tr_order_r <- dendro_tr_plot_r$data$label[match(c(1:length(phenotypes_tr_dendro_r$tip.label)), dendro_tr_plot_r$data$y)]

# Order the levels according to their position in the cluster
phenotypes_tr_long$six_letter_name_Ronco <- factor(x = phenotypes_tr_long$six_letter_name_Ronco, levels = phenotypes_tr_order, ordered = TRUE)

phenotypes_tr_long$measurement <- factor(x = phenotypes_tr_long$measurement, levels = phenotypes_tr_order_r, ordered = TRUE)


heatmap_plot_tr <- ggplot(data = phenotypes_tr_long, aes(x = measurement, y = six_letter_name_Ronco)) + geom_tile(aes(fill = value)) + theme(axis.text = element_text(size = 6), axis.text.x = element_text(angle = 90), axis.title = element_blank()) + scale_fill_distiller(palette = "Greys", direction = 1) #+ scale_fill_gradient2(midpoint = 0.5)

# Preview the heatmap
# print(heatmap_plot_tr)




### OK above works well, can I add the phenotype information?

dendro_dnd_plot + heatmap_plot_dnd + plot_spacer() + dendro_dnd_plot_r + plot_layout(guides = "collect", widths = c(5,10), heights = c(10,5))
dendro_tr_plot + heatmap_plot_tr + plot_spacer() + dendro_tr_plot_r + plot_layout(guides = "collect", widths = c(5,10), heights = c(10,5))










