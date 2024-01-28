library(caper)
library(pls)
library(geomorph)
library(ggplot2)
library(gsheet)
library(ggrepel)
library(dplyr)
library(patchwork)
library(viridis)
library(ape)
library(ggtree)

## Still on cluster, no worky
# source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")


########################################################################################################
################ Load and prep input data    ###########################################################
########################################################################################################

## Meta data has tribe, diet, habitat info
url <- 'https://docs.google.com/spreadsheets/d/1m6cV5W-UjDdIA6_ladE0H4iCwQU_xcbWwnQ5w2-kgbM/edit?usp=sharing'
meta_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

## Import raw normalized speed data and transpose it (so that row names are species)
speed_data <- read.csv("pca_input_zscore.csv", row.names = "time_of_day")
speed_data <- t(speed_data)

## Load consistancy metrics
intra_cons <- read.csv("pheno_data/intra-individual_corr_species_mean.csv")
inter_cons <- read.csv("pheno_data/inter-individual_corr_weekly_species_mean.csv")

## load in old and new data, so I can have total_rest, as well as Annika's pc values
old_data <- read.csv("pheno_data/combined_cichlid_data_2022-11-17_new_peaks.csv")
sleepData <- read.csv("cichild_corr_data.csv")
sleepData$total_rest <- old_data$total_rest[match(sleepData$species, old_data$six_letter_name_Ronco)]
sleepData$tribe <- meta_data$tribe[match(sleepData$species, meta_data$six_letter_name_Ronco)]

sleepData$inter_coef <- inter_cons$corr_coef[match(sleepData$species, inter_cons$species)]
sleepData$intra_coef <- intra_cons$corr_coef[match(sleepData$species, intra_cons$species)]

rownames(sleepData) <- sleepData$species

## Add or keep shell dwellers? Or other subsetting before tree subsetting
sleepData <- sleepData[!(sleepData$species %in% c("Neomul", "Neobre")),]


## Load and trim tree

phylo_lt <- read.nexus("scripts/05_BEAST_RAxML.tre")
cichlidTree <- keep.tip(phylo_lt, tip = sleepData$species[sleepData$species %in% phylo_lt$tip.label])

## Double trim with input data
sleepData <- sleepData[sleepData$species %in% cichlidTree$tip.label,]
speed_data <- speed_data[row.names(speed_data) %in% cichlidTree$tip.label,]


sleepData$tips <- sleepData$species
sleepData$x <- 11
d <- fortify(cichlidTree)
sleepData$tips <- factor(sleepData$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
sleepData$y <- as.numeric(sleepData$tips)

########################################################################################################
################ Reconstruct Ancestral States   ########################################################
########################################################################################################

## Probably want to do PC1, PC2, total_rest, and inter/intra coefs

pc1 <- sleepData$pc1
names(pc1) <- sleepData$species
pc2 <- sleepData$pc2
names(pc2) <- sleepData$species
total_rest <- sleepData$total_rest
names(total_rest) <- sleepData$tips
inter_coef <- sleepData$inter_coef
names(inter_coef) <- sleepData$species
intra_coef <- sleepData$intra_coef
names(intra_coef) <- sleepData$species

fitBM_pc1 <- ace(pc1, cichlidTree, method = "ML", model = "BM", type = "continuous")
fitBM_pc2 <- ace(pc2, cichlidTree, method = "ML", model = "BM", type = "continuous")
fitBM_total_rest <- ace(total_rest, cichlidTree, method = "ML", model = "BM", type = "continuous")
fitBM_inter <- ace(inter_coef, cichlidTree, method = "ML", model = "BM", type = "continuous")
fitBM_intra <- ace(intra_coef, cichlidTree, method = "ML", model = "BM", type = "continuous")

lik.anc <- data.frame(PC1 = fitBM_pc1$ace, PC2 = fitBM_pc2$ace, Total_Rest = fitBM_total_rest$ace, Inter_Coef = fitBM_inter$ace, Intra_Coef = fitBM_intra$ace)
lik.anc$node <- (1:nrow(lik.anc)) + length(cichlidTree$tip.label)

node.data <- data.frame(PC1 = pc1, PC2 = pc2, Total_Rest = total_rest, Inter_Coef = inter_coef, Intra_Coef = intra_coef, node = match(names(pc1), cichlidTree$tip.label))
node.data <- rbind(node.data, lik.anc)

################################################################################################################################################################################################
#### Identify nodes for tribes  ################################################################################################################################################################
################################################################################################################################################################################################

ggtree(as.phylo(cichlidTree), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab() + geom_text(aes(label = node))

lamps <- 57
ectos <- 92
trops <- 102
cyps <- 99
erets <- 47
# limnos <- 42
cyphos <- 36
bouls <- 55
haplos <- 53

################################################################################################################################################################################################
#### Make Plots  ###############################################################################################################################################################################
################################################################################################################################################################################################

phylo.plot <- ggtree(as.phylo(cichlidTree), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 2.5)

pc1_plot <- phylo.plot + geom_tile(data = sleepData, aes(y=y, x=x, fill = pc1), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu") #+ scale_color_discrete(na.value = 'transparent')
pc2_plot <- phylo.plot + geom_tile(data = sleepData, aes(y=y, x=x, fill = pc2), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "BrBG") #+ scale_color_discrete(na.value = 'transparent')
tr_plot <- phylo.plot + geom_tile(data = sleepData, aes(y=y, x=x, fill = total_rest), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "PRGn") #+ scale_color_discrete(na.value = 'transparent')
inter_plot <- phylo.plot + geom_tile(data = sleepData, aes(y=y, x=x, fill = inter_coef), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "YlOrRd", limits = c(0,1)) #+ scale_color_discrete(na.value = 'transparent')
intra_plot <- phylo.plot + geom_tile(data = sleepData, aes(y=y, x=x, fill = intra_coef), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "YlOrRd", limits = c(0,1)) #+ scale_color_discrete(na.value = 'transparent')

# # These don't work without the custom script (That is on basel cluster)
# p3 <- p3 + new_scale("size") + new_scale("colour") + geom_point(data = sleepData, aes(y=y, x=x+1.75, size = total_rest, colour = total_rest)) + scale_colour_distiller(palette = "PRGn", direction = 1) 
# p3 <- p3 + new_scale("fill") + new_scale("colour") + geom_tile(data = sleepData, aes(y=y, x=x+3, fill = pc2), width = 1, inherit.aes = FALSE) + scale_fill_distiller(palette = "BrBG", direction = 1)
# # p3 <- p3 + new_scale("fill") + geom_tile(data = sleepData, aes(y=y, x=x+3, fill = intra_coef), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "YlOrRd")
# # p3 <- p3 + xlim(c(0,17))
# 
# pdf("Summary_figure_phylogram_circular.pdf", width = 10, height = 10)
# p3
# dev.off()
# 
# png("Summary_figure_phylogram_circular.png", width = 10, height = 10, units = "in", res = 500)
# p3
# dev.off()



##### Plot ancestral reconstructions #######

ancestral.plot.pc1 <- ggtree(cichlidTree, layout = "circular", size = 2) %<+% node.data + aes(color = PC1) + scale_color_distiller(palette = "RdBu", direction = 1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC1), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu", direction = 1)
ancestral.plot.total_rest <- ggtree(cichlidTree, layout = "circular", size = 2) %<+% node.data + aes(color = Total_Rest) + scale_color_distiller(palette = "PRGn", direction = 1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Total_Rest), shape = 16, size = 6) + scale_color_distiller(palette = "PRGn", direction = 1)
ancestral.plot.pc2 <- ggtree(cichlidTree, layout = "circular", size = 2) %<+% node.data + aes(color = PC2) + scale_color_distiller(palette = "BrBG", direction = -1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC2), shape = 16, size = 6) + scale_color_distiller(palette = "BrBG", direction = -1)
ancestral.plot.inter <- ggtree(cichlidTree, layout = "circular", size = 2) %<+% node.data + aes(color = Inter_Coef) + scale_color_viridis(direction = 1, limits = c(0,1)) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Inter_Coef), shape = 16, size = 6) + scale_color_viridis(direction = 1, limits = c(0,1))
ancestral.plot.intra <- ggtree(cichlidTree, layout = "circular", size = 2) %<+% node.data + aes(color = Intra_Coef) + scale_color_viridis(direction = 1, limits = c(0,1)) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Intra_Coef), shape = 16, size = 6) + scale_color_viridis(direction = 1, limits = c(0,1))

ancestral.plot <- ancestral.plot.pc1 + ancestral.plot.pc2 + ancestral.plot.total_rest + ancestral.plot.intra + plot_layout(nrow = 2)

################################################################################################################################################################################################
#### Save Plots  ###############################################################################################################################################################################
################################################################################################################################################################################################

pdf("Ancestral_phylograms_NEW.pdf", width = 20, height = 15)
ancestral.plot
dev.off()

png("Ancestral_phylograms_NEW.png", width = 20, height = 15, units = "in", res = 500)
ancestral.plot
dev.off()





########################################################################################################
################ Reconstruct Ancestral States for ALL PC   #############################################
########################################################################################################

all_pcs <- lapply(c(2:11), function(x) {
  fitBM <- ace(sleepData[,x], cichlidTree, method = "ML", model = "BM", type = "continuous")
  out <- fitBM$ace
  return(out)
})

all_pcs <- Reduce(cbind, all_pcs)
colnames(all_pcs) <- colnames(sleepData)[2:11]

pdf("cichlid_phylogram_nodes.pdf", width = 10, height = 10)
ggtree(as.phylo(cichlidTree), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 0.5) + geom_text(aes(label = node))
dev.off()

write.csv(all_pcs, file = "reconstructed_pcs.csv")
