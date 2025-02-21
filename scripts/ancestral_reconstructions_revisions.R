library(ape)
library(ggtree)
library(ggplot2)
library(viridis)
library(here)
library(patchwork)
library(mvSLOUCH)
library(geiger)
library(phytools)

setwd(here())

# Source this gist (From "https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32"), which allows using multiple of the same aes for the same plot
# Here I use it to make multiple heatmaps for Annika's figure (but should work better than the super hack I used for Figure 5 of the Cavefish single-cell paper)

# setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid_sleep_videos/_analysis2/R_plots/")
source("/Volumes/BZ/Home/gizevo30/R_Projects/Plot-multiple-scales-same-aes-ggplot2.R")

cichlids.diel <- read.csv("pheno_data/cichlid_pc-loadings_eco-morph_rest_full.csv")
lt_phylo <- read.nexus("scripts/05_BEAST_RAxML.tre")

cichlids.diel$tips <- cichlids.diel$six_letter_name_Ronco

cichlids.diel$x <- 11

subset <- keep.tip(lt_phylo, tip = cichlids.diel$tips[cichlids.diel$tips %in% lt_phylo$tip.label])
d <- fortify(subset)
cichlids.diel <- cichlids.diel[cichlids.diel$tips %in% subset$tip.label,]
cichlids.diel$tips <- factor(cichlids.diel$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
cichlids.diel$y <- as.numeric(cichlids.diel$tips)

################################################################################################################################################################################################
#### Test phylogenetic models w/ geiger  #######################################################################################################################################################
################################################################################################################################################################################################
## Possible models c("BM","OU","EB","rate_trend","lambda","kappa","delta","mean_trend","white")

pc1 <- cichlids.diel$pc1
names(pc1) <- cichlids.diel$tips

pc1_q <- pc1 + pc1^2

pc2 <- cichlids.diel$pc2
names(pc2) <- cichlids.diel$tips

tr <- cichlids.diel$total_rest
names(tr) <- cichlids.diel$tips


### Run using geiger
PC1_BM <- fitContinuous(subset, pc1, model = "BM")
PC1_OU <- fitContinuous(subset, pc1, model = "OU")
PC1_EB <- fitContinuous(subset, pc1, model = "EB")


aic.vals.pc1 <- setNames(c(PC1_BM$opt$aicc,PC1_OU$opt$aicc,PC1_EB$opt$aicc),
                         c("BM","OU","EB"))

PC2_BM <- fitContinuous(subset, pc2, model = "BM")
PC2_OU <- fitContinuous(subset, pc2, model = "OU")
PC2_EB <- fitContinuous(subset, pc2, model = "EB")


aic.vals.pc2 <- setNames(c(PC2_BM$opt$aicc,PC2_OU$opt$aicc,PC2_EB$opt$aicc),
                         c("BM","OU","EB"))


## Run OU ancestral reconstruction using phytools
PT_PC1_OU <- anc.ML(subset, pc1, model = "OU")
PT_PC2_OU <- anc.ML(subset, pc2, model = "OU")

################################################################################################################################################################################################
##### make plots comparing PC1 vs PC2 for extant and internal nodes ############################################################################################################################
################################################################################################################################################################################################

# Determine ages of each node
node_heights <- nodeHeights(subset)
node.age <- node_heights[match(node.data$node, subset$edge[,2]),2] # This extracts all but one that is missing from phylo_tree$edge[,2], and generates an NA
node.age[node.data$node[is.na(node.age)]] <- 0 # This is the root
node.data$node.age <- node.age


ancestral.plot.pc1.brownian <- ggtree(subset, layout = "circular", size = 1) %<+% node.data + aes(color = PT_PC1_OU) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_tippoint(aes(color = PT_PC1_OU), shape = 16, size = 1) + scale_color_distiller(palette = "RdBu", direction = 1)

ancestral.plot.pc2.brownian <- ggtree(subset, layout = "circular", size = 1) %<+% node.data + aes(color = PT_PC2_OU) + scale_color_distiller(palette = "BrBG") + geom_tiplab(color = "black", size = 2, offset = 0.5) + geom_tippoint(aes(color = PT_PC2_OU), shape = 16, size = 1) + scale_color_distiller(palette = "BrBG", direction = 1)

### Identify pairs of points (parent -> daughter nodes)

edges <- subset$edge
new_df_test <- lapply(seq_along(c(1:118)), function(x) {
  
  df <- node.data[node.data$node %in% c(edges[x, 1], edges[x,2]),]
  df$pair <- paste(edges[x, 1], edges[x, 2], sep = "_")
  df$parent_daughter <- c("daughter", "parent")
  df$node.age.2 <- min(df$node.age)
  return(df)
})

new_df_test <- Reduce(rbind, new_df_test)

### Annotate this with some more information for plotting
new_df_test$internal <- ifelse(row.names(new_df_test) %in% subset$tip.label, "tip", "internal_node")
new_df_test$label <- ifelse(new_df_test$internal == "tip", row.names(new_df_test), NA)
new_df_test$highlight <- factor(ifelse(new_df_test$node %in% c(11,12,73,31,32,94,54,55,115,44,45,106,29,30,92, 14,15,76,50,51,111,33,34,96), "highlight", "lowlight"), levels = c("lowlight", "highlight"))


#### Make plots
nest_plot <- ggplot(new_df_test, aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age.2, fill = internal)) + geom_line(aes(group = pair), alpha = 0.25) + geom_point(shape = 21, size = 1, colour = "transparent") + scale_color_viridis() + theme_classic() + xlab("PC1 (diurnal - nocturnal preference)") + ylab("PC2 (crepuscular preference)")

# labelled_nest_plot <- ggplot(new_df_test, aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age, fill = internal, label = label)) + geom_line(aes(group = pair), size = 1, alpha = 0.75) + geom_point(shape = 21, size = 3, colour = "transparent") + geom_text_repel(colour = "black") + scale_color_viridis() + theme_classic()

subset_plot <- ggplot(new_df_test[new_df_test$highlight == "highlight",], aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age.2, fill = internal, label = label)) + geom_line(aes(group = pair), size = 1, alpha = 0.75) + geom_point(shape = 21, size = 1, colour = "transparent") + geom_text_repel(colour = "black", size = 2) + scale_color_viridis_c(limits = c(0, 9.080887)) + theme_classic()
subset_plot <- subset_plot + xlim(c(min(labelled_nest_plot$data$PT_PC1_OU),max(labelled_nest_plot$data$PT_PC2_OU))) + ylim(c(min(labelled_nest_plot$data$PC2_brownian),max(labelled_nest_plot$data$PC2_brownian)))  + xlab("PC1 (diurnal - nocturnal preference)") + ylab("PC2 (crepuscular preference)")

################################################################################################################################################################################################
########## output the figures here #############################################################################################################################################################
################################################################################################################################################################################################

final_plot <- ancestral.plot.pc1.brownian + ancestral.plot.pc2.brownian + nest_plot + subset_plot + plot_layout(ncol = 4, guides = "collect", width = unit(c(35), "mm"), height = unit(c(35), "mm"))
final_plot

pdf("figure_ancestral_reconstructions.pdf", width = 20, height = 15)
final_plot 
dev.off()
