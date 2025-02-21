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


#colnames(cichlids.diel) <- c("y", "tips", "species", "peak_amplitude", "peak", "day_night_dif")

cichlids.diel$tips <- cichlids.diel$six_letter_name_Ronco

cichlids.diel$x <- 11

subset <- keep.tip(lt_phylo, tip = cichlids.diel$tips[cichlids.diel$tips %in% lt_phylo$tip.label])
d <- fortify(subset)
cichlids.diel <- cichlids.diel[cichlids.diel$tips %in% subset$tip.label,]
cichlids.diel$tips <- factor(cichlids.diel$tips, levels = rev(with(d, label[order(y, decreasing=T)])))
cichlids.diel$y <- as.numeric(cichlids.diel$tips)

################################################################################################################################################################################################
#### Reconstruct ancestral states  #############################################################################################################################################################
################################################################################################################################################################################################

day_night_dif <- cichlids.diel$pc1
names(day_night_dif) <- cichlids.diel$tips
total_rest <- cichlids.diel$total_rest
names(total_rest) <- cichlids.diel$tips
peak <- cichlids.diel$pc2
names(peak) <- cichlids.diel$tips

fitBM_day_night_dif <- ace(day_night_dif, subset, method = "ML", model = "BM", type = "continuous")
fitBM_total_rest <- ace(total_rest, subset, method = "ML", model = "BM", type = "continuous")
fitBM_peak <- ace(peak, subset, method = "ML", model = "BM", type = "continuous")

lik.anc <- data.frame(Day_vs_Night_activity = fitBM_day_night_dif$ace, Total_Rest = fitBM_total_rest$ace, Peak_Percentage = fitBM_peak$ace)
lik.anc$node <- (1:nrow(lik.anc)) + length(subset$tip.label)

node.data <- data.frame(Day_vs_Night_activity = day_night_dif, Total_Rest = total_rest, Peak_Percentage = peak, node = match(names(day_night_dif), subset$tip.label))
node.data <- rbind(node.data, lik.anc)

################################################################################################################################################################################################
#### Reconstruct ancestral states w/ geiger  ###################################################################################################################################################
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
PC1_rate_trend <- fitContinuous(subset, pc1, model = "rate_trend")
PC1_lambda <- fitContinuous(subset, pc1, model = "lambda")
PC1_kappa <- fitContinuous(subset, pc1, model = "kappa")
PC1_delta <- fitContinuous(subset, pc1, model = "delta")
PC1_mean_trend <- fitContinuous(subset, pc1, model = "mean_trend")
PC1_white <- fitContinuous(subset, pc1, model = "white")

aic.vals.pc1 <- setNames(c(PC1_BM$opt$aicc,PC1_OU$opt$aicc,PC1_EB$opt$aicc,PC1_rate_trend$opt$aicc,PC1_lambda$opt$aicc, PC1_kappa$opt$aicc,PC1_delta$opt$aicc,PC1_mean_trend$opt$aicc,PC1_white$opt$aicc),
                   c("BM","OU","EB","rate_trend","lambda","kappa","delta","mean_trend","white"))

PC2_BM <- fitContinuous(subset, pc2, model = "BM")
PC2_OU <- fitContinuous(subset, pc2, model = "OU")
PC2_EB <- fitContinuous(subset, pc2, model = "EB")
PC2_rate_trend <- fitContinuous(subset, pc2, model = "rate_trend")
PC2_lambda <- fitContinuous(subset, pc2, model = "lambda")
PC2_kappa <- fitContinuous(subset, pc2, model = "kappa")
PC2_delta <- fitContinuous(subset, pc2, model = "delta")
PC2_mean_trend <- fitContinuous(subset, pc2, model = "mean_trend")
PC2_white <- fitContinuous(subset, pc2, model = "white")

aic.vals.pc2 <- setNames(c(PC2_BM$opt$aicc,PC2_OU$opt$aicc,PC2_EB$opt$aicc,PC2_rate_trend$opt$aicc,PC2_lambda$opt$aicc, PC2_kappa$opt$aicc,PC2_delta$opt$aicc,PC2_mean_trend$opt$aicc,PC2_white$opt$aicc),
                         c("BM","OU","EB","rate_trend","lambda","kappa","delta","mean_trend","white"))


## Run using phytools
# PT_PC1_BM <- anc.ML(subset, pc1, model = "BM")
PT_PC1_OU <- anc.ML(subset, pc1, model = "OU")
# PT_PC1_EB <- anc.ML(subset, pc1, model = "EB")

# PT_PC2_BM <- anc.ML(subset, pc2, model = "BM")
PT_PC2_OU <- anc.ML(subset, pc2, model = "OU")
# PT_PC2_EB <- anc.ML(subset, pc2, model = "EB")

## Can also do a variable rates BM model?

# test <- anc.ML(subset, pc1, model = "EB")
# test2 <- anc.ML(subset, pc1, model = "BM")

# mrBM <- multirateBM(subset, pc2)
# mrBM2 <- multirateBM(subset, pc2, lambda = 0.1)
# mrBM3 <- multirateBM(subset, pc2, lambda = 0.01)
# mrBM4 <- multirateBM(subset, pc1, lambda = 0.5)
# mrBM5 <- multirateBM(subset, pc2, lambda = 10)

# pc1_mr <- list(mrBM, mrBM2, mrBM3, mrBM4, mrBM5)
# pc2_mr <- list(mrBM, mrBM2, mrBM3)

# saveRDS(list(pc1_mr, pc2_mr), "phytools_multirate_objects.rds")

mrBM_models <- readRDS("phytools_multirate_objects.rds")
pc1_mr <- mrBM_models[[1]]
pc2_mr <- mrBM_models[[2]]

### make a function that takes the variable rates, modifies the branch lengths, then runs the BM model

fitModelfromMR <- function(phylo = phylo, trait_vector = trait_vector, mr_model = mr_model) {
  
  rates <- mr_model$sig2
  names(rates) <- 1:length(rates)
  
  ## I need to multiply the rate and the branch.lengths, based on the matching second column in edges
  rates_2 <- rates[match(phylo$edge[,2], names(rates))]
  
  phylo2 <- phylo
  branch.lengths <- phylo$edge.length*rates_2
  branch.lengths <- (branch.lengths/max(branch.lengths) ) * max(phylo$edge.length)
  phylo2$edge.length <- branch.lengths
  
  fitBM <- ace(trait_vector, phylo2, method = "ML", model = "BM", type = "continuous")
  
  return(fitBM)
  
}


####################################################################################################


## Compare with PC2 reconstruction
## The name of rates_2 is the tip that that branch point ends in
pc_2 <- c(pc2[subset$tip.label], fitBM_pc2$ace)


new_df <- data.frame(branch_rates = rates, nodes = names(rates), pc2 = pc_2)

## Can I now link parental node PC2 values to descending edge branch rates? Some PC2 values will be attached to multiple branches
## I can use the tree$edge to map them (parental -> daughter connections)
## will only need internal node PC2 values (fitBM_pc2$ace)

new_df2 <- data.frame(parent = subset$edge[,1], daughter = subset$edge[,2])
new_df2$pc2 <- new_df$pc2[match(new_df2$parent, new_df$nodes)]
new_df2$branch_rates <- new_df$branch_rates[match(new_df2$daughter, new_df$nodes)]

####################################################################################################

# subset2 <- subset
# branch.lengths <- subset$edge.length*rates_2
# branch.lengths <- (branch.lengths/max(branch.lengths) ) * max(subset$edge.length)
# subset2$edge.length <- branch.lengths
# 
# fitBM_pc1_2 <- ace(peak, subset2, method = "ML", model = "BM", type = "continuous")
# fitBM_pc2_2 <- ace(peak, subset2, method = "ML", model = "BM", type = "continuous")



fitBM_pc1_lmb1 <- fitModelfromMR(phylo = subset, trait_vector = day_night_dif, mr_model = pc1_mr[[1]])
fitBM_pc1_lmb0.1 <- fitModelfromMR(phylo = subset, trait_vector = day_night_dif, mr_model = pc1_mr[[2]])
fitBM_pc1_lmb0.01 <- fitModelfromMR(phylo = subset, trait_vector = day_night_dif, mr_model = pc1_mr[[3]])
fitBM_pc1_lmb0.5 <- fitModelfromMR(phylo = subset, trait_vector = day_night_dif, mr_model = pc1_mr[[4]])

fitBM_pc2_lmb1 <- fitModelfromMR(phylo = subset, trait_vector = peak, mr_model = pc2_mr[[1]])
fitBM_pc2_lmb0.1 <- fitModelfromMR(phylo = subset, trait_vector = peak, mr_model = pc2_mr[[2]])
fitBM_pc2_lmb0.01 <- fitModelfromMR(phylo = subset, trait_vector = peak, mr_model = pc2_mr[[3]])
fitBM_pc2_lmb0.5 <- fitModelfromMR(phylo = subset, trait_vector = peak, mr_model = pc2_mr[[4]])




## Test plotting

lik.anc <- data.frame(PC1_brownian = fitBM_day_night_dif$ace, Total_Rest = fitBM_total_rest$ace, PC2_brownian = fitBM_peak$ace, PC1_mrBM = fitBM_pc1_lmb1$ace, PC1_mrBM2 = fitBM_pc1_lmb0.1$ace, PC1_mrBM3 = fitBM_pc1_lmb0.01$ace, PC1_mrBM4 = fitBM_pc1_lmb0.5$ace, PT_PC1_OU = PT_PC1_OU$ace, PC2_mrBM = fitBM_pc2_lmb1$ace, PC2_mrBM2 = fitBM_pc2_lmb0.1$ace, PC2_mrBM3 = fitBM_pc2_lmb0.01$ace, PC2_mrBM4 = fitBM_pc2_lmb0.5$ace, PT_PC2_OU = PT_PC2_OU$ace)
lik.anc$node <- (1:nrow(lik.anc)) + length(subset$tip.label)

node.data <- data.frame(PC1_brownian = pc1, Total_Rest = total_rest, PC2_brownian = pc2, PC1_mrBM = pc1, PC1_mrBM2 = pc1, PC1_mrBM3 = pc1, PC1_mrBM4 = pc1, PT_PC1_OU = pc1,PC2_mrBM = pc2, PC2_mrBM2 = pc2, PC2_mrBM3 = pc2, PC2_mrBM4 = pc2, PT_PC2_OU = pc2, node = match(names(day_night_dif), subset$tip.label))
node.data <- rbind(node.data, lik.anc)

# ancestral.plot.pc1 <- ggtree(subset, layout = "rectangular", size = 2) %<+% node.data + aes(color = PC1_mrBM) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC1_mrBM), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu", direction = 1)
# ancestral.plot.pc1_2 <- ggtree(subset, layout = "rectangular", size = 2) %<+% node.data + aes(color = PC1_mrBM2) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC1_mrBM2), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu", direction = 1)
# ancestral.plot.pc1_3 <- ggtree(subset, layout = "rectangular", size = 2) %<+% node.data + aes(color = PC1_mrBM3) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC1_mrBM3), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu", direction = 1)
# ancestral.plot.pc1.brownian <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Day_vs_Night_activity) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Day_vs_Night_activity), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu", direction = 1)
# 
# 
# ancestral.plot.pc1.brownian + ancestral.plot.pc1 + ancestral.plot.pc1_3 + ancestral.plot.pc1_2 + plot_layout(ncol = 4)
# 
# 
# 
# ancestral.plot.pc2 <- ggtree(subset, layout = "rectangular", size = 2) %<+% node.data + aes(color = PC2_mrBM) + scale_color_distiller(palette = "BrBG") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC2_mrBM), shape = 16, size = 6) + scale_color_distiller(palette = "BrBG", direction = 1)
# ancestral.plot.pc2_2 <- ggtree(subset, layout = "rectangular", size = 2) %<+% node.data + aes(color = PC2_mrBM2) + scale_color_distiller(palette = "BrBG") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC2_mrBM2), shape = 16, size = 6) + scale_color_distiller(palette = "BrBG", direction = 1)
# ancestral.plot.pc2_3 <- ggtree(subset, layout = "rectangular", size = 2) %<+% node.data + aes(color = PC2_mrBM3) + scale_color_distiller(palette = "BrBG") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = PC2_mrBM3), shape = 16, size = 6) + scale_color_distiller(palette = "BrBG", direction = 1)
# ancestral.plot.pc2.brownian <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Peak_Percentage) + scale_color_distiller(palette = "BrBG") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Peak_Percentage), shape = 16, size = 6) + scale_color_distiller(palette = "BrBG", direction = 1)
# 
# 
# ancestral.plot.pc2.brownian + ancestral.plot.pc2 + ancestral.plot.pc2_3 + ancestral.plot.pc2_2 + plot_layout(ncol = 4)
# 
# 
# 
# 
# ancestral.plot.pc1.brownian + ancestral.plot.pc1 + ancestral.plot.pc1_2 + ancestral.plot.pc2.brownian + ancestral.plot.pc2 + ancestral.plot.pc2_2 + plot_layout(ncol = 3)




# Determine ages of each node
node_heights <- nodeHeights(subset)

# ancestor <- round(ancestral_states$recon_states[Ntip(subset)+1])

node.age <- node_heights[match(node.data$node, subset$edge[,2]),2] # This extracts all but one that is missing from phylo_tree$edge[,2], and generates an NA
node.age[node.data$node[is.na(node.age)]] <- 0 # This is the root
# ancestral_states$node.age[is.na(ancestral_states$node.age)] <- 0
#node.age <- (node.age - max(node.age))*-1

node.data$node.age <- node.age




##### make plots comparing PC1 vs PC2 for extant and internal nodes

# ggplot(node.data, aes(x = Day_vs_Night_activity, y = Peak_Percentage, colour = node.age)) + geom_point()
# ggplot(node.data, aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age)) + geom_point() # Lowest AIC scores (should also figure out AICc)

### Make paired plot
### Can colour by parent/daughter, but just need a column with individual pair labels (1, 2, 3, 4, etc), then do a geom_line(aes(group = grouping))
### BUT, they are triplets - one parent and two daughters - I could put each triplet in a group, but the lines might not connect in the correct order (daughter-parent-daughter)
### or alternatively, duplicate parent nodes in the data (many versions of ancestral nodes)
### 118 edges means 118 pairs, or 236 entries in node.data (where there are 119 now)

## can I seq along subset$edge, and add rows to a new df by extracting them from node.data?
edges <- subset$edge
new_df_test <- lapply(seq_along(c(1:118)), function(x) {
  
  df <- node.data[node.data$node %in% c(edges[x, 1], edges[x,2]),]
  df$pair <- paste(edges[x, 1], edges[x, 2], sep = "_")
  df$parent_daughter <- c("daughter", "parent")
  return(df)
})

new_df_test <- Reduce(rbind, new_df_test)
new_df_test$internal <- ifelse(row.names(new_df_test) %in% subset$tip.label, "tip", "internal_node")
new_df_test$label <- ifelse(new_df_test$internal == "tip", row.names(new_df_test), NA)
new_df_test$highlight <- factor(ifelse(new_df_test$node %in% c(11,12,73,31,32,94,54,55,115,44,45,106,29,30,92, 14,15,76,50,51,111,33,34,96), "highlight", "lowlight"), levels = c("lowlight", "highlight"))

nest_plot <- ggplot(new_df_test, aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age, fill = internal)) + geom_line(aes(group = pair)) + geom_point(shape = 21, size = 3, colour = "transparent") + scale_color_viridis() + theme_classic()

labelled_nest_plot <- ggplot(new_df_test, aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age, fill = internal, label = label)) + geom_line(aes(group = pair), size = 1, alpha = 0.75) + geom_point(shape = 21, size = 3, colour = "transparent") + geom_text_repel(colour = "black") + scale_color_viridis() + theme_classic()
# ggplot(new_df_test, aes(x = PC1_mrBM, y = PC2_mrBM, colour = node.age)) + geom_point() + geom_line(aes(group = pair)) + scale_color_viridis() + theme_classic()
# ggplot(new_df_test, aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age)) + geom_point() + geom_line(aes(group = pair)) + scale_color_viridis() + theme_classic()

subset_plot <- ggplot(new_df_test[new_df_test$highlight == "highlight",], aes(x = PT_PC1_OU, y = PT_PC2_OU, colour = node.age, fill = internal, label = label)) + geom_line(aes(group = pair), size = 1, alpha = 0.75) + geom_point(shape = 21, size = 3, colour = "transparent") + geom_text_repel(colour = "black") + scale_color_viridis() + theme_classic()
subset_plot <- subset_plot + xlim(c(min(labelled_nest_plot$data$PT_PC1_OU),max(labelled_nest_plot$data$PT_PC2_OU))) + ylim(c(min(labelled_nest_plot$data$PC2_brownian),max(labelled_nest_plot$data$PC2_brownian)))

########## OK I can output the figures here ##############

ancestral.plot.pc1.brownian + ancestral.plot.pc2.brownian + nest_plot + subset_plot + plot_layout(ncol = 4, guides = "collect")






#### Can I plot pairwise delta's in PC1, relative to the parent node's PC2?
#### in new_df_test, take the 4 entries with the same first node pair # is the same, delta the diff between PC1 for those where parent_daughter == "daughter", take PC2 value where parent_daughter == "parent"
#### The below works, but I'm not sure that it's helpful - all of hte ancestral nodes are likely to be considered 'crepuscular', so I'm not sure how to intrepret it
#### I think just highlighting examples would be helpful enough

new_df_test$first_node <- sapply(strsplit(new_df_test$pair,"_"), `[`, 1)

output <- lapply(unique(new_df_test$first_node), function(x) {
  
  df <- new_df_test[new_df_test$first_node == x,]
  
  if(any(row.names(df) %in% subset$tip.label)) {
    pc1_delta <- abs(diff(df[df$parent_daughter == "daughter", "PC1_brownian"]))
  } else {
    pc1_delta <- abs(diff(df[df$parent_daughter == "parent", "PC2_brownian"]))
  }
  
  pc2 <- df$PT_PC2_OU[df$node == x][1]
  node.age <- min(df$node.age)
  
  # if(any(row.names(df) %in% subset$tip.label)) {
  #   return(setNames(c(pc1_delta, pc2, node.age), c("pc1_delta", "pc2_value", "node_age")))
  # }
  
  return(setNames(c(pc1_delta, pc2, node.age), c("pc1_delta", "pc2_value", "node_age")))
  
})

output <- as.data.frame(Reduce(rbind, output))
ggplot(output[], aes(y = pc1_delta, x = pc2_value, colour = node_age)) + geom_point(size = 5)




## What about comparing PC2 values with the rate of evolution in PC1?
## This doens't show much, but then again, I'm not sure that this matters
## Two ways of looking at this: 1) transitions are faster through a crepuscular bridge, 2) transtision must go through a bridge (regardless of speed)
## It doesn't seem like transitions are faster through teh bridge, but that maybe they must go through it (since we don't have strong evidence for cathemeral species)

pc1_mr[[1]]$sig2[rownames(node.data)]

plot(node.data$PC2_mrBM, pc1_mr[[1]]$sig2[rownames(node.data)])
plot(node.data$PC2_mrBM3, pc1_mr[[3]]$sig2[rownames(node.data)])

plot(node.data$PC1_mrBM2, pc1_mr[[2]]$sig2[rownames(node.data)])










node.data$rate <- rates[match(node.data$node, names(rates))]

age_vs_pc1 <- ggplot(node.data, aes(x = node.age, y = PC1_mrBM)) + geom_point()
age_vs_pc1_2 <- ggplot(node.data, aes(x = node.age, y = PC1_mrBM2)) + geom_point()
age_vs_brown <- ggplot(node.data, aes(x = node.age, y = Day_vs_Night_activity)) + geom_point()
age_vs_rates <- ggplot(node.data, aes(x = node.age, y = rates)) + geom_point()

name_vector <- node.data$PC1_mrBM
names(name_vector) <- row.names(node.data)
twoD_pc1 <- ggtree(subset, yscale = "label", yscale_mapping = name_vector)

name_vector <- node.data$PC1_mrBM2
names(name_vector) <- row.names(node.data)
twoD_pc1_2 <- ggtree(subset, yscale = "label", yscale_mapping = name_vector)

name_vector <- node.data$Day_vs_Night_activity
names(name_vector) <- row.names(node.data)
twoD_brown <- ggtree(subset, yscale = "label", yscale_mapping = name_vector)

twoD_brown + twoD_pc1 + twoD_pc1_2 + plot_layout(ncol = 1)

## From this comparison (BM vs variable rates BM adjusted version), that it seems like there is some early burst (into diurnal and nocturnal), followed by late bursts?

## OK, can I derive the disparity? I first need to bin the data, then calculate the disparity per bin, and can plot it, and look at the slope? Lol how did they do that??

brown_disp <- lapply(seq_along(1:nrow(node.data)), function(i) {
  first <- i-1
  last <- i+1
  node.data2 <- node.data[order(node.data$node.age),]
  brown <- max(node.data2$Day_vs_Night_activity[1:i]) - min(node.data2$Day_vs_Night_activity[1:i])
  return(brown)
})

mrBM_disp <- lapply(seq_along(1:nrow(node.data)), function(i) {
  first <- i-1
  last <- i+1
  node.data2 <- node.data[order(node.data$node.age),]
  mrBM <- max(node.data2$PC1_mrBM[1:i]) - min(node.data2$PC1_mrBM[1:i])
  return(mrBM)
})

mrBM2_disp <- lapply(seq_along(1:nrow(node.data)), function(i) {
  first <- i-1
  last <- i+1
  node.data2 <- node.data[order(node.data$node.age),]
  mrBM2 <- max(node.data2$PC1_mrBM2[1:i]) - min(node.data2$PC1_mrBM2[1:i])
  return(mrBM2)
})

df <- data.frame(brown_disp = unlist(brown_disp), mrBM_disp = unlist(mrBM_disp), mrBM2_disp = unlist(mrBM2_disp), node.age = node.data$node.age[order(node.data$node.age)], group = "group")
df$mrBM_dif <- df$mrBM_disp-df$brown_disp
df$mrBM2_dif <- df$mrBM2_disp-df$brown_disp



### Can I run the simulations under brownian motion? If I can't compare slopes, I can at least do what I did above

simulation_pc1 <- replicate(500, rTraitCont(phy = subset, model = "BM", sigma2 = pc1_mr[[1]]$opt$sigsq, root.value = fitBM_pc1$ace[1], ancestor = T))
simulation_pc2 <- replicate(500, rTraitCont(phy = subset, model = "BM", sigma2 = pc2_mr[[1]]$opt$sigsq, root.value = fitBM_pc1$ace[1], ancestor = T))

# I think the problem is, that I don't care about the average values across simulations - these are nonsense
# I care about the average disparity increase, which needs to be calculated per simulation (there is no reason why Boumic always becomes diurnal or nocturnal, but the average should look like the ancestral state!)

row.names(simulation_pc1) <- c(names(sim_avg)[1:60], 61:(60+59))

node.data2 <- node.data[order(node.data$node.age),]

simulation_pc1 <- simulation_pc1[row.names(node.data2),]

## need to run this across all columns separately

brown_sim_disp <- lapply(seq_along(1:ncol(simulation_pc1)), function(y) {
  
  disp_col <- lapply(seq_along(1:nrow(simulation_pc1)), function(i) {
    disp <- max(simulation_pc1[1:i,y]) - min(simulation_pc1[1:i,y])
    return(disp)
  })
  
  return(unlist(disp_col))
})

brown_sim_disp <- Reduce(cbind, brown_sim_disp)
row.names(brown_sim_disp) <- row.names(node.data2)

node.data2$brown_sim_disp <- rowMeans(brown_sim_disp)


df$brown_sim_disp <- rowMeans(brown_sim_disp)

df$mrBM_dif2 <- df$mrBM_disp - df$brown_sim_disp
df$mrBM2_dif2 <- df$mrBM2_disp - df$brown_sim_disp



brown_disp_plot <- ggplot(df, aes(x = node.age, y = brown_disp, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))
brow_sim_disp_plot <- ggplot(df, aes(x = node.age, y = brown_sim_disp, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))
mrBM_disp_plot <- ggplot(df, aes(x = node.age, y = mrBM2_disp, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))
mrBM2_disp_plot <- ggplot(df, aes(x = node.age, y = mrBM_disp, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))
mrBM_dif_plot <- ggplot(df, aes(x = node.age, y = mrBM2_dif, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))
mrBM2_dif_plot <- ggplot(df, aes(x = node.age, y = mrBM_dif, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))

## This is comparing it to the brownian model for ancestral reconstruction, but in Ronco et al, they compare it to 500 simulations under brownian evolution
## Not sure why the above, and not this, but it wouldn't be difficult to run those simulations
brown_disp_plot + mrBM_disp_plot + mrBM2_disp_plot + plot_layout(nrow = 3)
mrBM_dif_plot + mrBM2_dif_plot + plot_layout(nrow = 2)

mrBM_dif_plot2 <- ggplot(df, aes(x = node.age, y = mrBM_dif2, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))
mrBM2_dif_plot2 <- ggplot(df, aes(x = node.age, y = mrBM2_dif2, group = group)) + geom_point() + geom_line() + scale_x_continuous(breaks = round(seq(min(df$node.age), max(df$node.age), by = 0.5),1))


### OK, so the way Fabrizia did it, was that she split the data into 0.15 million year bins, then calculated the slope for sliding windows of 3-consecutive timeslices, and compared those to each other
### She also does this for each of the 500 simulations, and takes the mean of the slope differences (across these 500 comparisons)
### Should be relatively easy to do, but this script is becoming a mess

## OK, so maybe run this # timepoints - 3
## First need to bin the data, then calculate the slope, then calculate the diff, and return the dif

## Make a new node data df
nnd <- node.data2[,c(1,2,3,4,5,6,7)]
nnd$bin <- cut_interval(nnd$node.age, length = 0.15)

test <- sapply(strsplit(as.character(nnd$bin),","), `[`, 2)
test <- as.numeric(str_remove(test, "]"))
nnd$bin <- test
## use the simulations
simulation <- simulation
# simulation$bin <- nnd$bin

# For each bin, calculate the diff between min/max values
nnd_mrBM <- lapply(unique(nnd$bin), function(i) {
  
  disp <- max(nnd$PC1_mrBM[nnd$bin == i]) - min(nnd$PC1_mrBM[nnd$bin == i])
  return(disp)
  
})
nnd_mrBM <- cumsum(unlist(nnd_mrBM))
names(nnd_mrBM) <- unique(nnd$bin)

## Do the same for the simulation
brown_sim_disp <- lapply(seq_along(1:ncol(simulation)), function(y) {
  
  disp_col <- lapply(unique(nnd$bin), function(i) {
    disp <- max(simulation[nnd$bin == i,y]) - min(simulation[nnd$bin == i,y])
    return(disp)
  })
  disp_col <- cumsum(unlist(disp_col))
  names(disp_col) <- unique(nnd$bin)
  return(disp_col)
})

plot(nnd_mrBM)

## OK, now I need to seq along the timescales, and generate the slopes and compare them, for each of the simulations

nnd_mrBM <- data.frame(node.age = as.numeric(names(nnd_mrBM)), mrBM = nnd_mrBM)

brown_sim_disp <- lapply(brown_sim_disp, function(x) {
  
  df <- data.frame(node.age = as.numeric(names(x)), sim_disp = x)
  
  return(df)
})


ggplot(nnd_mrBM, aes(x = node.age, y = mrBM)) + geom_point()
ggplot(brown_sim_disp[[1]], aes(x = node.age, y = sim_disp)) + geom_point()

## Calculate slopes
## This seems to work, but

smo=3

Dif_slopeBMrange = matrix(NA,(length(unique(nnd$bin))-smo), 500 )

for ( i in 1:length(brown_sim_disp)) {
  
  for ( j in 1:(length(unique(nnd$bin))-smo) ){
    
    tmpSIM_mrBM <- nnd_mrBM[j:(j+smo),]
    tmpSIM_mrBM$mrBM <- tmpSIM_mrBM$mrBM/max(max(tmpSIM_mrBM$mrBM), 0.0000000000000001)
    mrBM_slope <- lm(tmpSIM_mrBM$mrBM ~ tmpSIM_mrBM$node.age)$coefficients[2]
    
    tmpSIM_brown <- brown_sim_disp[[i]][j:(j+smo),]
    tmpSIM_brown$sim_disp <- tmpSIM_brown$sim_disp/max(max(tmpSIM_brown$sim_disp), 0.0000000000000001)
    brown_slope <- lm(tmpSIM_brown$sim_disp ~ tmpSIM_brown$node.age)$coefficients[2]
    # tmpSIM= cbind(statsOtime$time, normalize(rangeSim[,i+1])*100)
    #BM_slope= lm( tmpSIM[j:(j+smo),2] ~ tmpSIM[j:(j+smo),1])$coefficients[2]
    # OBS_slope= lm( LstatsOtime[j:(j+smo),2] ~ LstatsOtime[j:(j+smo),1])$coefficients[2]
    
    Dif_slopeBMrange[j,i] <- mrBM_slope - brown_slope
    
  }
  
}

MEAN_DIF_SLOPE_RANGE = c(apply(Dif_slopeBMrange, 1, mean), NA, NA, NA)
STDEV_DIF_SLOPE_RANGE = c(apply(Dif_slopeBMrange, 1, sd), NA, NA, NA)

new_df <- data.frame(mdsr = MEAN_DIF_SLOPE_RANGE, node.age = (unique(nnd$bin)-max(nnd$bin))*-1, group = "group", stdev_pos = MEAN_DIF_SLOPE_RANGE + STDEV_DIF_SLOPE_RANGE, stdev_neg = MEAN_DIF_SLOPE_RANGE - STDEV_DIF_SLOPE_RANGE)
new_df2 <- pivot_longer(new_df, cols = c(mdsr, stdev_pos, stdev_neg))

ggplot(new_df, aes(x = node.age, y = mdsr, group = group, ymin = stdev_neg, ymax = stdev_pos)) + geom_point() + geom_line() + geom_ribbon(alpha = 0.25) + scale_x_reverse()
ggplot(new_df2, aes(x = node.age, y = value, group = name, colour = name)) + geom_point() + geom_line()




################################################################################################################################################################################################
#### Identify nodes for tribes  ################################################################################################################################################################
################################################################################################################################################################################################

ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab() + geom_text(aes(label = node))

lamps <- 63
ectos <- 104
trops <- 114
cyps <- 110
erets <- 52
limnos <- 42
cyphos <- 41
bouls <- 61
haplos <- 53

################################################################################################################################################################################################
#### Make Plots  ###############################################################################################################################################################################
################################################################################################################################################################################################

phylo.plot <- ggtree(as.phylo(subset), layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 5)

p3 <- phylo.plot + geom_tile(data = cichlids.diel, aes(y=y, x=x, fill = day_night_dif), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu") + scale_color_discrete(na.value = 'transparent')


# Plot and save the tree with ancestral reconstruction colouring the edges (which are matched to a node)
#ancestral.plot <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Day_vs_Night_activity) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Day_vs_Night_activity), shape = 16, size = 6) + scale_color_distiller(palette = "YlOrRd", direction = 1)


##### THIS NEEDS TO BE FIXED SO THE SCALES ARE CENTERED #######

ancestral.plot.day_night_dif <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Day_vs_Night_activity) + scale_color_distiller(palette = "RdBu") + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Day_vs_Night_activity), shape = 16, size = 6) + scale_color_distiller(palette = "RdBu", direction = 1)
ancestral.plot.total_rest <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Total_Rest) + scale_color_distiller(palette = "PRGn", direction = 1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Total_Rest), shape = 16, size = 6) + scale_color_distiller(palette = "PRGn", direction = 1)
ancestral.plot.peak <- ggtree(subset, layout = "circular", size = 2) %<+% node.data + aes(color = Peak_Percentage) + scale_color_distiller(palette = "BrBG", direction = 1) + geom_tiplab(color = "black", size = 6, offset = 0.5) + geom_tippoint(aes(color = Peak_Percentage), shape = 16, size = 6) + scale_color_distiller(palette = "BrBG", direction = -1)

ancestral.plot <- ancestral.plot.day_night_dif + ancestral.plot.total_rest + ancestral.plot.peak + plot_layout(nrow = 2)

################################################################################################################################################################################################
#### Save Plots  ###############################################################################################################################################################################
################################################################################################################################################################################################

pdf("Ancestral_phylograms_NEW.pdf", width = 20, height = 15)
ancestral.plot
dev.off()

png("Ancestral_phylograms_NEW.png", width = 20, height = 15, units = "in", res = 500)
ancestral.plot
dev.off()


################################################################################################################################################################################################
#### Plot Phenotypes on Big Tree  ##############################################################################################################################################################
################################################################################################################################################################################################

phylo.plot <- ggtree(as.phylo(lt_phylo), layout = "circular", right = T) + theme_tree(bgcolor = NA) + geom_tiplab(offset = 3, size = 2.5)

cichlids.diel.2 <- data.frame(X = 1:length(lt_phylo$tip.label), six_letter_name_Ronco = lt_phylo$tip.label, total_rest = NA, peak_amplitude = NA, peak = NA, day_night_dif = NA, cluster = NA, tips = as.character(lt_phylo$tip.label), x = 10.5, y = NA)
cichlids.diel.2$tips <- factor(cichlids.diel.2$tips, levels = rev(with(phylo.plot$data, label[order(y, decreasing=T)])))
cichlids.diel.2$y <- as.numeric(cichlids.diel.2$tips)
cichlids.diel.2$day_night_dif <- cichlids.diel$day_night_dif[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]
cichlids.diel.2$total_rest <- cichlids.diel$total_rest[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]
cichlids.diel.2$peak <- cichlids.diel$peak[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]
cichlids.diel.2$peak_amplitude <- cichlids.diel$peak_amplitude[match(cichlids.diel.2$tips, cichlids.diel$tips, nomatch = NA)]

head(cichlids.diel.2)

p4 <- phylo.plot + geom_tile(data = cichlids.diel.2, aes(y=y, x=x, fill = day_night_dif), width = 1, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "RdBu", na.value = 'transparent', limits = c(-0.75,0.75))
p4 <- p4 + new_scale("size") + new_scale("colour") + geom_point(data = cichlids.diel.2, aes(y=y, x=x+1.25, size = total_rest, colour = total_rest)) + scale_colour_distiller(palette = "PRGn", direction = 1, na.value = 'transparent') 
p4 <- p4 + new_scale("colour") + geom_point(data = cichlids.diel.2, aes(y=y, x=x+2, colour = peak, size = 10), shape = 18) + scale_colour_distiller(palette = "BrBG", direction = 1, na.value = 'transparent', limits = c(0,1))
# p4 <- p4 + new_scale("fill") + geom_tile(data = cichlids.diel.2, aes(y=y, x=x+3, fill = peak), width = 2, inherit.aes = FALSE, color = "white") + scale_fill_distiller(palette = "YlOrRd")
# p4 <- p4 + xlim(c(0,17))

pdf("Ancestral_phylograms_all_species.pdf", width = 12, height = 12)
p4
dev.off()

png("Ancestral_phylograms_all_species.png", width = 12, height = 12, units = "in", res = 750)
p4
dev.off()

################################################################################################################################################################################################
#### Plot ancestral nodes + extant nodes  ######################################################################################################################################################
################################################################################################################################################################################################

pc1 <- c(cichlids.diel$pc1, rowMeans(fitBM_day_night_dif$CI95))
pc2 <- c(cichlids.diel$pc2, rowMeans(fitBM_peak$CI95))
df <- data.frame(pc1 = pc1, pc2 = pc2, nodes = c(cichlids.diel$tips, row.names(fitBM_day_night_dif$CI95)))

ggplot(df, aes(x = pc1, y = pc2, colour = nodes)) + geom_point()

## This works, and seems interesting, but it is only the BM model, which may not be likely
## I think I'm better off looking at early vs late bursts and comparing those models

################################################################################################################################################################################################
#### Test mvSLOUCH for multivariate  ###########################################################################################################################################################
################################################################################################################################################################################################

subset <- keep.tip(lt_phylo, tip = cichlids.diel$tips[cichlids.diel$tips %in% lt_phylo$tip.label])

cichlid.data <- as.matrix(cichlids.diel[,c(2:3,22)])
row.names(cichlid.data) <- cichlids.diel$tips
cichlid.data <- cichlid.data[row.names(cichlid.data) %in% subset$tip.label,]


BMestim <- mvSLOUCH::BrownianMotionModel(subset, cichlid.data)

# Ancestral state values
# These are really close to the confidence intervals from the individual BM ace models - not sure if that's all this is doing at the moment
BMestim$ParamsInModel$vX0

OU1BM <- mvSLOUCH::mvslouchModel(subset, cichlid.data, kY = 2, predictors = c(2))

OU1OU <- mvSLOUCH::ouchModel(subset, cichlid.data, predictors = c(2))

cbind(BMestim$ParamSummary$dof,
      OU1BM$FinalFound$ParamSummary$dof,
      OU1OU$FinalFound$ParamSummary$dof)

## OK so I think this isn't doing what I want, or not testing what I want to test
## This is asking whether or not the optima of one trait influences the optima of another - I don't think that's what I'm looking for
## My traits are not independent - they are dependent, which is weird and different
## As far as I can tell from this, there is evidence that the optima don't influence each other, but IDK

################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################

# crepuscular_ids <- c("Neosav", "Neopul", "Neooli", "Neokom", "Neohel", "Neobri")
# cathemeral_ids <- c("Neogra", "Neocyg",  "Neomar", "NeofaM", "Neocra", "Neofal")
# 
# fv3 <- cichlids.diel[cichlids.diel$six_letter_name_Ronco %in% c(crepuscular_ids, cathemeral_ids), c("six_letter_name_Ronco", "day_night_dif", "total_rest")]
# 
# isotope_data_3 <- isotope_data[isotope_data$SpeciesID %in% c(crepuscular_ids, cathemeral_ids),]
# isotope_data_3$type <- ifelse(isotope_data_3$SpeciesID %in% crepuscular_ids, "crepuscular", "non-crepuscular")
# 
# ggplot(isotope_data_3, aes(x = type, y = d15N)) + geom_jitter() + geom_boxplot() + theme_classic()

