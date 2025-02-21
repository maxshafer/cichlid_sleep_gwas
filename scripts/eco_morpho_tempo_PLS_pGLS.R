## OK, I need to run the pgls for all possible comparisons of traits
## Alternatively, I can do what Fabrizia did, and run a PLS, then a pgls on the loadings
## PLS is like PCA, but instead finds the loading that maximizes the relaionship between two samples (finds the relationship between two matrices)
## For this I think I would need the raw body/LPJ/isotope data, not the PC values
## I could then run PLS on this vs the same data Annika used to get PC values

## Ideally I would generate a df, with a column for responder and a column for predictor variable, and columns for Rsquared, pvalue, etc, which could be indexed or plotted (heatmap of pvalues)

library(caper)
library(pls)
library(geomorph)
library(ggplot2)
library(gsheet)
library(ggrepel)
library(dplyr)
library(patchwork)
library(viridis)
library(phytools)
library(ggtree)
library(ggnewscale)
library(ggrain)
library(geiger)

theme_set(theme_minimal())
theme_update(axis.text = element_text(colour = "black"))
theme_set(theme_classic())
theme_update(axis.text = element_text(colour = "black"))

########################################################################################################
################ Load and prep input data    ###########################################################
########################################################################################################

## Meta data has tribe, diet, habitat info
url <- 'https://docs.google.com/spreadsheets/d/1m6cV5W-UjDdIA6_ladE0H4iCwQU_xcbWwnQ5w2-kgbM/edit?usp=sharing'
meta_data <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

## Import raw normalized speed data and transpose it (so that row names are species)
speed_data <- read.csv("pca_input_zscore.csv", row.names = "time_of_day")
speed_data <- t(speed_data)

sizes <- read.csv("pheno_data/experimental_fish_sizes.csv")

## load in old and new data, so I can have total_rest, as well as Annika's pc values
sleepData <- read.csv("pheno_data/cichlid_pc-loadings_eco-morph_rest_full.csv")
sleepData$tribe <- meta_data$tribe[match(sleepData$six_letter_name_Ronco, meta_data$six_letter_name_Ronco)]
sleepData$max_size <- meta_data$size_male[match(sleepData$six_letter_name_Ronco, meta_data$six_letter_name_Ronco)]
sleepData$category <- meta_data$category_paper[match(sleepData$six_letter_name_Ronco, meta_data$six_letter_name_Ronco)]
sleepData$size_mean <- sizes$size_mean[match(sleepData$six_letter_name_Ronco, sizes$six_letter_name_Ronco)]
sleepData$size_max <- sizes$size_max[match(sleepData$six_letter_name_Ronco, sizes$six_letter_name_Ronco)]

rownames(sleepData) <- sleepData$six_letter_name_Ronco

## Load and trim tree
phylo_lt <- read.nexus("scripts/05_BEAST_RAxML.tre")
cichlidTree <- keep.tip(phylo_lt, tip = sleepData$six_letter_name_Ronco[sleepData$six_letter_name_Ronco %in% phylo_lt$tip.label])

## Double trim with input data
sleepData <- sleepData[sleepData$six_letter_name_Ronco %in% cichlidTree$tip.label,]
speed_data <- speed_data[row.names(speed_data) %in% cichlidTree$tip.label,]

########################################################################################################
################ Make tree figures    ##################################################################
########################################################################################################

phylo.plot <- ggtree(as.phylo(cichlidTree), layout = "rectangular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 6)

sleep_data <- sleepData
sleep_data$x <- 11
sleep_data$y <- phylo.plot$data$y[match(sleep_data$six_letter_name_Ronco, phylo.plot$data$label)]
sleep_data$diet <- meta_data$diet[match(sleep_data$six_letter_name_Ronco, meta_data$six_letter_name_Ronco)]

phylo.plot <- ggtree(as.phylo(cichlidTree), layout = "rectangular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 0)

# crepuscular': '#26D97A', 'nocturnal': '#40A9BF', 'diurnal': '#CED926', 'cathemeral': '#737F8C'

#geom_point(data = cichlids.diel.2, aes(y=y, x=x+2, colour = peak, size = 10), shape = 18)

p3 <- phylo.plot 
p3 <- p3 + geom_point(data = sleep_data, aes(colour = category, y = y, x = x+1.5), size = 2) + scale_colour_manual(values = c('#737F8C', '#26D97A', '#CED926', '#40A9BF'))
p3 <- p3 + geom_tile(data = sleep_data, aes(y=y, x=x+2.75, fill = pc1), width = 1, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", direction = 1) #+ scale_color_discrete(na.value = 'transparent')
p3 <- p3 + new_scale("fill") + new_scale("size") + geom_tile(data = sleep_data, aes(y=y, x=x+4, fill = pc2), width = 1, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "BrBG", direction = -1)
p3 <- p3 + new_scale("size") + new_scale("fill") + geom_tile(data = sleep_data, aes(y=y, x=x+5.75, fill = total_rest, width = total_rest/9), height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "PRGn", direction = -1) 

p4 <- p3 + plot_layout(guides = "collect", width = unit(c(100), "mm"), height = unit(c(100), "mm"))

########################################################################################################
################ Make quadratic model figure    ########################################################
########################################################################################################

#create a new variable for pc1_2
sleep_data$pc1_2 <- sleep_data$pc1^2
sleep_data$pc2_2 <- sleep_data$pc2^2
sleep_data$total_rest_2 <- sleep_data$total_rest^2

#fit quadratic regression model
quadraticModel <- lm(pc2 ~ pc1 + pc1_2, data=sleep_data) # p-value: 2.245e-11
# quadraticModel <- lm(pc2 ~ total_rest + total_rest_2, data=sleep_data) # p-value: 0.2893
# quadraticModel <- lm(pc1 ~ pc2 + pc2_2, data=sleep_data) # p-value: 0.5013
# quadraticModel <- lm(pc1 ~ total_rest + total_rest_2, data=sleep_data) # p-value: 0.663
# quadraticModel <- lm(total_rest ~ pc1 + pc1_2, data=sleep_data) # p-value: 0.5532
# quadraticModel <- lm(total_rest ~ pc2 + pc2_2, data=sleep_data) # p-value: 0.7891

# create sequence of hour values
Values <- seq(-0.2, 0.2, 0.001)

## Make a trend line for the quadratic fit
#create list of predicted happines levels using quadratic model
Predict <- predict(quadraticModel,list(pc1=Values, pc1_2=Values^2))

# PC2 = -4.626(pc1)^2 + -0.026(hours) - 0.191

## Run a PGLS
compD = comparative.data(phy = cichlidTree, data=sleep_data, six_letter_name_Ronco)
fitmod <- pgls(pc2 ~ pc1 + pc1_2, compD)
summary_fitmod <- summary(fitmod)

# PC2 = -4.083814(pc1)^2 + 0.035057(hours) - 0.186

Predict_pgls <- predict(fitmod,list(pc1=Values, pc1_2=Values^2))

#create scatterplot of original data values
trend_line_data <- data.frame(pc1 = Values, pc2 = Predict_pgls)

#add predicted lines based on quadratic regression model
sleep_data$tribe <- factor(sleep_data$tribe, levels = c("Boulengerochromini", "Cyphotilapiini", "Cyprichromini", "Ectodini", "Eretmodini", "Lamprologini", "Limnochromini", "Trophenini"))
pca_plot <- ggplot(sleep_data, aes(x = pc1, y = pc2, colour = tribe, label = six_letter_name_Ronco)) + geom_point(size = 4) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("PC2 (crepuscular preference)") + xlab("PC1 (diurnal - nocturnal preference)") #+ geom_text_repel(force = 2) 
pca_plot <- pca_plot + geom_line(data = trend_line_data, aes(x=pc1, y = pc2), inherit.aes = FALSE, size = 1.5, alpha = 1, colour = "black")

## Update #s based on summary(quandraticModel) and summary_fitmod)
pc1_vs_pc2_plot <- pca_plot + annotate("text", x=-0.1, y=0.01, label= "Quadratic \n r^2 = 0.6334 \n p-value = 3.787e-13 \n Quadratic (PGLS) \n r^2 = 0.5302 \n p-value = 2.023e-09", hjust = 0)


########################################################################################################
################ Make diet figures    ##################################################################
########################################################################################################

# 'Invertivore': 'tomato', 'Piscivore': 'steelblue', 'Zooplanktivore': 'sandybrown', 'Algivore': 'mediumseagreen'
diet.colours <- c('Invertivore' = "#EF4026", 'Piscivore' = "#4682b4", 'Zooplanktivore' = "#f4a460", 'Algivore' = "#3cb371") 

sleep_data$diet <- factor(sleep_data$diet, levels = c("Invertivore", "Piscivore", "Zooplanktivore", "Algivore"))
tr_box <- ggplot(sleep_data, aes(x = diet, y = total_rest)) + geom_boxplot(alpha = 0.25, colour = "black", fill = "grey", outlier.colour = "transparent") + geom_jitter(colour = "black", width = 0.15) + xlab("") + theme_classic() + coord_flip() #+ scale_colour_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371")) + scale_fill_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371"))
pc1_box <- ggplot(sleep_data, aes(x = diet, y = pc1)) + geom_boxplot(alpha = 0.25, colour = "black", fill = "grey", outlier.colour = "transparent") + geom_jitter(colour = "black", width = 0.15) + xlab("") + theme_classic() + coord_flip() #+ scale_colour_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371")) + scale_fill_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371"))
pc2_box <- ggplot(sleep_data, aes(x = diet, y = pc2)) + geom_boxplot(alpha = 0.25, colour = "black", fill = "grey", outlier.colour = "transparent") + geom_jitter(colour = "black", width = 0.15) + xlab("") + theme_classic() + coord_flip() #+ scale_colour_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371")) + scale_fill_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371"))

### Below is for doing ANOVAs, but there is no phylogenetically significant result
### I can use either phytools or geiger's implementation, which only differ because they run simulations
### I'll use gieger, because it gives me both the non-phylogenetically correct and phylo corrected pvalues, which phytools doens't

test <- aov(total_rest ~ diet, sleep_data)
summary.lm(test)

sleep_data_diet <- sleep_data$diet[match(cichlidTree$tip.label, sleep_data$six_letter_name_Ronco)]
names(sleep_data_diet) <- cichlidTree$tip.label
sleep_data_pc1 <- sleep_data$pc1[match(cichlidTree$tip.label, sleep_data$six_letter_name_Ronco)]
names(sleep_data_pc1) <- cichlidTree$tip.label
sleep_data_pc2 <- sleep_data$pc2[match(cichlidTree$tip.label, sleep_data$six_letter_name_Ronco)]
names(sleep_data_pc2) <- cichlidTree$tip.label
sleep_data_tr <- sleep_data$total_rest[match(cichlidTree$tip.label, sleep_data$six_letter_name_Ronco)]
names(sleep_data_tr) <- cichlidTree$tip.label
sleep_data_d13c <- sleep_data$d13C[match(cichlidTree$tip.label, sleep_data$six_letter_name_Ronco)]
names(sleep_data_d13c) <- cichlidTree$tip.label

summary(aov.phylo(sleep_data_pc1 ~ sleep_data_diet, cichlidTree, nsim = 1000, test = "Wilks")) # Pr(>F) = 0.86647, Pr(>F) given phy = 0.9451 ## Phylo signal is random, due to the nsim
summary(aov.phylo(sleep_data_pc2 ~ sleep_data_diet, cichlidTree, nsim = 1000, test = "Wilks")) # Pr(>F) = 0.28647, Pr(>F) given phy = 0.6344
summary(aov.phylo(sleep_data_tr ~ sleep_data_diet, cichlidTree, nsim = 1000, test = "Wilks")) # Pr(>F) = 0.86112, Pr(>F) given phy = 0.9471

########################################################################################################
################ Compute phylogenetic signal    ########################################################
########################################################################################################

## Pagel’s lambda is a measure of phylogenetic ‘signal’ in which the degree to which shared history of taxa 
## has driven trait distributions at tips. In this model, internal branch lengths are transformed by the lambda 
## parameter value. When the parameter lambda equals 1, branches are transformed by multiplying by 1 and so 
## the model is equal to Brownian motion (high phylogenetic signal). Values of lambda under 1 suggest there has 
## been less influence of shared history on trait values at the tips. Finally, a lambda value of 0 indicates no 
## phylogenetic influence on trait distributions, and is equivalent to a ‘star phylogeny’ with no shared branch lengths.

trait.data.pc1 <- sleepData$pc1[match(cichlidTree$tip.label, sleepData$six_letter_name_Ronco)]
names(trait.data.pc1) <- cichlidTree$tip.label

pc1_l <-phylosig(cichlidTree, trait.data.pc1, method = "lambda", test = T) # 0.630368
pc1_K <- phylosig(cichlidTree, trait.data.pc1, method = "K", test = T) # 0.47297

trait.data.pc2 <- sleepData$pc2[match(cichlidTree$tip.label, sleepData$six_letter_name_Ronco)]
names(trait.data.pc2) <- cichlidTree$tip.label

pc2_l <- phylosig(cichlidTree, trait.data.pc2, method = "lambda", test = T) # 6.6107e-05
pc2_K <- phylosig(cichlidTree, trait.data.pc2, method = "K", test = T) # 0.438309

trait.data.tr <- sleepData$total_rest[match(cichlidTree$tip.label, sleepData$six_letter_name_Ronco)]
names(trait.data.tr) <- cichlidTree$tip.label

tr_l <- phylosig(cichlidTree, trait.data.tr, method = "lambda", test = T) # 0.608944
tr_K <- phylosig(cichlidTree, trait.data.tr, method = "K", test = T) # 0.52516

## pc1 and tr are significantly different from the null hypothesis of 0, pc2 is not (because there is literally no signal)
## I wonder if pc1 and tr and significantly different from 1, but that doesn't mean anything? Either way, these new numbers with more data have lower lambda values

########################################################################################################
################ OK run pairwise pgls on activity features    ##########################################
########################################################################################################

variables <- c("pc1", "pc2", "total_rest", "max_size", "size_mean", "size_max")
comp.data <- comparative.data(cichlidTree, sleepData, six_letter_name_Ronco, vcv = TRUE)

## This gives me the full grid, and also keeps them as the correct class objects
grid <- expand.grid(variables, variables)

output.grid <- apply(grid, 1, function(x) {
  mod <- pgls(as.formula(paste(x[1], "~", x[2])), comp.data, lambda = "ML")
  return(mod)
})

## make the data frame as before

lmp <- function (modelobject) {
  if (class(modelobject) != "pgls") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

out <- data.frame(namex = unlist(lapply(output.grid, function(x) x$varNames[2])),
                  namey = unlist(lapply(output.grid, function(x) x$varNames[1])),
                  r.squared = unlist(lapply(output.grid, function(x) summary(x)$r.squared)),
                  p.value = unlist(lapply(output.grid, function(x) lmp(x))))

out <- out[!(is.na(out$namex)),]

## Plot
pc1_tr <- ggplot(sleep_data, aes(x = pc1, y = total_rest, colour = tribe, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Total rest (hrs)") + xlab("PC1 (diurnal - nocturnal preference)") #+ geom_text_repel(force = 2) 
pc1_tr <- pc1_tr + annotate("text", x=-0.05, y=16, label = paste("Linear (pGLS) \n r^2 = ", round(out$r.squared[out$namex == "pc1" & out$namey == "total_rest"], 4), "\n p-value = ", round(out$p.value[out$namex == "pc1" & out$namey == "total_rest"], 4), sep = ""), hjust = 0)

pc2_tr <- ggplot(sleep_data, aes(x = pc2, y = total_rest, colour = tribe, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Total rest (hrs)") + xlab("PC2 (crepuscular preference)") #+ geom_text_repel(force = 2) 
pc2_tr <- pc2_tr + annotate("text", x=-0.05, y=17, label = paste("Linear (pGLS) \n r^2 = ", round(out$r.squared[out$namex == "pc2" & out$namey == "total_rest"], 4), "\n p-value = ", round(out$p.value[out$namex == "pc2" & out$namey == "total_rest"], 4), sep = ""), hjust = 0)

pc1_pc2 <- ggplot(sleep_data, aes(x = pc1, y = pc2, colour = tribe, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("PC2 (crepuscular preference)") + xlab("PC1 (diurnal - nocturnal preference)") #+ geom_text_repel(force = 2) 
pc1_pc2 <- pc1_pc2 + annotate("text", x=-0.05, y=0, label = paste("Linear (pGLS) \n r^2 = ", round(out$r.squared[out$namex == "pc1" & out$namey == "pc2"], 4), "\n p-value = ", round(out$p.value[out$namex == "pc1" & out$namey == "pc2"], 4), sep = ""), hjust = 0)


###### Compare sizes to total rest
tr_ms <- ggplot(sleep_data, aes(y = total_rest, x = max_size, colour = diet, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Total rest (hrs)") + xlab("Max size (literature)") #+ geom_text_repel(force = 2) 
tr_ms <- tr_ms + annotate("text", x=40, y=10, label = paste("Linear (pGLS) \n r^2 = ", round(out$r.squared[out$namex == "total_rest" & out$namey == "max_size"], 4), "\n p-value = ", round(out$p.value[out$namex == "total_rest" & out$namey == "max_size"], 4), sep = ""), hjust = 0)

tr_sm <- ggplot(sleep_data, aes(y = total_rest, x = size_mean, colour = diet, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Total rest (hrs)") + xlab("Mean size (measured)") #+ geom_text_repel(force = 2) 
tr_sm <- tr_sm + annotate("text", x=6.5, y=15, label = paste("Linear (pGLS) \n r^2 = ", round(out$r.squared[out$namex == "total_rest" & out$namey == "size_mean"], 4), "\n p-value = ", round(out$p.value[out$namex == "total_rest" & out$namey == "size_mean"], 4), sep = ""), hjust = 0)

tr_smax <- ggplot(sleep_data, aes(y = total_rest, x = size_max, colour = diet, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Total rest (hrs)") + xlab("Max size (measured)") #+ geom_text_repel(force = 2) 
tr_smax <- tr_smax + annotate("text", x=10, y=11, label = paste("Linear (pGLS) \n r^2 = ", round(out$r.squared[out$namex == "total_rest" & out$namey == "size_max"], 4), "\n p-value = ", round(out$p.value[out$namex == "total_rest" & out$namey == "size_max"], 4), sep = ""), hjust = 0)

## make plots for comparing sizes
ms_sm <- ggplot(sleep_data, aes(y = max_size, x = size_mean, colour = diet, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Max size (literature)") + xlab("Mean size (measured)") #+ geom_text_repel(force = 2) 
ms_smax <- ggplot(sleep_data, aes(y = max_size, x = size_max, colour = diet, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Max size (literature)") + xlab("Max size (measured)") #+ geom_text_repel(force = 2) 
sm_smax <- ggplot(sleep_data, aes(y = size_mean, x = size_max, colour = diet, label = six_letter_name_Ronco)) + geom_point(size = 2) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("Mean size (measured)") + xlab("Max size (measured)") #+ geom_text_repel(force = 2) 

# tr_ms + tr_sm + tr_smax + ms_sm + ms_smax + sm_smax + plot_layout(nrow = 2, guides = "collect")

########################################################################################################
################ OK run pairwise PLS / pGLS    #########################################################
########################################################################################################

ecomorpho.tree <- drop.tip(cichlidTree, c("Astbur", "Neodev", "Neolon", "Telluf"))


###### OK, so I ran the scripts in ronco_et_al to generate the input data for her PLS analysis, and can now use it against our data!
## Try the geomorph two-block pls (it's somehow like cca?)
## This works, and is not significant, again I worry about the non-correlated nature of the paired PC scores...

## Load and prep body morpho
load("~/Documents/R_Projects/ronco_et_al/trait_evolution/02_Morpho_Eco/body/sp_means_body.Rdata")
body= sp_means_body[,, dimnames(sp_means_body)[[3]] %in% sleepData$six_letter_name_Ronco ]
body= body[,, dimnames(body)[[3]] %in% ecomorpho.tree$tip.label]

## Load and prep UOJ morpho
## This is a subset of the whole body morpho data above
load("~/Documents/R_Projects/ronco_et_al/trait_evolution/02_Morpho_Eco/oral/sp_means_OJ.Rdata")
uoj= sp_means_OJ[,, dimnames(sp_means_OJ)[[3]] %in% sleepData$six_letter_name_Ronco ]
uoj= uoj[,, dimnames(uoj)[[3]] %in% ecomorpho.tree$tip.label]

## Load and prep LPJ morpho
load("~/Documents/R_Projects/ronco_et_al/trait_evolution/02_Morpho_Eco/LPJ/sp_means_LPJ.Rdata")
lpj= sp_means_LPJ[,, dimnames(sp_means_LPJ)[[3]] %in% sleepData$six_letter_name_Ronco ]
lpj= lpj[,, dimnames(lpj)[[3]] %in% ecomorpho.tree$tip.label]

# Prep isotope data
isotopes <- sleepData[,c("d15N", "d13C")]
isotopes <- isotopes[!is.na(isotopes$d13C),]

# Prep speed and total rest data
speed <- as.data.frame(speed_data[row.names(isotopes),])
rest <- sleepData[row.names(isotopes),c("total_rest"), drop = FALSE]
size <- sleepData[row.names(isotopes),c("max_size"), drop = FALSE]

### Calculate multivariate K using geomorph

kmulti <- physignal(A = as.matrix(speed), phy = ecomorpho.tree, iter = 999, print.progress = T) # Kmulti = 0.4114, p-value = 0.006

## Run two block pls for speed and rest vs isotopes, body, uoj, and lpj
variables.pls <- c("speed", "rest", "body", "uoj", "lpj", "isotopes", "size")

## This gives me the full grid, and also keeps them as the correct class objects
grid.pls <- expand.grid(variables.pls, variables.pls)

output.grid.pls <- apply(grid.pls, 1, function(x) {
  
  if (x[1] == x[2]) {
    return ("NA")
  } 
  
  if ("rest" %in% x & "size" %in% x) {
    return ("NA")
  } else {

    f <- two.b.pls(eval(as.name(paste(x[1]))), eval(as.name(paste(x[2]))))
    f1 <- data.frame( "XScores"=f$ XScores[,1], "YScores"=f$ YScores[,1] , "species"= row.names(speed))
    
    compset = comparative.data(ecomorpho.tree, f1, species)
    # Response ~ Predictor
    mod <- pgls(YScores ~ XScores, compset, lambda= 'ML')
    return(mod)
  }
  
})

# ## To get deformed grids for body shape that correlates with high or low rest amounts
# pls_plot <- plot(output.grid.pls[[10]])
# picknplot.shape(pls_plot)
body_morph_data <- output.grid.pls[[10]]$data$data
body_morph_data$total_rest <- sleep_data$total_rest[match(rownames(body_morph_data), sleep_data$six_letter_name_Ronco)]
body_morph_data$tribe <- sleep_data$tribe[match(rownames(body_morph_data), sleep_data$six_letter_name_Ronco)]
body_morph_plot <- ggplot(body_morph_data, aes(x = total_rest, y = XScores, colour = tribe)) + geom_point() + xlab("Total Rest (hrs)") + ylab("Body morphology PLS1 values") + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773"))

output.grid.pls <- output.grid.pls[!(is.na(output.grid.pls))]

out.pls <- data.frame(namex = grid.pls$Var1,
                  namey = grid.pls$Var2,
                  r.squared = unlist(lapply(output.grid.pls, function(x) ifelse(x[1] == "NA", NA, summary(x)$r.squared))),
                  adj.r.squared = unlist(lapply(output.grid.pls, function(x) ifelse(x[1] == "NA", NA, summary(x)$adj.r.squared))),
                  p.value = unlist(lapply(output.grid.pls, function(x) ifelse(x[1] == "NA", NA, lmp(x)))))

out.pls$Bonf.pvalue <- out.pls$p.value*nrow(out.pls)

out.pls$label <- ifelse(out.pls$Bonf.pvalue < 0.0001, "****", ifelse(out.pls$Bonf.pvalue < 0.001, "***", ifelse(out.pls$Bonf.pvalue < 0.01, "**", ifelse(out.pls$Bonf.pvalue < 0.05, "*", "")) ) )

# Maybe remove body vs uoj, since they have overlapping data (16 and 21)
out.pls[out.pls$namex == "uoj" & out.pls$namey == "body",c(3:5)] <- NA
out.pls[out.pls$namex == "body" & out.pls$namey == "uoj",c(3:5)] <- NA

## Plot
# Y ~ X, meaning, does X predict Y
pls.grid.plot.arsq <- ggplot(out.pls, aes(x = namex, y = namey, fill = adj.r.squared)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = paste(round(adj.r.squared, 2), "\n", label, sep = "")), colour = "white")


##### Redo, to make plots
output.grid.data <- apply(grid.pls, 1, function(x) {
  
  if (x[1] == x[2]) {
    return (NA)
  } 
  
  if ("rest" %in% x & "size" %in% x) {
    return ("NA")
  } else {
    
    f <- two.b.pls(eval(as.name(paste(x[1]))), eval(as.name(paste(x[2]))))
    f1 <- data.frame( "XScores"=f$ XScores[,1], "YScores"=f$ YScores[,1] , "species"= row.names(speed))
    
    f1$tribe <- sleepData$tribe[match(row.names(isotopes), sleepData$six_letter_name_Ronco)]
    f1$d13C <- sleepData$d13C[match(row.names(isotopes), sleepData$six_letter_name_Ronco)]
    f1$d15N <- sleepData$d15N[match(row.names(isotopes), sleepData$six_letter_name_Ronco)]
    f1$pc1 <- sleepData$pc1[match(row.names(isotopes), sleepData$six_letter_name_Ronco)]
    f1$tr <- sleepData$total_rest[match(row.names(isotopes), sleepData$six_letter_name_Ronco)]
    f1$comp <- paste(x[1], x[2], sep = "~")
    
    return(f1)
  }
  
})

## Then do a lapply to make ggplots
tribe.plots <- lapply(output.grid.data, function(x) {

  if (is.na(x)) {
    plot <- ggplot() + theme_void()
  } else {
    plot <- ggplot(x, aes(x=XScores, y=YScores, colour = tribe)) + geom_point(size = 0.25) + theme_classic() + ggtitle(label = x$comp[1]) + theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8), title = element_text(size = 8))
    plot <- plot + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773"))
    return(plot)
  }

})

ggplotRegression <- function(fit, data) {
  
  plot <- ggplot(data, aes(x = XScores, y = YScores, colour = tribe)) + geom_point() + theme_classic() + stat_smooth(method = "glm", col = "red")
  plot + labs(title = paste("Lambda =", round(summary(fit)$param.CI$lambda$opt, 3),
                              "Adj R2 = ",round(summary(fit)$adj.r.squared, 2),
                              " P =",round(summary(fit)$coef[2,4], 7)))
         
}
  
tribe.plots.lm <- lapply(seq_along(output.grid.data), function(x) {
  
  if (is.na(output.grid.data[[x]])) {
    plot <- ggplot() + theme_void()
  } else {
    plot <- ggplotRegression(fit = output.grid.pls[[x]], data = output.grid.data[[x]])
    return(plot)
  }
  
  
})


########################################################################################################
################ Save out plots    #####################################################################
########################################################################################################

pdf("Extended_data_Fig5.pdf", width = 15, height = 15)
(wrap_plots(tribe.plots) + plot_layout(guides = "collect", width = unit(c(15), "mm"), height = unit(c(15), "mm"))) / (pc1_pc2 + pc1_tr + pc2_tr + plot_layout(guides = "collect", width = unit(c(45), "mm"), height = unit(c(45), "mm")))
dev.off()


# design = "
# AAAAAABE
# AAAAAACF
# AAAAAADG
# HHHIIIJJ
# HHHIIIJJ
# HHHIIIJJ"

design = "
AAAAABE
AAAAACE
AAAAADE
FFGHHII
FFGHHII"

# A should be 150 wide, and 75 tall (2x as wide as tall), B-D (and E-G) should be same wide as tall, H-I should be same wide as tall

# patchwork <- p3 + pc1_box + pc2_box + tr_box + plot_spacer() + plot_spacer() + plot_spacer() + pls.grid.plot.arsq + pc1_vs_pc2_plot + plot_spacer() + plot_layout(guides = "collect", design = design, width = unit(c(35), "mm"), height = unit(c(35), "mm"))

patchwork <- p3 + pc1_box + pc2_box + tr_box + plot_spacer() + pls.grid.plot.arsq + body_morph_plot + tr_sm + pc1_vs_pc2_plot + plot_layout(guides = "collect", design = design, width = unit(c(35), "mm"), height = unit(c(35), "mm"))

pdf("figure_2_mockup_rev.pdf", width = 20, height = 15)
patchwork + plot_annotation(tag_levels = 'A')
dev.off()









### Finding distances between points (eye size)

test <- two.d.array(body)

xenbat <- as.data.frame(body[,, dimnames(body)[[3]] %in% "Xenbat"])
colnames(xenbat) <- c("X", "Y")
xenbat$labels <- c(1:20)
## Plot to see if I got it right - I didn't, because #16 was removed, all of these numbers are shifted (I think there was)
ggplot(xenbat, aes(x = X, y = Y, label = labels)) + geom_point() + geom_text_repel()

lmks <- matrix(c(16,18,17,19,1,7), ncol=2, byrow=TRUE, dimnames = list(c("eyeW", "eyeH", "SL"),c("start", "end")))

lineardists <- as.data.frame(interlmkdist(body, lmks))

lineardists$ratio_width <- lineardists[,1] / lineardists[,3]
lineardists$ratio_height <- lineardists[,2] / lineardists[,3]
lineardists$ratio_eyes <- lineardists[,1] / lineardists[,2]

sleep_data$eye_ratio_w <- lineardists$ratio_width[match(sleep_data$six_letter_name_Ronco, row.names(lineardists))]
sleep_data$eye_ratio_h <- lineardists$ratio_height[match(sleep_data$six_letter_name_Ronco, row.names(lineardists))]
sleep_data$eye_ratio_eyes <- lineardists$ratio_eyes[match(sleep_data$six_letter_name_Ronco, row.names(lineardists))]


pc1_erw <- ggplot(sleep_data, aes(x = pc1, y = eye_ratio_w, colour = tribe)) + geom_point()
pc1_erh <- ggplot(sleep_data, aes(x = pc1, y = eye_ratio_h, colour = tribe)) + geom_point()
pc2_erw <- ggplot(sleep_data, aes(x = pc2, y = eye_ratio_w, colour = tribe)) + geom_point()
pc2_erh <- ggplot(sleep_data, aes(x = pc2, y = eye_ratio_h, colour = tribe)) + geom_point()
tr_erw <- ggplot(sleep_data, aes(x = total_rest, y = eye_ratio_w, colour = tribe)) + geom_point()
tr_erh <- ggplot(sleep_data, aes(x = total_rest, y = eye_ratio_h, colour = tribe)) + geom_point()

## Focus on width, which has the best relationship (but still not phylogenetically significant)
pc1_erw <- pc1_erw + xlab("PC1 (diurnal - nocturnal preference)") + ylab("Relative eye width") + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773"))
pc1_erw <- pc1_erw + annotate("text", x=-0.1, y=0.14, label= "Linear (PGLS) \n p-value = 0.4301", hjust = 0)

eye_size_category_plot <- ggplot(sleep_data, aes(x = category, y = eye_ratio_w, group = category)) + geom_boxplot(alpha = 0.25, colour = "black", fill = "grey", outlier.colour = "transparent") + geom_jitter(colour = "black", width = 0.15) + xlab("") + theme_classic()
eye_size_category_plot <- eye_size_category_plot + annotate("text", x = 2, y=0.14, label= "Phylogenetic ANOVA \n p-value = 0.09291", hjust = 0)

design2 <- "
AAAB"
eye_plot <- eye_size_category_plot + pc1_erw + plot_layout(guides = "collect", design = design2, width = unit(c(35), "mm"), height = unit(c(35), "mm"))

pdf("figure_eye_plot_rev.pdf", width = 20, height = 15)
eye_plot 
dev.off()





ggplot(sleep_data, aes(x = pc1, y = eye_ratio_w, colour = tribe, label = six_letter_name_Ronco)) + geom_point() + scale_colour_viridis() + geom_text_repel()
ggplot(sleep_data, aes(x = category, y = eye_ratio_w, colour = pc1, label = six_letter_name_Ronco)) + geom_boxplot() + geom_jitter(width = 0.1, size = 3) + scale_colour_viridis() + geom_text_repel()

test <- aov(eye_ratio_w ~ category, sleep_data)
summary.lm(test)
eye_ratio_w <- sleep_data$eye_ratio_w
names(eye_ratio_w) <- sleep_data$six_letter_name_Ronco
eye_ratio_w <- eye_ratio_w[!is.na(eye_ratio_w)]

sleep_data_category <- factor(sleep_data$category, levels = c("nocturnal", "diurnal", "crepuscular", "cathemeral"))
names(sleep_data_category) <- sleep_data$six_letter_name_Ronco
sleep_data_category <- sleep_data_category[names(eye_ratio_w)]

summary(aov.phylo(eye_ratio_w ~ sleep_data_category, cichlidTree, nsim = 1000, test = "Wilks"))

comp.data <- comparative.data(cichlidTree, sleep_data, six_letter_name_Ronco, vcv = TRUE)
mod <- pgls(eye_ratio_w ~ pc1, comp.data, lambda = "ML")
summary(mod)




#create a new variable for pc1_2
sleep_data$pc1_2 <- sleep_data$pc1^2
sleep_data$pc2_2 <- sleep_data$pc2^2
sleep_data$total_rest_2 <- sleep_data$total_rest^2

#fit quadratic regression model
quadraticModel <- lm(eye_ratio_w ~ pc1 + pc1_2, data=sleep_data)


# create sequence of hour values
Values <- seq(-0.2, 0.2, 0.001)

## Make a trend line for the quadratic fit
#create list of predicted happines levels using quadratic model
Predict <- predict(quadraticModel,list(pc1=Values, pc1_2=Values^2))

# PC2 = -4.626(pc1)^2 + -0.026(hours) - 0.191

compD = comparative.data(phy = cichlidTree, data=sleep_data, six_letter_name_Ronco)
fitmod <- pgls(eye_ratio_w ~ pc1 + pc1_2, compD)
summary_fitmod <- summary(fitmod)

# PC2 = -4.083814(pc1)^2 + 0.035057(hours) - 0.186

Predict_pgls <- predict(fitmod,list(pc1=Values, pc1_2=Values^2))

#create scatterplot of original data values
trend_line_data <- data.frame(pc1 = Values, eye_ratio_w = Predict_pgls)

#add predicted lines based on quadratic regression model
sleep_data$tribe <- factor(sleep_data$tribe, levels = c("Boulengerochromini", "Cyphotilapiini", "Cyprichromini", "Ectodini", "Eretmodini", "Lamprologini", "Limnochromini", "Trophenini"))
test_plot <- ggplot(sleep_data, aes(x = pc1, y = eye_ratio_w, colour = tribe, label = six_letter_name_Ronco)) + geom_point(size = 4) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("PC2 (crepuscular preference)") + xlab("PC1 (diurnal - nocturnal preference)") #+ geom_text_repel(force = 2) 
test_plot <- test_plot + geom_line(data = trend_line_data, aes(x=pc1, y = eye_ratio_w), inherit.aes = FALSE, size = 1.5, alpha = 1, colour = "black")

## Update #s based on summary(quandraticModel) and summary_fitmod)
pc1_vs_pc2_plot <- pca_plot + annotate("text", x=-0.1, y=0.01, label= "Quadratic \n r^2 = 0.6334 \n p-value = 3.787e-13 \n Quadratic (PGLS) \n r^2 = 0.5302 \n p-value = 2.023e-09", hjust = 0)





library(sp)

x1 <- rnorm(100, 0.8, 0.3)
y1 <- rnorm(100, 0.8, 0.3)
hpts <- chull(x = xenbat$X, y = xenbat$Y)
hpts <- c(hpts, hpts[1])
xy.coords <- cbind(x1, y1)
chull.coords <- xy.coords[hpts,]
chull.poly <- Polygon(chull.coords, hole=F)
chull.area <- chull.poly@area

test <- lapply(ecomorpho.tree$tip.label, function(species) {
  data <- as.data.frame(body[,, dimnames(body)[[3]] %in% species])
  colnames(data) <- c("X", "Y")
  data$labels <- c(1:20)
  x <- data[,1]
  y <- data[,2]
  hpts <- chull(x = x, y = y)
  hpts <- c(hpts, hpts[1])
  xy.coords <- cbind(x, y)
  chull.coords <- xy.coords[hpts,]
  chull.poly <- Polygon(chull.coords, hole=F)
  chull.area <- chull.poly@area
  
  # eyes
  
  hpts <- chull(x = x[16:19], y = y[16:19])
  hpts <- c(hpts, hpts[1])
  xy.coords <- cbind(x, y)
  chull.coords <- xy.coords[hpts,]
  chull.poly <- Polygon(chull.coords, hole=F)
  chull.area.eyes <- chull.poly@area
  
  return(chull.area.eyes/chull.area)
  
})

df <- data.frame(names = ecomorpho.tree$tip.label, ratio = unlist(test))
sleep_data$eye_area_ratio <- df$ratio[match(sleep_data$six_letter_name_Ronco, df$names)]

test <- aov(eye_area_ratio ~ category, sleep_data)
summary.lm(test)

eye_area_ratio <- sleep_data$eye_area_ratio
names(eye_area_ratio) <- sleep_data$six_letter_name_Ronco
eye_ratio_eyes <- sleep_data$eye_ratio_eyes
names(eye_ratio_eyes) <- sleep_data$six_letter_name_Ronco
eye_ratio_w <- sleep_data$eye_ratio_w
names(eye_ratio_w) <- sleep_data$six_letter_name_Ronco
eye_ratio_h <- sleep_data$eye_ratio_h
names(eye_ratio_h) <- sleep_data$six_letter_name_Ronco

aov.test <- aov.phylo(eye_ratio_w ~ sleep_data_category, cichlidTree, nsim = 1000, test = "Wilks")
summary(aov.test)

