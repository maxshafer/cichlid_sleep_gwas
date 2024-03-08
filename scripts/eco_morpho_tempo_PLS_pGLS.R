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

########################################################################################################
################ Make tree figures    ##################################################################
########################################################################################################

ggtree(cichlidTree, layout = "circular") + theme_tree(bgcolor = NA) + geom_tiplab() + geom_text(aes(label = node))

lamps <- 57
ectos <- 93
trops <- 102
cyps <- 99
erets <- 47
limnos <- 37
cyphos <- 36
bouls <- 1
#haplos <- 53


phylo.plot <- ggtree(as.phylo(cichlidTree), layout = "rectangular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 6)

sleep_data <- sleepData
sleep_data$x <- 11
sleep_data$y <- phylo.plot$data$y[match(sleep_data$species, phylo.plot$data$label)]
sleep_data$diet <- meta_data$diet[match(sleep_data$species, meta_data$six_letter_name_Ronco)]


# p3 <- phylo.plot + geom_tile(data = sleep_data, aes(y=y, x=x, fill = pc1), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", direction = 1) #+ scale_color_discrete(na.value = 'transparent')
# p3 <- p3 + new_scale("fill") + new_scale("size") + geom_tile(data = sleep_data, aes(y=y, x=x+2.25, fill = pc2), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "BrBG", direction = -1)
# p3 <- p3 + new_scale("size") + new_scale("fill") + geom_tile(data = sleep_data, aes(y=y, x=x+4.5, fill = total_rest), width = 2, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "PRGn", direction = -1) 
# p3

phylo.plot <- ggtree(as.phylo(cichlidTree), layout = "rectangular") + theme_tree(bgcolor = NA) + geom_tiplab(offset = 0)

p3 <- phylo.plot + geom_tile(data = sleep_data, aes(y=y, x=x+1.5, fill = pc1), width = 1, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "RdBu", direction = 1) #+ scale_color_discrete(na.value = 'transparent')
# p3 <- p3 + new_scale("colour") + new_scale("size") + geom_point(data = sleep_data, aes(y=y, x=x+0.75, colour = pc2), size = 6) + scale_colour_distiller(palette = "BrBG", direction = -1)
p3 <- p3 + new_scale("fill") + new_scale("size") + geom_tile(data = sleep_data, aes(y=y, x=x+2.75, fill = pc2), width = 1, height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "BrBG", direction = -1)
p3 <- p3 + new_scale("size") + new_scale("fill") + geom_tile(data = sleep_data, aes(y=y, x=x+4.5, fill = total_rest, width = total_rest/9), height = 0.9, size = 1.5, inherit.aes = FALSE) + scale_fill_distiller(palette = "PRGn", direction = -1) 



# ggplot(sleep_data, aes(x = pc1, y = total_rest, colour = total_rest, size = total_rest, label = species)) + geom_point() + scale_colour_viridis() + new_scale("size") + new_scale("colour") + geom_text_repel() + theme_classic() + ylab("Total Rest") + xlab("PC1 (diurnal - nocturnal preference)")
# ggplot(sleep_data, aes(y = pc2, x = total_rest, colour = total_rest, size = total_rest, label = species)) + geom_point() + scale_colour_viridis() + new_scale("size") + new_scale("colour") + geom_text_repel() + theme_classic() + xlab("Total Rest") + ylab("PC2 (crepuscular preference)")

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

#view model summary
summary(quadraticModel)

#create sequence of hour values
Values <- seq(-0.2, 0.2, 0.001)
#create list of predicted happines levels using quadratic model
Predict <- predict(quadraticModel,list(pc1=Values, pc1_2=Values^2))

# PC2 = -4.626(pc1)^2 + -0.026(hours) - 0.191

compD = comparative.data(phy = cichlidTree, data=sleep_data, species)
fitmod <- pgls(pc2 ~ pc1 + pc1_2, compD)

summary_fitmod <- summary(fitmod)

# PC2 = -4.083814(pc1)^2 + 0.035057(hours) - 0.186

Predict_pgls <- predict(fitmod,list(pc1=Values, pc1_2=Values^2))

#create scatterplot of original data values
trend_line_data <- data.frame(pc1 = Values, pc2 = Predict_pgls)
#add predicted lines based on quadratic regression model

# tribe	color
# bathy	"#242626"
# benth	"#AE262A"
# boule	"#59595C" ##
# cypho	"#FDDF13" ##
# cypri	"#F04D29" ##
# ectod	"#9AB9D9" ##
# eretm	"#682E7A" ##
# oreoc	"#845F25"
# lampr	"#C588BB" ##
# limno	"#535CA9" ##
# peris	"#FBAC43"
# tyloc	"#949598"
# trema	"#959170"
# troph	"#86C773" ##

sleep_data$tribe <- factor(sleep_data$tribe, levels = c("Boulengerochromini", "Cyphotilapiini", "Cyprichromini", "Ectodini", "Eretmodini", "Lamprologini", "Limnochromini", "Trophenini"))

pca_plot <- ggplot(sleep_data, aes(x = pc1, y = pc2, colour = tribe, label = species)) + geom_point(size = 4) + scale_colour_manual(values = c("#59595C", "#FDDF13", "#F04D29", "#9AB9D9", "#682E7A", "#C588BB", "#535CA9","#86C773")) + theme_classic() + ylab("PC2 (crepuscular preference)") + xlab("PC1 (diurnal - nocturnal preference)") #+ geom_text_repel(force = 2) 
pca_plot <- pca_plot + geom_line(data = trend_line_data, aes(x=pc1, y = pc2), inherit.aes = FALSE, size = 1.5, alpha = 1, colour = "black")

pc1_vs_pc2_plot <- pca_plot + annotate("text", x=-0.1, y=0.01, label= "Quadratic \n r^2 = 0.6027 \n p-value = 2.245e-11 \n Quadratic (PGLS) \n r^2 = 0.5145 \n p-value = 3.723e-09", hjust = 0)


p4 <- p3 + plot_layout(guides = "collect", width = unit(c(100), "mm"), height = unit(c(100), "mm"))


# tr_rain <- ggplot(sleep_data, aes(x = diet, y = total_rest, fill = diet)) + geom_rain(alpha = 0.5, rain.side = "r") + xlab("") + theme_classic() + coord_flip()
# pc1_rain <- ggplot(sleep_data, aes(x = diet, y = pc1, fill = diet)) + geom_rain(alpha = 0.5, rain.side = "r") + xlab("") + theme_classic() + coord_flip()
# pc2_rain <- ggplot(sleep_data, aes(x = diet, y = pc2, fill = diet)) + geom_rain(alpha = 0.5, rain.side = "r") + xlab("") + theme_classic() + coord_flip()

diet.colours <- c('Invertivore' = "#EF4026", 'Piscivore' = "#4682b4", 'Zooplanktivore' = "#f4a460", 'Algivore' = "#3cb371") 

# 'Invertivore': 'tomato', 'Piscivore': 'steelblue', 'Zooplanktivore': 'sandybrown',
# 'Algivore': 'mediumseagreen'

sleep_data$diet <- factor(sleep_data$diet, levels = c("Invertivore", "Piscivore", "Zooplanktivore", "Algivore"))

tr_box <- ggplot(sleep_data, aes(x = diet, y = total_rest, fill = diet)) + geom_boxplot(alpha = 0.25, colour = "black") + geom_jitter(aes(colour = diet), width = 0.15) + xlab("") + theme_classic() + coord_flip() + scale_colour_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371")) + scale_fill_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371"))
pc1_box <- ggplot(sleep_data, aes(x = diet, y = pc1, fill = diet)) + geom_boxplot(alpha = 0.25, colour = "black") + geom_jitter(aes(colour = diet), width = 0.15) + xlab("") + theme_classic() + coord_flip() + scale_colour_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371")) + scale_fill_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371"))
pc2_box <- ggplot(sleep_data, aes(x = diet, y = pc2, fill = diet)) + geom_boxplot(alpha = 0.25, colour = "black") + geom_jitter(aes(colour = diet), width = 0.15) + xlab("") + theme_classic() + coord_flip() + scale_colour_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371")) + scale_fill_manual(values = c("#EF4026", "#4682b4", "#f4a460", "#3cb371"))

# tr_box + pc1_box + pc2_box + plot_layout(ncol = 3, guides = "collect")

### Below is for doing ANOVAs, but there is no phylogenetically significant result

# test <- aov(d13C ~ diet, sleep_data)
# summary.lm(test)
# 
# sleep_data_diet <- sleep_data$diet[match(cichlidTree$tip.label, sleep_data$species)]
# names(sleep_data_diet) <- cichlidTree$tip.label
# sleep_data_pc1 <- sleep_data$y[match(cichlidTree$tip.label, sleep_data$species)]
# names(sleep_data_pc1) <- cichlidTree$tip.label
# 
# test <- phylANOVA(cichlidTree, sleep_data_diet, sleep_data_pc1, nsim=1000, posthoc=TRUE, p.adj="holm")

########################################################################################################
################ compute phylogenetic signal    ########################################################
########################################################################################################

## Pagel’s lambda is a measure of phylogenetic ‘signal’ in which the degree to which shared history of taxa 
## has driven trait distributions at tips. In this model, internal branch lengths are transformed by the lambda 
## parameter value. When the parameter lambda equals 1, branches are transformed by multiplying by 1 and so 
## the model is equal to Brownian motion (high phylogenetic signal). Values of lambda under 1 suggest there has 
## been less influence of shared history on trait values at the tips. Finally, a lambda value of 0 indicates no 
## phylogenetic influence on trait distributions, and is equivalent to a ‘star phylogeny’ with no shared branch lengths.

trait.data.pc1 <- sleepData$pc1[match(cichlidTree$tip.label, sleepData$species)]
names(trait.data.pc1) <- cichlidTree$tip.label

pc1_l <-phylosig(cichlidTree, trait.data.pc1, method = "lambda", test = T)
pc1_K <- phylosig(cichlidTree, trait.data.pc1, method = "K", test = T)

trait.data.pc2 <- sleepData$pc2[match(cichlidTree$tip.label, sleepData$species)]
names(trait.data.pc2) <- cichlidTree$tip.label

pc2_l <- phylosig(cichlidTree, trait.data.pc2, method = "lambda", test = T)
pc2_K <- phylosig(cichlidTree, trait.data.pc2, method = "K", test = T)

trait.data.tr <- sleepData$total_rest[match(cichlidTree$tip.label, sleepData$species)]
names(trait.data.tr) <- cichlidTree$tip.label

tr_l <- phylosig(cichlidTree, trait.data.tr, method = "lambda", test = T)
tr_K <- phylosig(cichlidTree, trait.data.tr, method = "K", test = T)




########################################################################################################
################ OK run pairwise pgls    ###############################################################
########################################################################################################

## Maybe I can use combn()?
## sleepData has all data for this

variables <- c("pc1", "pc2", "total_rest", "body_PC1", "body_PC2", "LPJ_PC1", "LPJ_PC2", "oral_PC1", "oral_PC2", "d15N", "d13C")

comp.data <- comparative.data(cichlidTree, sleepData, species, vcv = TRUE)


# ## OK, this doesn't give me all possible combinations, which doesn't make for an interesting grid
# output <- combn(variables, 2, function(x) {
# 
# mod <- pgls(as.formula(paste(x[1], "~", x[2])), comp.data, lambda = "ML")
# 
# return(mod) })
# 
# ## Need to add back names to output, using a call to pgls
# mod <- pgls(as.formula(paste(variables[1], "~", variables[2])), comp.data, lambda = "ML")
# row.names(output) <- names(mod)
# 
# # this now works
# mod <- output[,1] %>% "class<-"("pgls")
# summary(mod)$r.squared


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

grid.plot.pval <- ggplot(out, aes(x = namex, y = namey, fill = log(p.value)*-1)) + geom_tile() + scale_fill_viridis(direction = 1)
grid.plot.rsq <- ggplot(out, aes(x = namex, y = namey, fill = r.squared)) + geom_tile() + scale_fill_viridis(direction = 1)




########################################################################################################
################ OK run pairwise PLS / pGLS    #########################################################
########################################################################################################


###### OK, so I ran the scripts in ronco_et_al to generate the input data for her PLS analysis, and can now use it against our data!
## Try the geomorph two-block pls (it's somehow like cca?)
## This works, and is not significant, again I worry about the non-correlated nature of the paired PC scores...

## Load and prep body morpho

load("~/Documents/R_Projects/ronco_et_al/trait_evolution/02_Morpho_Eco/body/sp_means_body.Rdata")
body= sp_means_body[,, dimnames(sp_means_body)[[3]] %in% sleepData$species ]
body= body[,, dimnames(body)[[3]] %in% cichlidTree$tip.label]

## Load and prep UOJ morpho
## This is a subset of the whole body morpho data above

load("~/Documents/R_Projects/ronco_et_al/trait_evolution/02_Morpho_Eco/oral/sp_means_OJ.Rdata")
uoj= sp_means_OJ[,, dimnames(sp_means_OJ)[[3]] %in% sleepData$species ]
uoj= uoj[,, dimnames(uoj)[[3]] %in% cichlidTree$tip.label]

## Load and prep LPJ morpho

load("~/Documents/R_Projects/ronco_et_al/trait_evolution/02_Morpho_Eco/LPJ/sp_means_LPJ.Rdata")
lpj= sp_means_LPJ[,, dimnames(sp_means_LPJ)[[3]] %in% sleepData$species ]
lpj= lpj[,, dimnames(lpj)[[3]] %in% cichlidTree$tip.label]


# Prep isotope data
isotopes <- sleepData[,c("d15N", "d13C")]

# Prep speed and total rest data
speed <- as.data.frame(speed_data[row.names(isotopes),])
rest <- sleepData[,c("total_rest"), drop = FALSE]

# coefs <- sleepData[,c("intra_coef", "inter_coef"), drop = FALSE]





## Run two block pls for speed and rest vs isotopes, body, uoj, and lpj

variables.pls <- c("speed","rest", "body", "uoj", "lpj", "isotopes")

# comp.data <- comparative.data(cichlidTree, sleepData, species, vcv = TRUE)

## This gives me the full grid, and also keeps them as the correct class objects
grid.pls <- expand.grid(variables.pls, variables.pls)

output.grid.pls <- apply(grid.pls, 1, function(x) {
  
  if (x[1] == x[2]) {
    return ("NA")
  } else {

    f <- two.b.pls(eval(as.name(paste(x[1]))), eval(as.name(paste(x[2]))))
    f1 <- data.frame( "XScores"=f$ XScores[,1], "YScores"=f$ YScores[,1] , "species"= sleepData$species)
    
    compset = comparative.data(cichlidTree, f1, species)
    # Response ~ Predictor
    mod <- pgls(YScores ~ XScores, compset, lambda= 'ML')
    return(mod)
  }
  
})

output.grid.pls <- output.grid.pls[!(is.na(output.grid.pls))]



out.pls <- data.frame(namex = grid.pls$Var1,
                  namey = grid.pls$Var2,
                  r.squared = unlist(lapply(output.grid.pls, function(x) ifelse(x[1] == "NA", NA, summary(x)$r.squared))),
                  adj.r.squared = unlist(lapply(output.grid.pls, function(x) ifelse(x[1] == "NA", NA, summary(x)$adj.r.squared))),
                  p.value = unlist(lapply(output.grid.pls, function(x) ifelse(x[1] == "NA", NA, lmp(x)))))

out.pls$label <- ifelse(out.pls$p.value < 0.0001, "****", ifelse(out.pls$p.value < 0.001, "***", ifelse(out.pls$p.value < 0.01, "**", ifelse(out.pls$p.value < 0.05, "*", "")) ) )

# Maybe remove body vs uoj, since they have overlapping data (16 and 21)

out.pls[out.pls$namex == "uoj" & out.pls$namey == "body",c(3:5)] <- NA
out.pls[out.pls$namex == "body" & out.pls$namey == "uoj",c(3:5)] <- NA


## Plot

# Y ~ X, meaning, does X predict Y
pls.grid.plot.pval <- ggplot(out.pls, aes(x = namex, y = namey, fill = log(p.value)*-1)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = label), colour = "white")
pls.grid.plot.rsq <- ggplot(out.pls, aes(x = namex, y = namey, fill = r.squared)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = round(r.squared, 2)), colour = "white")
# pls.grid.plot.arsq <- ggplot(out.pls, aes(x = namex, y = namey, fill = adj.r.squared)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = round(adj.r.squared, 2)), colour = "white")
pls.grid.plot.arsq <- ggplot(out.pls, aes(x = namex, y = namey, fill = adj.r.squared)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = paste(round(adj.r.squared, 2), "\n", label, sep = "")), colour = "white")


##### Redo, to make plots


output.grid.data <- apply(grid.pls, 1, function(x) {
  
  if (x[1] == x[2]) {
    return (NA)
  } else {
    
    f <- two.b.pls(eval(as.name(paste(x[1]))), eval(as.name(paste(x[2]))))
    f1 <- data.frame( "XScores"=f$ XScores[,1], "YScores"=f$ YScores[,1] , "species"= sleepData$species)
    
    f1$tribe <- sleepData$tribe
    f1$d13C <- sleepData$d13C
    f1$d15N <- sleepData$d15N
    f1$pc1 <- sleepData$pc1
    f1$tr <- sleepData$total_rest
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
    return(plot)
  }

})


# wrap_plots(tribe.plots) + plot_layout(guides = "collect")



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


pdf("pairwise_eco_morpho_PLS_pGLS_plots_pvalue_rsquared_grid.pdf", width = 14, height = 10)
grid.plot.pval + grid.plot.rsq + pls.grid.plot.pval + pls.grid.plot.rsq + plot_layout(ncol = 2)
dev.off()

pdf("pairwise_eco_morpho_PLS_pGLS_plots_xvsy.pdf", width = 15, height = 15)
wrap_plots(tribe.plots) + plot_layout(guides = "collect", width = unit(c(15), "mm"), height = unit(c(15), "mm"))
dev.off()




corr_data <- read.csv("cichild_corr_data.csv")
pheno <- read.csv("pheno_data/combined_cichlid_data_2022-11-17_new_peaks.csv")



design = "
AAAABFF
AAAABFF
AAAACFF
AAAACGG
AAAADGG
AAAADGG
EEEEEEE"


patchwork <- p3 + pc1_box + pc2_box + tr_box + plot_spacer() + pls.grid.plot.rsq + pc1_vs_pc2_plot + plot_layout(guides = "collect", design = design, width = unit(c(50), "mm"), height = unit(c(25), "mm"))

pdf("figure_3_mockup.pdf", width = 20, height = 15)
patchwork + plot_annotation(tag_levels = 'A')
dev.off()




