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
################ OK run pairwise pgls    ###############################################################
########################################################################################################

## Maybe I can use combn()?
## sleepData has all data for this

variables <- c("pc1", "pc2", "total_rest", "body_PC1", "body_PC2", "LPJ_PC1", "LPJ_PC2", "oral_PC1", "oral_PC2", "d15N", "d13C", "inter_coef", "intra_coef")

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

coefs <- sleepData[,c("intra_coef", "inter_coef"), drop = FALSE]





## Run two block pls for speed and rest vs isotopes, body, uoj, and lpj

variables.pls <- c("speed","rest", "coefs", "body", "uoj", "lpj", "isotopes")

comp.data <- comparative.data(cichlidTree, sleepData, species, vcv = TRUE)

## This gives me the full grid, and also keeps them as the correct class objects
grid.pls <- expand.grid(variables.pls, variables.pls)

output.grid.pls <- apply(grid.pls, 1, function(x) {
  
  if (x[1] == x[2]) {
    return ("NA")
  } else {
    
    f <- two.b.pls(eval(as.name(paste(x[1]))), eval(as.name(paste(x[2]))))
    f1 <- data.frame( "XScores"=f$ XScores[,1], "YScores"=f$ YScores[,1] , "species"= sleepData$species)
    
    compset = comparative.data(cichlidTree, f1, species)
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

out.pls[16,c(3:5)] <- NA
out.pls[21,c(3:5)] <- NA


## Plot

pls.grid.plot.pval <- ggplot(out.pls, aes(x = namex, y = namey, fill = log(p.value)*-1)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = label), colour = "white")
pls.grid.plot.rsq <- ggplot(out.pls, aes(x = namex, y = namey, fill = r.squared)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = round(r.squared, 2)), colour = "white")
pls.grid.plot.arsq <- ggplot(out.pls, aes(x = namex, y = namey, fill = adj.r.squared)) + geom_tile() + scale_fill_viridis(direction = 1) + geom_text(aes(label = round(adj.r.squared, 2)), colour = "white")







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
    plot <- ggplot(x, aes(x=XScores, y=YScores, colour = tribe)) + geom_point() + theme_classic() + ggtitle(label = x$comp[1])
    return(plot)
  }

})


wrap_plots(tribe.plots) + plot_layout(guides = "collect")



ggplotRegression <- function(fit, data) {
  
  plot <- ggplot(data, aes(x = XScores, y = YScores, colour = tribe)) + geom_point() + theme_classic() + stat_smooth(method = "glm", col = "red")
  plot + labs(title = paste("Lambda =", round(summary(fit)$param.CI$lambda$opt, 3),
                              "Adj R2 = ",round(summary(fit)$adj.r.squared, 2),
                              " P =",round(summary(fit)$coef[2,4], 7)))
         
}
  
tribe.plots <- lapply(seq_along(output.grid.data), function(x) {
  
  if (is.na(output.grid.data[[x]])) {
    plot <- ggplot() + theme_void()
  } else {
    plot <- ggplotRegression(fit = output.grid.pls[[x]], data = output.grid.data[[x]])
    return(plot)
  }
  
  
})


pdf("pairwise_eco_morpho_PLS_pGLS_plots_pvalue_rsquared_grid.pdf", width = 14, height = 10)
grid.plot.pval + grid.plot.rsq + pls.grid.plot.pval + pls.grid.plot.rsq + plot_layout(ncol = 2)
dev.off()

pdf("pairwise_eco_morpho_PLS_pGLS_plots_xvsy.pdf", width = 30, height = 30)
wrap_plots(tribe.plots) + plot_layout(guides = "collect")
dev.off()






