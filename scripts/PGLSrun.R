########################################################################
### Script by Milan Malinsky and Carolin Sommer-Trembo
### Last edit: 17th June 2021
### This is an executable R script for running GWAS with correction for phylogeny (but not gene-flow) using the PGLS method
### GWAS means associating species-level allele frequencies with the behaviour median values for those species
### Takes some time (definitely many hours) to run for each chromosome
### Example to execute this script on the Sanger Farm:
# bsub -G rdgroup -e getDP_e_%J -o getDP_o_%J -q long -R'select[mem>60000] rusage[mem=60000]' -M60000 /software/R-4.0.3/bin/Rscript PGLSrun.R medians_for_GWAS_withoutfilter_speciesCodesAndGenomeIDs_wOretan_TanganyikaGWAS_wOreTan_NC_031965_fromProbabilities_AF.txt.gz 

########################################################################
#######   Load packages  ######################################
########################################################################

library(caper)


########################################################################
#######   Defining trailing arguments  #################################
########################################################################

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
	stop("must have 3 arguments for Rscript command")
}

AFfile <- args[1]
Phenofile <- args[2]
pheno_to_test <- as.character(args[3])

########################################################################
#######   Defining functions that are used below  ######################
########################################################################

# Obviously - runs the PGLS
doPGSL <- function(x) {
	ourDF <- data.frame(pheno= y[,pheno_to_test],geno=as.numeric(x))
	rownames(ourDF) <- y$species_six
	ourDF$species_six = row.names(ourDF)
	compD = comparative.data(phy = tanTree2, data=ourDF, species_six)
	if (length(unique(x)) == 1) { print(x); return(1.0); }
	else {
		fitmod= pgls( pheno ~ geno, compD)
		f <- summary(fitmod)$fstatistic
		p <- pf(f[1],f[2],f[3],lower.tail=F)
		attributes(p) <- NULL
		return(p)
	}
}

########################################################################
##  The script to be executed                 ##########################
########################################################################

tanTree <- read.nexus("/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/05_BEAST_RAxML.tre")

############
### Read and process the Allele Frequencies from the file supplied by the command line argument
############

AF <- read.table(AFfile,header=T)
# -1 stands originally for missing data; then subset only the Allele Frequency columns
AF[AF==-1] <- NA; AFfull <- AF; AF <- AF[,5:ncol(AF)]   # ncol(AF) is the final column of the AF object - needs to change if using a different number of species

# Remove those that aren't in the tree file
AF <- AF[,colnames(AF) %in% tanTree$tip.label]
# Scale the actual Allele Frequncies to be on the interval [-1,1]  
AFscaled <- 2*((AF - min(AF,na.rm=T))/(max(AF,na.rm=T)-min(AF,na.rm=T))) - 1

############
### Read phenotype vector  
############ 
y <- read.table(Phenofile,header=T,sep=",")
y <- y[match(colnames(AFscaled), y$species_six),]

##########
### Remove SNPs with missing data (i.e. those where the allele frequency could not be calculated)
#########
N <- apply(AFscaled,1, function(x) length(which(is.na(x == T))))        # This calculates the number of missing data per SNP
AFscaled.noMiss <- AFscaled[which(N == 0),];                            # This subsets the AF data to remove SNPs with missing data 

##########
# The way things are done here, we actually still need to run the standard LM GWAS: 
# y = X_a*beta_a + epsilon 
# in order to find some SNPs for which for some reason related to missingness we can't get the slope of the lm
# and then we exclude these SNPs from the PGLS runs below
##########
GWASvector.noMiss <- apply(AFscaled.noMiss, 1, function(x) lm(y[,pheno_to_test] ~ x))
naVector <- sapply(GWASvector.noMiss, function(x) is.na(x$coefficients[2]));    ### this object is used below in the PGLS to exclude the SNPs 

##########
### Here we do the actual PGLS!!!!!!! takes many hours
### We load the phenotype vector again and process it slightly differently for this purpose - this is not optimal, but it works so let's not mess with it for no good reason
#########

tanTree2 <- drop.tip(tanTree, tanTree$tip.label[! tanTree$tip.label %in% y$species_six]);
pglsGWASvector <- apply(AFscaled.noMiss[-which(naVector == TRUE),], 1, doPGSL)
pglsGWASvector.withLoc <- cbind(AFfull[names(pglsGWASvector),1:2],pglsGWASvector)
write.table(pglsGWASvector.withLoc,file=paste0(substr(AFfile, 1, nchar(AFfile)-7),"_",args[3],"_PGLSpiVals_corrected.txt"),quote=F,sep="\t")

