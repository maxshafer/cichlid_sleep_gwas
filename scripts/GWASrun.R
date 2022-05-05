########################################################################
### Script by Milan Malinsky and Carolin Sommer-Trembo
### Last edit: 17th June 2021
### This is an executable R script for running the standard linear-model GWAS 
### i.e. associating species-level allele frequencies with the behaviour median values for those species
### Takes some time (possibly hours) to run 
### Example to execute this script on the Sanger Farm:
# bsub -G rdgroup -e getDP_e_%J -o getDP_o_%J -R'select[mem>40000] rusage[mem=40000]' -M40000 /software/R-4.0.3/bin/Rscript GWASrun.R medians_for_GWAS_withoutfilter_speciesCodesAndGenomeIDs_wOretan_TanganyikaGWAS_wOreTan_NC_031965_fromProbabilities_AF.txt.gz  

### Modified by MERS 27.04.2022

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

# Get p-values from the GWAS object 
lmp <- function (modelobject) {
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}


########################################################################
##  The script to be executed                 ##########################
########################################################################

tanTree <- read.nexus("/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/05_BEAST_RAxML.tre")

############
### Read and process the Allele Frequencies from the file supplied by the command line argument
############

args = commandArgs(trailingOnly=TRUE)
AFfile <- args[1]
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
Phenofile <- args[2]
y <- read.table(Phenofile,header=T,sep=",")
y <- y[match(colnames(AFscaled), y$species_six),]


##########
### Remove SNPs with missing data (i.e. those where the allele frequency could not be calculated)
#########
N <- apply(AFscaled,1, function(x) length(which(is.na(x == T))))        # This calculates the number of missing data per SNP
AFscaled.noMiss <- AFscaled[which(N == 0),];				# This subsets the AF data to remove SNPs with missing data 

##########
# !!!!! Run the GWAS: y = X_a*beta_a + epsilon !!!!! - TAKES TIME
##########
pheno <- as.character(args[3])
GWASvector.noMiss <- apply(AFscaled.noMiss, 1, function(x) lm(y[,pheno] ~ x))

##########
# Calculate the p-values from the LM objects and write all the p-values out into a text file
##########
naVector <- sapply(GWASvector.noMiss, function(x) is.na(x$coefficients[2])); 
GWASvector.noMiss <- GWASvector.noMiss[-which(naVector == TRUE)]; 
ps <- sapply(GWASvector.noMiss, lmp); ps.withLoc <- cbind(AFfull[names(ps),1:2],ps)
write.table(ps.withLoc,file=paste0(substr(AFfile, 1, nchar(AFfile)-7),"_",args[3],"_piVals.txt"),quote=F,sep="\t")

##########
# Write out a subset of the LM R objects - reducing the overall size by choosing the ones with -log10(p-val) > 8 (unlikely to need to look in detail at the ones with less significant p-vals)
GWASvector.noMiss <- GWASvector.noMiss[which(-log10(ps) > 8)]
save(GWASvector.noMiss,file=paste0(substr(AFfile, 1, nchar(AFfile)-7),"_",args[3],".gwas"))


