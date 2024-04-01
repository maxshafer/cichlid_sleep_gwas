## This script is for generating random subsets of SNPs for GO analysis
## Can be run with arguments for the number of permutations, and the number of SNPs to select per permutation
## Outputs an RDS object of length *number of permutations*

library(stringr)
library(ggplot2)
library(data.table)
library(dplyr)

setwd("/scicore/home/schiera/gizevo30/projects/cichlids_2")

########################################################################
#######   Defining trailing arguments  #################################
########################################################################

# Arguments would be the alleleFreq table(s) and the snpEff table, as well as the pvalue cutoff
# Can I just give it the chromosome id, then use that to pull all the tables?
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("must have 2 arguments for Rscript command")
}

# This is the number of SNPs to randomly pick
numb_pick <- as.numeric(args[1])
# This is the number of times to pick (how many random subsets)
times_pick <- as.numeric(args[2])

########################################################################
#######   Defining  arguments  #########################################
########################################################################

files <- list.files("sra_reads_nobackup/dump/combined_ann/", pattern = "_pvals_ann.gz")

########################################################################
#######   Load/merge/subset per comparison  ############################
########################################################################

dfs_all <- lapply(files, function(x) fread(paste("sra_reads_nobackup/dump/combined_ann/", x, sep = ""), showProgress = T))
names(dfs_all) <- files

dfs_all <- lapply(dfs_all, function(x) {
  x <- x[, c(2:3,16:21)]
  colnames(x) <- c("CHROM","POS","REF","ALT","ANN_GENE","ANN_IMPACT","ANN_EFFECT","ANN_DISTANCE")
  x$location <- paste(x$CHROM, x$POS, sep = ":")
  return(x)
})
  
merged <- Reduce(rbind, dfs_all)
snp_locations <- unique(merged$location)

output <- lapply(1:times_pick, function(i) {
  
  set.seed(i)
  sample <- sample(snp_locations, numb_pick)
  sample.out <- merged[merged$location %in% sample,]
  return(sample.out)
})

saveRDS(output, file = paste("sra_reads_nobackup/dump/combined_ann/Permute_SNPs_perchr", times_pick, "permutations", numb_pick, "SNPs", sep = "_"))


