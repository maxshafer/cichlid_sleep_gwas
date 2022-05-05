#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


library(dplyr)

if (length(args)!=1) {
  stop("Only  one argument must be supplied, the intervals to use", call.=FALSE)
}

# args[1] should be the interval # to use, can be "40x", "80x", "1000x", "chr" 

# This is a few lines to make the index_array for either 40x or 1000x intervals and whichever samples I'm using (just 5 or all)

intervals <- list.files(paste("~/projects/cichlids_2/genome/intervals_ASM185804v2_", args[1], "/", sep = ""))

# samples <- list.files("~/ncbi/public/sra/")

samples <- read.csv("~/projects/cichlids_2/scripts/index_samples.csv", head = F, row.names = "V1")
samples <- samples$V2
test <- lapply(samples, function(x) data.frame(sample = x, interval = intervals))
test <- lapply(samples, function(x) data.frame(sample = x, interval = intervals))
test2 <- Reduce(rbind, test)

write.table(test2, paste("~/projects/cichlids_2/scripts/index_array_", args[1], ".csv", sep = ""), quote = FALSE, col.names = FALSE, sep = ",")

