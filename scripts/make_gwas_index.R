#!/usr/bin/env Rscript

library(dplyr)

# This is a few lines to make an index for chromosome by test, for use in the GWAS and PGLS scripts

intervals <- read.csv("/scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.chromosomes", head = F)
intervals <- intervals[1:23,1]

# tests <- read.csv("/scicore/home/schiera/gizevo30/projects/cichlids_2/pheno_data/GWAS_data_PCs_TR_2024-01-28.csv", row.names = "X")
# tests <- colnames(tests)[2:ncol(tests)]

tests <- c("pc1", "pc2", "total_rest")

test <- lapply(intervals, function(x) data.frame(test = tests, interval = x))
test2 <- Reduce(rbind, test)

write.table(test2, file = "/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/interval_test_array_chr.csv", quote = FALSE, col.names = FALSE, sep = ",")
