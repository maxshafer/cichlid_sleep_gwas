library(data.table)
library(stringr)

setwd("~/Documents/R_Projects/cichlid_sleep_gwas/")

########################################################################
#######   Defining trailing arguments  #################################
########################################################################

# Arguments would be the alleleFreq table(s) and the snpEff table, as well as the pvalue cutoff
# Can I just give it the chromosome id, then use that to pull all the tables?
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("must have 3 arguments for Rscript command")
}

# # Setting arguments for testing script
args <- c("NC_031965", "0.005")

chromosome <- args[1]
cutoff <- as.numeric(args[2])


files <- list.files("sra_reads_nobackup/alleleFreqs/", pattern = chromosome)
# AF_file <- files[grep("AF.txt", files)]
snpeff_file <- list.files("sra_reads_nobackup/final_snps_chr/", pattern = chromosome)

gwas_files <- files[grepl("_piVals_", files) & grepl(".txt", files)]
pgls_files <- files[grepl("_PGLSpiVals_", files) & grepl(".txt", files)]

gwas_names <- str_sub(gwas_files, start = 80, end = -5)
pgls_names <- str_sub(pgls_files, start = 80, end = -5)


########################################################################
#######   Load in files  ###############################################
########################################################################

# AF <- fread(paste("sra_reads_nobackup/alleleFreqs/", AF_file, sep = ""))

snpeff <- fread(paste("sra_reads_nobackup/final_snps_chr/", snpeff_file, sep = ""))

gwas <- lapply(gwas_files, function(x) fread(paste("sra_reads_nobackup/alleleFreqs/", x, sep = "")))

pgls <- lapply(pgls_files, function(x) fread(paste("sra_reads_nobackup/alleleFreqs/", x, sep = "")))


########################################################################
#######   Wrangle + Merge  #############################################
########################################################################

# Removes the ".1 from each chromosome identifier (though probably isn't necessary)
if (paste(chromosome, ".1", sep = "") %in% snpeff$CHROM[1]) {
  snpeff$CHROM <- str_remove(snpeff$CHROM, pattern = "[.]1")
}

pgls <- lapply(seq_along(pgls), function(x) {
  colnames(pgls[[x]]) <- c("V1", "CHROM", "POS", pgls_names[x])
  return(pgls[[x]][,2:4])
})

gwas <- lapply(seq_along(gwas), function(x) {
  colnames(gwas[[x]]) <- c("V1", "CHROM", "POS", gwas_names[x])
  return(gwas[[x]][,2:4])
})

# merge them together, should use the CHROM and POS columns
gwas <- Reduce(merge, gwas)
pgls <- Reduce(merge, pgls)
combined <- merge(gwas, pgls)
rm(gwas, pgls)
combined <- merge(combined, snpeff)

# save out
gz1 <- gzfile(paste("sra_reads_nobackup/combined_ann/combined_", chromosome, "snps_pvals_ann.gz", sep = ""), "w")
write.csv(combined, gz1)
close(gz1)




