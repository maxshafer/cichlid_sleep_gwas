library(data.table)
library(stringr)

setwd("/home/ayasha/scratch/temp_data/gwas/")

########################################################################
#######   Defining trailing arguments  #################################
########################################################################

# Arguments would be the alleleFreq table(s) and the snpEff table, as well as the pvalue cutoff
# Can I just give it the chromosome id, then use that to pull all the tables?
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("must have 1 arguments for Rscript command")
}

# # Setting arguments for testing script
# args <- c("NC_013663", "0.005")

chromosome <- args[1]
# cutoff <- as.numeric(args[2])
print(chromosome)

# for LG1 and 2, fails bc it gathers other files for LG10-19, etc
if (chromosome == "LG1" | chromosome == "LG2") {
   use_adjustment <- TRUE
   chromosome <- paste(chromosome, "_", sep = "")
} else { use_adjustment <- FALSE }

files <- list.files("alleleFreqs/", pattern = chromosome)
# AF_file <- files[grep("AF.txt", files)]

if (use_adjustment == TRUE) { chromosome <- str_replace(chromosome, "_", ".txt") }

snpeff_file <- list.files("final_snps_chr/", pattern = chromosome)
snpeff_file <- snpeff_file[grepl("extracted", snpeff_file)]

gwas_files <- files[grepl("_piVals_", files) & grepl(".txt", files)]
pgls_files <- files[grepl("_PGLSpiVals_", files) & grepl(".txt", files)]

gwas_names <- str_sub(gwas_files, start = 80, end = -5)
pgls_names <- str_sub(pgls_files, start = 80, end = -5)


########################################################################
#######   Load in files  ###############################################
########################################################################

# AF <- fread(paste("sra_reads_nobackup/alleleFreqs/", AF_file, sep = ""))

snpeff <- fread(paste("final_snps_chr/", snpeff_file, sep = ""))

gwas <- lapply(gwas_files, function(x) fread(paste("alleleFreqs/", x, sep = "")))

pgls <- lapply(pgls_files, function(x) fread(paste("alleleFreqs/", x, sep = "")))


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

########################################################################
#######   Summary statistics  ##########################################
########################################################################

## I will need to calculate the mean, min/max, and products per chromosome, then retain just these? No, I will probably
## want to look at the values in each, so I should save these out, then repeat a filtering step

comparison <- c("pref_sfi", "pref_max_cons", "dn_pref")
#comparison <- c("pref_sfi", "pref_max_cons")
columns <- unlist(lapply(comparison, function(x) grep(x, colnames(combined))))

combined$summary_mean <- rowMeans(combined[,..columns])
combined$summary_min <- apply(combined[,..columns], 1, function(x) min(x))
combined$summary_prod <- apply(combined[,..columns], 1, function(x) prod(x))

# save out
gz1 <- gzfile(paste("combined_ann/combined_", chromosome, "snps_pvals_ann.gz", sep = ""), "w")
write.csv(combined, gz1)
close(gz1)




