


## Make a script to read in all of the allelefreq files, combine them, plot the relationships, then subset by p-values
## This can then be modified, and re-run if we want other plots produced!

# Load arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("must have 3 arguments for Rscript command")
}

# # Setting arguments for testing script
# args <- c("../cichlid_sleep_gwas/pheno_data/combined_cichlid_data_2022-05-05.csv", 0.001, "peak_amplitude")


Phenofile <- args[1]
cutoff <- as.numeric(args[2])
pheno_to_test <- as.character(args[3])

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/alleleFreqs")

# Retrieve dataset names
# Should have files per chromosome, per trait, and per test (pi or PGLSpi)
datasets <- list.files()[grep(pheno_to_test, list.files())]

datasets.gwas <- datasets[grep("_piVals_", datasets)]
datasets.pgls <- datasets[grep("_PGLSpiVals_", datasets)]

# The above are all of the GWAS results per chromosome
# Read in and combine into one dataframe for plotting/subsetting

gwas.perchr <- lapply(datasets.gwas, function(x) read.table(x, header = T))
gwas.perchr <- Reduce(rbind, gwas.perchr) # maybe remove the list after this, since they are each ~3Gb

pgls.perchr <- lapply(datasets.pgls, function(x) read.table(x, header = T))
pgls.perchr <- Reduce(rbind, pgls.perchr)

snps <- intersect(which(gwas.perchr[,3] < 0.05), which(pgls.perchr[,3] < 0.05))

gwas.perchr <- gwas.perchr[snps,]
pgls.perchr <- pgls.perchr[snps,]

perchr <- merge(gwas.perchr[snps,], pgls.perchr[snps,])


plot(log(perchr$ps), log(perchr$pglsGWASvector))

ggplot(perchr, aes(x = log(ps), y = log(pglsGWASvector))) + geom_point() + theme_classic() + scale_x_reverse() + scale_y_reverse()


pheno_data <- read.table(Phenofile,header=T,sep=",", row.names = "X")
