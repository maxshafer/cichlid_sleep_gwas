library(GenomicFeatures)
library(ggrepel)
library(data.table)
library(dplyr)
library(Rgb)
library(readr)
library(patchwork)
library(Gviz)
library(biomaRt)

# list = ls()
# rm(list = list[!(list %in% c("allele_freq", "gtf", "gtf2", "gtf3", "gtf4"))])

# This should be run in the folder with the allele_frequncies
setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/sra_reads_nobackup/alleleFreqs")

## Make a script to read in all of the allelefreq files, combine them, plot the relationships, then subset by p-values
## This can then be modified, and re-run if we want other plots produced!

source('/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/scripts/variants_functions.R')

# Read in gtf file for finding genes
gtf <- read.gtf("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff")
gtf2 <- gtf[gtf$feature == "gene",]
gtf3 <- gtf[gtf$feature == "mRNA",]
gtf4 <- gtf[gtf$feature == "CDS",]
# rm(gtf)

# Setting arguments for testing script
# Need to make sure to specify the same phenotype file as was used in the gwas (if testing fewer species)
args <- c("../../pheno_data/combined_cichlid_data_2022-05-05.csv", 1, "day_night_dif")

Phenofile <- args[1]
cutoff <- as.numeric(args[2])
pheno_to_test <- as.character(args[3])

################################################################################################################################################################################################
#### Load the data  ############################################################################################################################################################################
################################################################################################################################################################################################

# Retrieve dataset names
# Should have files per chromosome, per trait, and per test (pi or PGLSpi)
if (pheno_to_test == "peak") {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
  datasets <- datasets[!(grepl("amplitude", datasets))]
  datasets <- datasets[grep(".txt", datasets)]
} 
if (pheno_to_test == "day_night_dif") {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
  datasets <- datasets[!(grepl("-species", datasets))]
  datasets <- datasets[grep(".txt", datasets)]
}  else {
  datasets <- list.files()[grep(pheno_to_test, list.files())]
  datasets <- datasets[grep(".txt", datasets)]
}

datasets.gwas <- datasets[grep("_piVals", datasets)]
datasets.pgls <- datasets[grep("_PGLSpiVals", datasets)]

gwas.names <- substr(datasets.gwas, start = 52, stop = 60)
pgls.names <- substr(datasets.pgls, start = 52, stop = 60)

# The above are all of the GWAS results per chromosome
# Read in and combine into one dataframe for plotting/subsetting

gwas.datasets <- lapply(datasets.gwas, function(x) fread(x, showProgress = F, drop = "V1"))
gwas.datasets <- lapply(gwas.datasets, function(x) x[x$ps < cutoff,]) # this keeps only the top ~10% of snps
names(gwas.datasets) <- gwas.names

pgls.datasets <- lapply(datasets.pgls, function(x) fread(x, showProgress = F, drop = "V1"))
pgls.datasets <- lapply(pgls.datasets, function(x) {
  colnames(x) <- c("chr", "coord", "PGLSps")
  return(x)
})
pgls.datasets <- lapply(pgls.datasets, function(x) x[x$PGLSps < cutoff,]) # this keeps only the top ~10% of snps
names(pgls.datasets) <- pgls.names

################################################################################################################################################################################################
#### Retrieve the info for gene lists  #########################################################################################################################################################
################################################################################################################################################################################################

genes_of_interest <- getLociEns(gene_list = list("kcnc2","LOC100698148", "LOC100701437"), window_padding = 0, gtf_file = gtf2, directory = "/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/sra_reads_nobackup/alleleFreqs/", gwas.datasets = gwas.datasets, pgls.datasets = pgls.datasets, pheno_to_test = pheno_to_test)


gene_locations <- gtf2[gtf2$gene %in% list("kcnc2","LOC100698148", "LOC100701437"),]
gene_locations <- gene_locations[order(gene_locations$seqname),] # I think this is not required (they are already in order)
chromosomes <- substr(gene_locations$seqname, start = 1, stop = 9)
gene_locations$chromosome <- chromosomes

plots <- lapply(genes_of_interest, function(x) {
  plot <- ggplot(x, aes(x = coord, y = (-1*log(ps))+(-1*log(PGLSps)))) + geom_point() + theme_classic() + ggtitle(names(x))
  return(plot)
})

plots <- lapply(genes_of_interest, function(x) {
  plot <- ggplot(x, aes(x = coord, y = (-1*log(PGLSps)))) + geom_point() + theme_classic() + ggtitle(names(x))
  return(plot)
})

plots[["LOC100698148"]]



################################################################################################################################################################################################
#### Plot something?  ##########################################################################################################################################################################
################################################################################################################################################################################################
options(ucscChromosomeNames=FALSE)
# Load mart objects for biomTrack and for extracting gene start and end positions
ensembl <- useMart("ensembl")
ensembl <- useDataset("oniloticus_gene_ensembl", mart = ensembl)

# Add FST data as datatrack
# This works now, must have been the same issue regarding the chromosome name (with the .1 added now)

gene_of_interest_2 <- genes_of_interest[[3]]

gtrack <- GenomeAxisTrack()
ranges1 <- GRanges(seqnames = paste(unique(gene_of_interest_2$chr), ".1", sep = ""), ranges = IRanges(start=gene_of_interest_2$coord, end = gene_of_interest_2$coord))
data1 <- DataTrack(range = ranges1, data = -1*log(gene_of_interest_2$ps), name = "GWAS p-values", cex = .5, col = "blue")
data2 <- DataTrack(range = ranges1, data = -1*log(gene_of_interest_2$PGLSps), name = "GWAS PGLS p-values", cex = .5, col = "red")
data3 <- DataTrack(range = ranges1, data = (-1*log(gene_of_interest_2$PGLSps))+(-1*log(gene_of_interest_2$PGLSps)), name = "GWAS PGLS p-values", cex = .5, col = "purple")

plotTracks(list(data3, gtrack))

# Hmm, this seems to work now, and uses the full data, then just plots part of it *shrug*

ranges2 <- GRanges(seqnames = paste(unique(gene_of_interest_2$chr), ".1", sep = ""), ranges = IRanges(start=gwas.datasets[[unique(gene_of_interest_2$chr)]]$coord, end = gwas.datasets[[unique(gene_of_interest_2$chr)]]$coord))
data1 <- DataTrack(range = ranges2, data = -1*log(gwas.datasets[[unique(gene_of_interest_2$chr)]]$ps), name = "GWAS p-values", cex = .5, col = "blue")
data2 <- DataTrack(range = ranges2, data = -1*log(pgls.datasets[[unique(gene_of_interest_2$chr)]]$PGLSps), name = "GWAS PGLS p-values", cex = .5, col = "red")
data3 <- DataTrack(range = ranges2, data = (-1*log(gwas.datasets[[unique(gene_of_interest_2$chr)]]$ps))+(-1*log(pgls.datasets[[unique(gene_of_interest_2$chr)]]$PGLSps)), name = "GWAS PGLS p-values", cex = .5, col = "purple")
# plotTracks(list(OverlayTrack(trackList = list(data1, data2), name = "OverlayTrack"), gtrack))


# Setup tracks
# This will never work, because ensembl doesn't have the version I use from Refseq
# biomTrack <- BiomartGeneRegionTrack(start = 1, end = 100000000, biomart = ensembl, chromosome = paste(unique(gene_of_interest_2$chr), ".1", sep = ""), name = paste(unique(gene_of_interest_2$chr), ".1", sep = ""), stacking = "squish")
#gtfTrack <- GeneRegionTrack("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff", chromosome = "NC_031972.1", start = min(gene_of_interest_2$coord), end = max(gene_of_interest_2$coord))
#strack <- SequenceTrack(biomart = ensembl, chromosome = unique(gene_of_interest_2$chr))

gtfTrack <- GeneRegionTrack("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff", stacking = "full")
chromosome(gtfTrack) <- paste(unique(gene_of_interest_2$chr), ".1", sep = "")
gtfTrack2 <- gtfTrack[feature(gtfTrack) %in% c("gene")]

shoulder <- 5000
plotTracks(c(data3, gtfTrack, gtrack), from = min(gene_of_interest_2$coord)-shoulder, to = max(gene_of_interest_2$coord)+shoulder, chromosome = paste(unique(gene_of_interest_2$chr), ".1", sep = ""), showID = TRUE, geneSymbol = TRUE, transcriptAnnotation = "feature", collapse = F)


txdbFromGFF <- makeTxDbFromGFF(file = "/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/O-niloticus_ASM185804v2/GCF_001858045.1_ASM185804v2_genomic.gff.gz") 
# customFromTxDb <- GeneRegionTrack(txdbFromGFF, chromosome = paste(unique(gene_of_interest_2$chr), ".1", sep = "")) 
customFromTxDb <- GeneRegionTrack(txdbFromGFF, chromosome = paste(unique(gene_of_interest_2$chr), ".1", sep = ""), from = min(gene_of_interest_2$coord), to = max(gene_of_interest_2$coord),) 

plotTracks(customFromTxDb, from = min(gene_of_interest_2$coord)-shoulder, to = max(gene_of_interest_2$coord)+shoulder, transcriptAnnotation="gene") 



# Plot tracks
shoulder <- 50000

plotTracks(c(gtrack, gtfTrack), type = c("p", "polygon"), fill.mountain = FALSE, transcriptAnnotation = "symbol", fontcolor.feature = "black", indel = "red", cex = 1, fontsize = 12, featureAnnotation = "id", just.feature = "right")

plotTracks(gtfTrack)



plotTracks(c(OverlayTrack(trackList = list(data1, data2, data3, data4), name = "OverlayTrack"), gtrack, biomTrack, indelTrack, strack), type = c("p", "polygon"), from = (gene.info[1,3]-shoulder), to = (gene.info[1,4]+shoulder), fill.mountain = FALSE, transcriptAnnotation = "symbol", fontcolor.feature = "black", indel = "red", cex = 1, fontsize = 12, featureAnnotation = "id", just.feature = "right")






