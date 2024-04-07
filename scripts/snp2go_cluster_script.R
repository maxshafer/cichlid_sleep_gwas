library(SNP2GO)
library(Rgb)
library(here)
library(stringr)

# ########################################################################
# #######   Defining GO file  ############################################
# ########################################################################
# 
# ## Create a GO file by replacing human gene names with orthologous oreochromis gene names
# ## Only needs to be done once
# 
# ## Load orthologs
# ortholog_df <- readRDS("GTF_human_orthologs_all.rds")
# 
# human_go_terms <- read.csv("mart_export_human_go_terms.txt", sep = "\t")
# human_go_terms <- human_go_terms[human_go_terms$Gene.name %in% ortholog_df$human_ortholog,]
# human_go_terms$Gene.name <- ortholog_df$oreochromis_gene[match(human_go_terms$Gene.name, ortholog_df$human_ortholog)]
# human_go_terms <- human_go_terms[,c(1,3)]
# write.table(human_go_terms, file = "mart_export_human_go_terms_oreochromis_orthos.txt", quote = F, row.names = F, sep = "\t")

## For a vcf file, the first colummn is the CHROM, the 2nd is the position, which I have in not a vcf file
## This should be easier than loading in some gigantic vcf file

########################################################################
#######   Load/merge per comparison  ###################################
########################################################################

files <- list.files("sra_reads_nobackup/dump/combined_ann/", pattern = "_pvals_ann.gz")

dfs_all <- lapply(files, function(x) fread(paste("sra_reads_nobackup/dump/combined_ann/", x, sep = ""), showProgress = T))
names(dfs_all) <- files

dfs_all <- lapply(dfs_all, function(x) {
  x <- x[, c(2:3)]
  colnames(x) <- c("CHROM","POS")
  x$location <- paste(x$CHROM, x$POS, sep = ":")
  return(x)
})

# This is the background
merged <- Reduce(rbind, dfs_all)

########################################################################
#######   Load significant SNP lists  ##################################
########################################################################
filter_snps <- readRDS(here("sra_reads_nobackup/dump/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

snps <- lapply(filter_snps, function(x) GRanges(seqnames=x[,1],ranges=IRanges(x[,2],x[,2])))

# Make GRanges object from all SNPs               
bkgrd_snps <- GRanges(seqnames=merged[,1],ranges=IRanges(merged[,2],merged[,2]))

########################################################################
#######   Run the analysis per candidate SNP list  E####################
########################################################################

out <- lapply(seq_along(snps), function(i) {
  
  # Candidate snps
  cand <- snps[[i]]
  # non-candidate snps
  noncand <- bkgrd_snps[!(merged$location %in% filter_snps[[i]]$location)]
  
  y <- snp2go(gtf="genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf",
              goFile="mart_export_human_go_terms_oreochromis_orthos.txt",
              candidateSNPs=cand,
              noncandidateSNPs=noncand,
              FDR=0.05,
              runs=10000,
              extension=50)
  return(y)
})

saveRDS(out, file = "snp2go_out.rds")



