library(stringr)
library(ggplot2)
library(data.table)

setwd("/scicore/home/schiera/gizevo30/projects/cichlids_2")

########################################################################
#######   Defining trailing arguments  #################################
########################################################################

# Arguments would be the alleleFreq table(s) and the snpEff table, as well as the pvalue cutoff
# Can I just give it the chromosome id, then use that to pull all the tables?
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("must have 1 arguments for Rscript command")
}

percentile <- as.numeric(args[1])

########################################################################
#######   Defining  arguments  #########################################
########################################################################

files <- list.files("sra_reads_nobackup/combined_ann/")

comparison <- c("55-species", "58-species", "spd_60-species", "peak_dawn", "peak_dusk", "total_rest")

# ########################################################################
# #######   Load in data  ################################################
# ########################################################################
# 
# # Generate a list of lists, upper level is per chromosome, lower level is per comparison
# per_chr <- lapply(files, function(file) {
#   
#   # read in the chromosome
#   df <- fread(paste("sra_reads_nobackup/combined_ann/", file, sep = ""), showProgress = T)
#   
#   comp_list <- lapply(comparison, function(comp) {
#     columns <- c(c(2:3,16:21, as.numeric(grep(comp, colnames(df)))))
#     comp_df <- df[, ..columns] 
#     colnames(comp_df) <- c("CHROM","POS","REF","ALT","ANN_GENE","ANN_IMPACT","ANN_EFFECT","ANN_DISTANCE","piVals","PGLSpiVals")
#     return(comp_df)
#   })
#   
#   names(comp_list) <- comparison
#   
#   return(comp_list)
#   
# })
# 
# ########################################################################
# #######   Merge and filter  ############################################
# ########################################################################
# 
# per_chr <- lapply(comparison, function(comp) {
#   
#   merged <- lapply(per_chr, function(x) x[[comp]])
#   merged <- Reduce(rbind, merged)
#   merged$location <- paste(merged$CHROM, merged$POS, sep = ":")
#   
#   
#   pi <- unique(merged[,c("piVals", "location")]) %>% slice_min(order_by = piVals, prop = percentile) %>% pull(location)
#   pgls <- unique(merged[,c("PGLSpiVals", "location")]) %>% slice_min(order_by = PGLSpiVals, prop = percentile) %>% pull(location)
#   
#   merged <- merged[merged$location %in% c(pi, pgls),]
#   return(merged)
# })
# 
# saveRDS(per_chr, file = paste("sra_reads_nobackup/combined_ann/filter_SNPs_perchr", percentile, "percentile", sep = "_"))

########################################################################
#######   Load/merge/filter per comparison  ############################
########################################################################

per_chr <- lapply(comparison, function(comp) {
  
  gwas.datasets <- lapply(files, function(x) {
    df <- fread(paste("sra_reads_nobackup/combined_ann/", x, sep = ""), showProgress = T)
    columns <- c(c(2:3,16:21, as.numeric(grep(comp, colnames(df)))))
    df <- df[, ..columns] # Change last two numbers to modify which comparison to keep and filter by
    colnames(df) <- c("CHROM","POS","REF","ALT","ANN_GENE","ANN_IMPACT","ANN_EFFECT","ANN_DISTANCE","piVals","PGLSpiVals")
    return(df)
  })
  
  merged <- Reduce(rbind, gwas.datasets)
  merged$location <- paste(merged$CHROM, merged$POS, sep = ":")
  
  pi <- unique(merged[,c("piVals", "location")]) %>% slice_min(order_by = piVals, prop = percentile) %>% pull(location)
  pgls <- unique(merged[,c("PGLSpiVals", "location")]) %>% slice_min(order_by = PGLSpiVals, prop = percentile) %>% pull(location)
  
  merged <- merged[merged$location %in% c(pi, pgls),]
  
  return(merged)
  
})


saveRDS(per_chr, file = paste("sra_reads_nobackup/combined_ann/filter_SNPs_perchr", percentile, "percentile", sep = "_"))


