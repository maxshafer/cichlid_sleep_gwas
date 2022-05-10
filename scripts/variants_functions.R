# getBM_MS retrieves ensembl gene symbols given a chromosome and Start Stop positions
# Modification of biomaRt get_BM function
# Used to retrieve genes associated with FST windows

getBM_MS <- function(x, window = 25000) {
	require(biomaRt)
	CHROM <- x[1]
	START <- as.numeric(x[2])-window
	STOP <- as.numeric(x[3])+window
	genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = c('chromosome_name', 'start', 'end'), values = list(CHROM, max(START, 0), STOP), mart = ensembl)
	if (nrow(genes) > 0)
		genes$SNP <- x[7]
	return(genes)
}

# Above is very slow, has to make new query every time
# getGenesEns takes genes (my_gene) from ensembl gene annotations (downloaded) (using the refGenome package)
# getGenesEns <- function(x, gene_table = gtf2, window = 25000, gene_table_2 = gtf3, windowed = TRUE) {
#   require(stringr)
# 	chrom <- as.character(x[1]) # change to character
# 	if (windowed == TRUE) {
# 	  start <- as.numeric(x[2])-window
# 	  stop <- as.numeric(x[3])+window
# 	  cols <- c(1,3:5,10,20:21,24)
# 	} 
# 	if (windowed == FALSE) {
# 	  start <- as.numeric(x[2])-window
# 	  stop <- as.numeric(x[2])+window
# 	  cols <- c(1,3:5,10,20:21,24) 
# 	}
# 	gene_table <- gene_table[which(gene_table[1] == chrom), ] # 1 not 2
# 	gene_table_2 <- gene_table_2[which(gene_table_2[1] == chrom), ] # 1 not 2
# 	gene_ids.1 <- gene_table[which(start > gene_table$start & start < gene_table$end), cols]
# 	gene_ids.2 <- gene_table[which(stop > gene_table$start & stop < gene_table$end), cols]
# 	gene_ids.3 <- gene_table[which(start < gene_table$start & stop > gene_table$end), cols]
# 	gene_ids <- rbind(gene_ids.1, gene_ids.2, gene_ids.3)
# 	gene_ids <- unique(gene_ids)
# 	if (nrow(gene_ids) > 0) {
# 	  gene_ids$product <- gene_table_2$product[unlist(lapply(gene_ids$Dbxref, function(y) grep(y, gene_table_2$Dbxref)[1]))]
# 	  gene_ids$Dbxref <- substr(gene_ids$Dbxref, start = 8, stop = 50)
# 	  gene_ids$SNP <- as.numeric(x[4])
# 	  gene_ids$CHROM <- as.character(x[1])
# 	  gene_ids$POS <- as.numeric(x[2])
# 	  gene_ids$WEIR_AND_COCKERHAM_FST <- as.numeric(x[3])
# 	  gene_ids$zscore <- as.numeric(x[4])
# 	  gene_ids$FisherZ <- as.numeric(x[5])
# 	  within_gene <- ifelse(between(gene_ids$POS,gene_ids$start, gene_ids$end), "within", ifelse(gene_ids$POS < gene_ids$start, "upstream","downstream"))
# 	  gene_ids$distance_to_gene <- ifelse(within_gene == "upstream", gene_ids$start - gene_ids$POS, ifelse(within_gene == "downstream", gene_ids$POS - gene_ids$end, 0))
# 	  gene_ids$closest <- ifelse(gene_ids$distance_to_gene == 0, "yes", ifelse(gene_ids$distance_to_gene == min(gene_ids$distance_to_gene), "yes", "no"))
# 	}
# 	return(gene_ids)
# }


# This one works for more situations, but needs to be run on the results table, rather than a row from that table (provided by an apply loop)
getGenesEns <- function(results.table, gene_table = gtf2, window = 25000, gene_table_2 = gtf3, windowed = TRUE, combined = FALSE) {
  require(stringr)
  colnames <- colnames(results.table)
  results.table <- as.data.frame(results.table)
  
  results <- bind_rows(apply(results.table, 1, function(x) {
    # names(x) <- colnames
    if (combined == FALSE) {
      chrom <- as.character(x[grep("CHROM|chr", colnames)]) # change to character
      if (any(grepl(".1", gene_table[,1]))) {
        chrom <- paste(chrom, ".1", sep = "")
      }
      if (windowed == TRUE) {
        start <- as.numeric(x[grep("POS|coord", colnames)])-window
        stop <- as.numeric(x[grep("POS|coord", colnames)])+window
        cols <- c("seqname", "feature", "start", "end", "Dbxref", "gene", "gene_biotype", "product")
      } 
      if (windowed == FALSE) {
        start <- as.numeric(x[grep("POS|coord", colnames)])-window
        stop <- as.numeric(x[grep("POS|coord", colnames)])+window
        cols <- c("seqname", "feature", "start", "end", "Dbxref", "gene", "gene_biotype", "product") 
      }
    }
    if (combined == TRUE) {
      chrom <- as.character(strsplit(x[grep("CHROM_POS", colnames)], ":")[[1]][[1]])
      pos <- as.numeric(strsplit(x[grep("CHROM_POS", colnames)], ":")[[1]][[2]])
      start <- as.numeric(pos)-window
      stop <- as.numeric(pos)+window
      cols <- c("seqname", "feature", "start", "end", "Dbxref", "gene", "gene_biotype", "product")
    }
    
    gene_table <- gene_table[which(gene_table[1] == chrom), ] # 1 not 2
    gene_table_2 <- gene_table_2[which(gene_table_2[1] == chrom), ] # 1 not 2
    gene_ids.1 <- gene_table[which(start > gene_table$start & start < gene_table$end), cols]
    gene_ids.2 <- gene_table[which(stop > gene_table$start & stop < gene_table$end), cols]
    gene_ids.3 <- gene_table[which(start < gene_table$start & stop > gene_table$end), cols]
    gene_ids <- rbind(gene_ids.1, gene_ids.2, gene_ids.3)
    gene_ids <- unique(gene_ids)
    if (nrow(gene_ids) > 0) {
      gene_ids$product <- gene_table_2$product[unlist(lapply(gene_ids$Dbxref, function(y) grep(y, gene_table_2$Dbxref)[1]))]
      gene_ids$Dbxref <- substr(gene_ids$Dbxref, start = 8, stop = 50)
      if ("WEIR_AND_COCKERHAM_FST" %in% colnames) {
        gene_ids$CHROM <- as.character(x[grep("CHROM", colnames)])
        gene_ids$POS <- as.numeric(x[grep("POS", colnames)])
        gene_ids$SNP <- as.numeric(x[4])
        gene_ids$WEIR_AND_COCKERHAM_FST <- as.numeric(x[grep("WEIR_AND_COCKERHAM_FST", colnames)])
        gene_ids$zscore <- as.numeric(x[4])
        gene_ids$FisherZ <- as.numeric(x[5])
      }
      if ("p_wald" %in% colnames) {
        gene_ids$CHROM <- as.character(x[grep("CHROM", colnames)])
        gene_ids$POS <- as.numeric(x[grep("POS", colnames)])
        gene_ids$p_wald <- as.numeric(x[grep("p_wald", colnames)])
      }
      if ("p_WALD" %in% colnames) {
        gene_ids$CHROM <- as.character(chrom)
        gene_ids$POS <- as.numeric(pos)
        gene_ids$p_WALD <- as.numeric(x[grep("p_WALD", colnames)])
        gene_ids$WEIR_FST <- as.numeric(x[grep("WEIR_FST", colnames)])
        gene_ids$CHROM_POS <- as.character(x[grep("CHROM_POS", colnames)])
      }
      if ("PGLSps" %in% colnames) {
        gene_ids$CHROM <- as.character(x[grep("chr", colnames)])
        gene_ids$POS <- as.numeric(x[grep("coord", colnames)])
        gene_ids$ps <- as.numeric(x[grep("ps", colnames)[1]])
        gene_ids$PGLSps <- as.numeric(x[grep("PGLSps", colnames)])
      }
      within_gene <- ifelse(data.table::between(gene_ids$POS,gene_ids$start, gene_ids$end), "within", ifelse(gene_ids$POS < gene_ids$start, "upstream","downstream"))
      gene_ids$distance_to_gene <- ifelse(within_gene == "upstream", gene_ids$start - gene_ids$POS, ifelse(within_gene == "downstream", gene_ids$POS - gene_ids$end, 0))
      gene_ids$closest <- ifelse(gene_ids$distance_to_gene == 0, "yes", ifelse(gene_ids$distance_to_gene == min(gene_ids$distance_to_gene), "yes", "no"))
    }
    return(gene_ids)
  }))
  
  return(results)
}

# getGenesEns <- function(x, gene_table = my_gene, window = 25000) {
	# require(refGenome)
	# chrom <- x[1]
	# start <- as.numeric(x[2])-window
	# stop <- as.numeric(x[3])+window
	# fst <- seq(start, stop)
	# gene_ids <- apply(gene_table, 1, function(x) any(seq(x[4], x[5]) %in% fst)) # 30 mins per fst window?!
	# return(gene_ids)
# }

# readFST returns a data.frame of FST values
# Can use an fst cutoff value (0.1 for indels, 0.2 for snps), or z-transform fst values

readFST <- function(x, fisherZ = TRUE, variants = "windowed", gemma = FALSE) {
	require(DescTools)
  if (gemma == TRUE) {
    fst.table <- read.table(x, header = T, stringsAsFactors = FALSE, sep = "\t")
    colnames <- colnames(fst.table)
    colnames[grep("chr", colnames)] <- "CHROM"
    colnames[grep("ps", colnames)] <- "POS"
    colnames(fst.table) <- colnames
    fst.table$SNP <- c(1:(nrow(fst.table)))
  } else {
    fst.table <- read.table(x, header = T, stringsAsFactors = FALSE)
    fst.table <- fst.table[complete.cases(fst.table),]
    fst.table$SNP <- c(1:(nrow(fst.table)))
    if (variants == "windowed") {
      fst.table$zscore <- (fst.table$MEAN_FST-mean(fst.table$MEAN_FST))/sd(fst.table$MEAN_FST)
      if (fisherZ == TRUE) {
        fst.table$FisherZWeighted <- FisherZ(fst.table$WEIGHTED_FST)
        fst.table$FisherZMean <- FisherZ(fst.table$MEAN_FST)
        fst.table <- fst.table[!(fst.table$FisherZMean == Inf),]
        fst.table <- fst.table[!(fst.table$FisherZMean == -Inf),]
      }
    }
    if (variants == "perbase") {
      fst.table <- fst.table[fst.table$WEIR_AND_COCKERHAM_FST > 0.1, ]
      fst.table$zscore <- (fst.table$WEIR_AND_COCKERHAM_FST-mean(fst.table$WEIR_AND_COCKERHAM_FST))/sd(fst.table$WEIR_AND_COCKERHAM_FST)
      if (fisherZ == TRUE) {
        fst.table$FisherZ <- FisherZ(fst.table$WEIR_AND_COCKERHAM_FST)
        fst.table <- fst.table[!(fst.table$FisherZ == Inf),]
        fst.table <- fst.table[!(fst.table$FisherZ == -Inf),]
      }
    }
  }
	fst.table$CHROM <- paste(fst.table$CHROM, ".1", sep = "")
	fst.table$CHROM_POS <- paste(fst.table$CHROM, fst.table$POS, sep = ":")
	return(fst.table)
}

# Determine fdr cutoffs for weir.fst tables (returned by readFST)

fdrcutoff <- function(x, cutoff = 0.05, statistic = "zscore") {
  if (!(statistic %in% c("zscore", "pvalue"))) {
    print("Statistic must be zscore or pvalue")
    stop()
  }
	require(fdrtool)
	require(dplyr)
  if (statistic == "zscore") {
    fdr <- fdrtool(x$zscore, statistic = "normal", plot = T, verbose = F)
  }
  if (statistic == "pvalue") {
    fdr <- fdrtool(x$p_wald, statistic = "pvalue", plot = T, verbose = F)
  }
	weir <- x[tibble(fdr[[3]]) < cutoff,]
	return(weir)
}