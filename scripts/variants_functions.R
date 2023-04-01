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
        start <- as.numeric(x[grep("start|POS|coord", colnames)])-window
        stop <- as.numeric(x[grep("end|POS|coord", colnames)])+window
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
    
    # the next two lines takes ~0.2 sec for every apply loop, but can be done outside of this function, or minimally outside of the apply
    #gene_table <- gene_table[which(gene_table[1] == chrom), ] # 1 not 2
    #gene_table_2 <- gene_table_2[which(gene_table_2[1] == chrom), ] # 1 not 2
    
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
      if ("pvalue" %in% colnames) {
        if (windowed == TRUE) {
          gene_ids$CHROM <- as.character(x[grep("chr", colnames)])
          gene_ids$START <- as.numeric(x[grep("start", colnames)])
          gene_ids$END <- as.numeric(x[grep("end", colnames)])
          gene_ids$ps <- as.numeric(x[grep("pvalue", colnames)[1]])
        } else {
          gene_ids$CHROM <- as.character(x[grep("chr", colnames)])
          gene_ids$POS <- as.numeric(x[grep("coord", colnames)])
          gene_ids$ps <- as.numeric(x[grep("pvalue", colnames)[1]])
        }
      }
      
      if (windowed == TRUE) {
        require(DescTools)
        within_gene <- apply(gene_ids, 1, function(x) ifelse( c(x["start"], x["end"]) %overlaps% c(x["START"], x["END"]), "within",  ifelse(x["END"] < x["start"], "upstream", ifelse(x["START"] > x["end"], "downstream", NA))) )
        gene_ids$distance_to_gene <- ifelse(within_gene == "upstream", gene_ids$start - gene_ids$START, ifelse(within_gene == "downstream", gene_ids$END - gene_ids$end, 0))
        gene_ids$closest <- ifelse(gene_ids$distance_to_gene == 0, "yes", ifelse(gene_ids$distance_to_gene == min(gene_ids$distance_to_gene), "yes", "no"))
      } else {
        within_gene <- ifelse(data.table::between(gene_ids$POS,gene_ids$start, gene_ids$end), "within", ifelse(gene_ids$POS < gene_ids$start, "upstream","downstream"))
        gene_ids$distance_to_gene <- ifelse(within_gene == "upstream", gene_ids$start - gene_ids$POS, ifelse(within_gene == "downstream", gene_ids$POS - gene_ids$end, 0))
        gene_ids$closest <- ifelse(gene_ids$distance_to_gene == 0, "yes", ifelse(gene_ids$distance_to_gene == min(gene_ids$distance_to_gene), "yes", "no"))
      }
    }
    return(gene_ids)
  }))
  
  return(results)
}

## This function annotates whether a coding change is in a wobble position

wobbleAA <- function(results.table = exon_genes[1:20,], gene_table = gtf4, window = 0) {
  require(stringr)
  colnames <- colnames(results.table)
  results.table <- as.data.frame(results.table)
  
  results <- bind_rows(apply(results.table, 1, function(x) {
    # names(x) <- colnames
    chrom <- as.character(x[grep("CHROM|chr", colnames)]) # change to character
    if (any(grepl(".1", gene_table[,1]))) {
      chrom <- paste(chrom, ".1", sep = "")
    }
    start <- as.numeric(x[grep("POS|coord", colnames)])-window
    stop <- as.numeric(x[grep("POS|coord", colnames)])+window
    cols <- c("seqname", "feature", "start", "end", "strand", "frame", "Dbxref", "gene", "gene_biotype", "product")
    
    gene_table <- gene_table[which(gene_table[1] == chrom), ] # 1 not 2
    # gene_table_2 <- gene_table_2[which(gene_table_2[1] == chrom), ] # 1 not 2
    gene_ids.1 <- gene_table[which(start > gene_table$start & start < gene_table$end), cols]
    gene_ids.2 <- gene_table[which(stop > gene_table$start & stop < gene_table$end), cols]
    gene_ids.3 <- gene_table[which(start < gene_table$start & stop > gene_table$end), cols]
    gene_ids <- rbind(gene_ids.1, gene_ids.2, gene_ids.3)
    gene_ids <- unique(gene_ids)
    if (nrow(gene_ids) > 0) {
      gene_ids$Dbxref <- substr(gene_ids$Dbxref, start = 8, stop = 50)
      if ("PGLSps" %in% colnames) {
        gene_ids$CHROM <- as.character(x[grep("CHROM|chr", colnames)])
        gene_ids$POS <- as.numeric(x[grep("POS|coord", colnames)])
        gene_ids$ps <- as.numeric(x[grep("ps", colnames)[1]])
        gene_ids$PGLSps <- as.numeric(x[grep("PGLSps", colnames)])
      }
      if ("ps" %in% colnames) {
        gene_ids$CHROM <- as.character(x[grep("CHROM|chr", colnames)])
        gene_ids$POS <- as.numeric(x[grep("POS|coord", colnames)])
        gene_ids$ps <- as.numeric(x[grep("ps", colnames)[1]])
      }
      gene_ids$location <- paste(gene_ids$CHROM, gene_ids$POS, sep = ":")
      gene_ids$mut_pos <- ifelse(gene_ids$strand == "+", gene_ids$POS - gene_ids$start, gene_ids$end - gene_ids$POS)
      gene_ids$wobble <- unlist(apply(gene_ids, 1, function(y) {
        if(y[5] == "+") {
          if(as.numeric(y[6]) == 0){
            w <- ((as.numeric(y[12]) - as.numeric(y[3]) + 1) / 3)%%1==0
          } else {
            w <- ((as.numeric(y[12]) - as.numeric(y[3]) + 1 - as.numeric(y[6])) / 3)%%1==0
          }
        }
        if(y[5] == "-") {
          if(as.numeric(y[6]) == 0){
            w <- ((as.numeric(y[4]) - as.numeric(y[12]) + 1) / 3)%%1==0
          } else {
            w <- ((as.numeric(y[4]) - as.numeric(y[12]) + 1 - as.numeric(y[6])) / 3)%%1==0
          }        }
        return(w)
      }))
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


# Pull all SNPS from a given window, or a window around a gene (upstream and downstream)
# Have it such that if I give it a gene, it looks in the gtf file, and adds the window_padding, otherwise uses the window values given
# Probably makes sense to fast read in the pvalues file (for whichever chromo), rather then giving it the full pvalues df
# This way I can also give it a list of genes, and for each gene run this (after first finding their locations and sorting by location)

getLociEns <- function(gene_list = list("kcnc2","LOC100698148", "LOC100701437"), window_padding = NULL, start = NULL, stop = NULL, gtf_file = gtf2, directory = "/Volumes/BZ/Scientific Data/RG-AS04-Data01/Cichlid-genomes/cichlid_sleep_gwas/sra_reads_nobackup/alleleFreqs/", gwas.datasets = NULL, pgls.datasets = NULL, pheno_to_test = "day_night_dif") {
  if (is.list(gene_list) != TRUE) {
    stop("gene_list must be a list of genes")
  }
  if (dir.exists(directory) != TRUE) {
    stop("directory does not exist")
  }
  
  gene_locations <- gtf_file[gtf_file$gene %in% gene_list,]
  gene_locations <- gene_locations[order(gene_locations$seqname),] # I think this is not required (they are already in order)
  chromosomes <- substr(gene_locations$seqname, start = 1, stop = 9)
  gene_locations$chromosome <- chromosomes
  
  # Should I load all of the chromosomes that are in gene_locations? Maybe look for them in the environment? Oh, I know, I can supply a list of chromosomes, OR load them
  
  if (is.null(gwas.datasets) == TRUE | is.null(pgls.datasets) == TRUE) {
    # Load the datasets that are dataset_type
    if (pheno_to_test == "peak") {
      datasets <- list.files(directory)[grep(pheno_to_test, list.files(directory))]
      datasets <- datasets[!(grepl("amplitude", datasets))]
      datasets <- datasets[grep(".txt", datasets)]
    } 
    if (pheno_to_test == "day_night_dif") {
      datasets <- list.files(directory)[grep(pheno_to_test, list.files(directory))]
      datasets <- datasets[!(grepl("-species", datasets))]
      datasets <- datasets[grep(".txt", datasets)]
    }  else {
      datasets <- list.files(directory)[grep(pheno_to_test, list.files(directory))]
      datasets <- datasets[grep(".txt", datasets)]
    }
    
    datasets <- datasets[grep(paste(chromosomes,collapse="|"), datasets)]
    datasets.gwas <- datasets[grep("_piVals", datasets)]
    datasets.pgls <- datasets[grep("_PGLSpiVals", datasets)]
    
    gwas.datasets <- lapply(datasets.gwas, function(x) fread(x, showProgress = F, drop = "V1"))
    names(gwas.datasets) <- chromosomes
    pgls.datasets <- lapply(datasets.pgls, function(x) fread(x, showProgress = F, drop = "V1"))
    names(pgls.datasets) <- chromosomes
  }
  
  snps_output <- apply(gene_locations, 1, function(x) {
    gwas.result <- gwas.datasets[[as.character(x[60])]]
    gwas.result <- gwas.result[as.numeric(gwas.result[["coord"]]) > as.numeric(x[4])-as.numeric(window_padding) & as.numeric(gwas.result[["coord"]]) < as.numeric(x[5])+as.numeric(window_padding)]
    
    pgls.result <- pgls.datasets[[as.character(x[60])]]
    pgls.result <- pgls.result[as.numeric(pgls.result[["coord"]]) > as.numeric(x[4])-as.numeric(window_padding) & as.numeric(pgls.result[["coord"]]) < as.numeric(x[5])+as.numeric(window_padding)]
    merged <- merge(gwas.result, pgls.result)
    return(merged)
  })
  names(snps_output) <- gene_locations$gene
  
  return(snps_output)
}




## This one needs to run a bunch of system commands to query the orthoDB database
## Takes a Dbxref ID as input, and outputs a list of vectors (matching orthologs)
## These are better then the above, but miss some and should be supplemented
## This can be l/apply'd across a list of gene_ids

generateOrthoDBorthologs <- function(dbxrefs = Dbxref, gene_ids = gene_id, species = c("zebrafish", "mouse")) {
  if (species == "zebrafish") {
    group <- "Actinopterygii"
    species <- "danio"
  }
  if (species == "mouse") {
    group <- "Vertebrata"
    species <- "Mus musculus"
  }
  
  tilapia_orthos <- list()
  # Search orthoDB for gene ID and return orthoID
  # This doesn't work if it is a LOC gene (eg LOC102076095) and returns null
  system(paste("curl 'https://data.orthodb.org/current/search?query=", dbxrefs, "' -L -o search.dat", sep = ""))
  orthoID_json <- fromJSON(file = "search.dat")
  
  # If the ID doesn't find anything, try LOC-ID
  if (is.null(orthoID_json[[4]][grep(group, orthoID_json[[4]])])) {
    system(paste("curl 'https://data.orthodb.org/current/search?query=LOC", dbxrefs, "' -L -o search.dat", sep = ""))
    orthoID_json <- fromJSON(file = "search.dat")
  }
  
  for (k in grep(group, orthoID_json[[4]])) {
    orthoID_json_act <- orthoID_json[[4]][k][[1]]$id
    
    # Use orthoID to search for gene IDs
    system(paste("curl 'https://data.orthodb.org/current/orthologs?id=", orthoID_json_act,"' -L -o orthologs.dat", sep = ""))
    orthologs <- fromJSON(file = "orthologs.dat")
    
    # this gives me just danio (could also do mus musculus)
    if (any(grepl(species, orthologs[[1]]))) {
      orthologs <- orthologs[[1]][grep(species, orthologs[[1]])]
      orthos <- list()
      for (j in 1:length(orthologs[[1]]$genes)) {
        orthos[j] <- orthologs[[1]]$genes[[j]]$gene_id$id
      }
      
      tilapia_orthos[[dbxrefs]] <- c(tilapia_orthos[[dbxrefs]], orthos)
      
    } else { # This is for the occasion where orthoDB can't find the LOC or Dbxref ID #, but the gene name is useful
      
      system(paste("curl 'https://data.orthodb.org/current/search?query=", gene_ids, "' -L -o search.dat", sep = ""))
      orthoID_json <- fromJSON(file = "search.dat")
      
      for (p in grep(group, orthoID_json[[4]])) {
        orthoID_json_act <- orthoID_json[[4]][p][[1]]$id
        
        system(paste("curl 'https://data.orthodb.org/current/orthologs?id=",orthoID_json_act,"' -L -o orthologs.dat", sep = ""))
        orthologs <- fromJSON(file = "orthologs.dat")
        
        if (any(grepl(species, orthologs[[1]]))) {
          orthologs <- orthologs[[1]][grep(species, orthologs[[1]])]
          orthos <- list()
          for (j in 1:length(orthologs[[1]]$genes)) {
            orthos[j] <- orthologs[[1]]$genes[[j]]$dbxrefs$id
          }
          tilapia_orthos[[dbxrefs]] <- c(tilapia_orthos[[dbxrefs]], orthos)
        }
      }
    }
  }
  return(tilapia_orthos)
}



## This one needs to recursively find all the paths
## Then load and extract the ortholog from the csv
## This output has duplicate entries (because of transcript variants)

generateNCBIorthologs <- function(directory = dir) {
  files <- list.files(directory, recursive = T, pattern = "data_table.tsv")
  output <- lapply(files, function(x) {
    df <- read.csv(paste(directory, x, sep = ""), header = T, sep = "\t")
    df$gene_id_oreo <- str_split(x, pattern = "/")[[1]][1]
    return(df)
  })
  return(output)
}



