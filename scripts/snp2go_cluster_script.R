library(SNP2GO)
library(Rgb)
library(here)
library(stringr)
library(data.table)

setwd("/scicore/home/schiera/gizevo30/projects/cichlids_2")

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("must have 2 arguments for Rscript command")
}

args[1] <- as.character(args[1])
args[2] <- as.character(args[2])

## args[[1]] should be mart file, args[[2]] should be phenotype
## args[1] can be "mart_export_human_go_terms_oreochromis_orthos.txt" (for everything)
## or "custom_gene_set_geneBygo_DISGENET.txt" etc

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
# 
# non_go_terms <- read.csv(file = "custom_gene_set_geneBygo.txt", header = F, sep = "\t")
# colnames(non_go_terms) <- c("Gene.name", "GO.term.accession")
# non_go_terms$GO.term.accession <- str_replace(non_go_terms$GO.term.accession, "'", "Q")
# 
# 
# human_go_terms_all <- rbind(human_go_terms, non_go_terms)
# write.table(human_go_terms, file = "mart_export_human_go_terms_oreochromis_orthos.txt", quote = F, row.names = F, sep = "\t")
# write.table(human_go_terms_all, file = "mart_export_human_go_terms_oreochromis_orthos_all.txt", quote = F, row.names = F, sep = "\t")
# 
# human_go_terms_all_sub <- human_go_terms_all[sample(1:nrow(human_go_terms_all), 5000),]
# write.table(human_go_terms_all_sub, file = "mart_export_human_go_terms_oreochromis_orthos_all_sub.txt", quote = F, row.names = F, sep = "\t")

## For a vcf file, the first colummn is the CHROM, the 2nd is the position, which I have in not a vcf file
## This should be easier than loading in some gigantic vcf file

########################################################################
#######   Modify the functions + add back to package  ##################
########################################################################

snp2goNew <- function (gtf, gff, goFile, FDR = 0.05, runs = 1e+05, candidateSNPs, 
                       noncandidateSNPs, extension = 0, verbose = TRUE, min.regions = 10) 
{
  if (!missing("gff")) {
    if (!missing("goFile")) {
      warning("Parameter 'goFile' is ignored in GFF mode")
    }
    if (!missing("gtf")) {
      warning("Paramter 'gtf' is ignored in GFF mode")
    }
    inData <- ProcessGFF(gff.file = gff, c.ranges = candidateSNPs, 
                         nc.ranges = noncandidateSNPs, updownregion = extension, 
                         verbose = verbose)
  }
  else if (!missing("gtf")) {
    if (!missing("gff")) {
      warning("Parameter 'gff' is ignored in GTF mode")
    }
    inData <- ProcessMartExport(gtf.file = gtf, mart.file = goFile, 
                                c.ranges = candidateSNPs, nc.ranges = noncandidateSNPs, 
                                region.type = "gene", updownregion = extension, 
                                verbose = verbose)
  }
  else {
    stop("No GTF or GFF file specified")
  }
  numSNPs <- list()
  numInformativeSNPs <- list()
  numSNPs[["candidates"]] <- length(candidateSNPs)
  numInformativeSNPs[["candidates"]] <- length(unique(subjectHits(inData$regions2candidates)))
  numSNPs[["noncandidates"]] <- length(noncandidateSNPs)
  numInformativeSNPs[["noncandidates"]] <- length(unique(subjectHits(inData$regions2noncandidates)))
  if (verbose) {
    message("Organizing GO categories by level")
  }
  gotermsByLevel <- ComputeLevels(inData$goterms)
  termlevel <- gotermsByLevel[["termlevel"]]
  if (verbose) {
    message("Calculating number of candidate SNPs and non-candidate SNPs\nassociated with GO terms from different ontology levels:")
  }
  snpsByLevel <- list()
  ### Only here do I need to change anything, to add the other category
  for (onto in c("MF", "BP", "CC", "other")) {
    snpsByLevel[[onto]] <- ComputeSNPsByLevel(gotermsByLevel = gotermsByLevel[[onto]], 
                                              gotermMappings = inData$go2region, numInformativeSNPs = numInformativeSNPs, 
                                              ontoName = onto, verbose = verbose)
  }
  goterm2ranges <- list()
  goterm2ranges[["regions"]] <- hash::hash()
  goterm2ranges[["candidates"]] <- hash::hash()
  goterm2ranges[["noncandidates"]] <- hash::hash()
  if (verbose) {
    message("Updating mappings of GO terms to SNPs and genomic regions\nfor inclusive analysis")
  }
  progEnvir <- new.env()
  progIncr <- GetProgressBarIncrementer(verbose = verbose, 
                                        environment = progEnvir, min = 0, max = length(inData$goterms), 
                                        style = 3)
  for (g in inData$goterms) {
    progIncr(progEnvir)
    x <- GetAllRegionsAndSNPs(g, inData$goterms, inData$go2region)
    goterm2ranges[["regions"]][[g]] <- x[["allCandidateRegions"]]
    goterm2ranges[["candidates"]][[g]] <- x[["allCandidateSNPs"]]
    goterm2ranges[["noncandidates"]][[g]] <- x[["allNoncandidateSNPs"]]
  }
  rm(progEnvir)
  if (verbose) {
    message("Testing GO terms for over-representation of candidate SNPs")
  }
  progEnvir <- new.env()
  progIncr <- GetProgressBarIncrementer(verbose = verbose, 
                                        environment = progEnvir, min = 0, max = length(inData$goterms), 
                                        style = 3)
  enrichmentPvals <- lapply(inData$goterms, function(g) {
    progIncr(progEnvir)
    if (is.na(goterm2ranges[["regions"]][[g]][1])) {
      return(NA)
    }
    GetStatistic(goterm = g, gotermsByLevel = gotermsByLevel, 
                 snpsByLevel = snpsByLevel, cc = length(goterm2ranges[["candidates"]][[g]]), 
                 cn = length(goterm2ranges[["noncandidates"]][[g]]))
  })
  rm(progEnvir)
  enrichmentPvals <- enrichmentPvals[!is.na(enrichmentPvals)]
  enrichmentPvalsAdjusted <- p.adjust(unlist(sapply(enrichmentPvals, 
                                                    "[", 2)), method = "BH")
  isSignificant <- enrichmentPvalsAdjusted < FDR
  signEnrichmentPvals <- enrichmentPvals[isSignificant]
  signGOTerms <- unlist(sapply(signEnrichmentPvals, "[", 1))
  signEnrichmentPvalsAdjusted <- enrichmentPvalsAdjusted[isSignificant]
  minRegions <- min.regions
  hasMinRegions <- unlist(sapply(signGOTerms, function(g) {
    length(reduce(inData$ranges[goterm2ranges[["regions"]][[g]]])) >= 
      minRegions
  }))
  signGOTerms <- signGOTerms[hasMinRegions]
  signEnrichmentPvals <- signEnrichmentPvals[hasMinRegions]
  signEnrichmentPvalsAdjusted <- signEnrichmentPvalsAdjusted[hasMinRegions]
  if (length(signGOTerms) == 0) {
    if (verbose) {
      message("No significant GO terms found")
    }
    result <- data.frame(GO = character(0), P = numeric(0), 
                         FDR = numeric(0), p.L = numeric(0), p.G = numeric(0), 
                         g = numeric(0), G = numeric(0), nc = numeric(0), 
                         mc = numeric(0), GO.def = character(0), child.GOs = character(0))
  }
  else {
    if (verbose) {
      message(paste("Found", length(signGOTerms), "significant GO terms."))
      message("Testing significant GO terms for local effects")
    }
    progEnvir <- new.env()
    progIncr <- GetProgressBarIncrementer(verbose = verbose, 
                                          environment = progEnvir, min = 0, max = length(signGOTerms), 
                                          style = 3)
    samplingValues <- lapply(signGOTerms, function(sig) {
      progIncr(progEnvir)
      SamplingTest(goterm = sig, iterations = runs, ranges = inData$ranges, 
                   c.ranges = candidateSNPs, nc.ranges = noncandidateSNPs, 
                   regions2candidates = inData$regions2candidates, 
                   regions2noncandidates = inData$regions2noncandidates, 
                   go2allRegions = goterm2ranges[["regions"]], 
                   go2allCandidateSNPs = goterm2ranges[["candidates"]], 
                   go2allNoncandidateSNPs = goterm2ranges[["noncandidates"]])
    })
    rm(progEnvir)
    samplingDataFrame <- data.frame(matrix(unlist(samplingValues), 
                                           byrow = TRUE, nrow = length(samplingValues)))
    colnames(samplingDataFrame) <- c("p.L", "p.G", "g")
    signData <- lapply(signGOTerms, function(g) {
      list(numTotalRegions = length(reduce(inData$ranges[goterm2ranges[["regions"]][[g]]])), 
           numCandidateSNPs = length(unique(goterm2ranges[["candidates"]][[g]])), 
           numNoncandidateSNPs = length(unique(goterm2ranges[["noncandidates"]][[g]])), 
           definition = Definition(g), significantOffspringTerms = AddOffspringInformation(g, 
                                                                                           significantGOTerms = signGOTerms))
    })
    signDataFrame <- matrix(unlist(signData), byrow = TRUE, 
                            nrow = length(signGOTerms))
    colnames(signDataFrame) <- c("G", "nc", "mc", "GO.def", 
                                 "child.GOs")
    result <- data.frame(GO = signGOTerms, P = unlist(sapply(signEnrichmentPvals, 
                                                             "[", 2)), FDR = unlist(signEnrichmentPvalsAdjusted), 
                         samplingDataFrame, signDataFrame)
    for (col in c("G", "nc", "mc")) {
      result[[col]] <- as.numeric(levels(result[[col]]))[result[[col]]]
    }
  }
  return(list(FDR = FDR, runs = runs, extension = extension, 
              candidateSNPs = numSNPs[["candidates"]], noncandidateSNPs = numSNPs[["noncandidates"]], 
              informative.candidateSNPs = numInformativeSNPs[["candidates"]], 
              informative.noncandidateSNPs = numInformativeSNPs[["noncandidates"]], 
              enriched = result, goterms = inData$goterms, regions = inData$ranges, 
              go2ranges = goterm2ranges, termlevel = termlevel))
}


environment(snp2goNew) <- asNamespace('SNP2GO')
assignInNamespace("snp2go", snp2goNew, ns = "SNP2GO")

ComputeLevelsNEW <- function (goterms) 
{
  count <- 1
  GO_Levels <- list()
  termlevel <- list()
  
  ## The following extracts all terms associated with MF, BF, and CC, and creates nested lists based on their level
  ## This is basically to speed up analysis later on I think, and to make the structure work better?
  for (onto in c("MF", "BP", "CC")) {
    Levels <- list()
    repeat {
      # getGOLevel is from clusterProfiler
      levelTerms <- getGOLevel(onto, count)
      if (length(levelTerms) == 0) {
        count <- 1
        break
      }
      Levels[[onto]] <- c(Levels[[onto]], list(levelTerms))
      count <- count + 1
    }
    # This removes terms from upper levels (1 etc) that also appear in lower levels?
    Levels[[onto]] <- RemoveFromUpperLevels(Levels[[onto]])
    for (i in 1:length(Levels[[onto]])) {
      Levels[[onto]][[i]] <- Levels[[onto]][[i]][Levels[[onto]][[i]] %in% 
                                                   goterms]
    }
    # This then adds all the lower level terms to the level above (5 -> 4, then 4 -> 3, etc, so that level 1 has all the genes)
    # This is then returned for each of the categories (MF, CC, BP)
    oo <- paste(onto, ".new", sep = "")
    Levels[[oo]] <- list()
    for (i in 1:length(Levels[[onto]])) {
      Levels[[oo]][[i]] <- character()
      for (j in i:length(Levels[[onto]])) {
        Levels[[oo]][[i]] <- c(Levels[[oo]][[i]], Levels[[onto]][[j]])
      }
      Levels[[oo]][[i]] <- unique(Levels[[oo]][[i]])
      ## This is a list, for each GO term, with it's associated level
      termlevel[Levels[[oo]][[i]]] <- i
    }
    GO_Levels[[onto]] <- Levels[[oo]]
  }
  
  ## I need to now take all 'goterms' (from inData$goterms) that don't appear in termlevel
  ## then assign them a level of 1 (doesn't matter I suppose, but probably best). I then
  ## need to created a nested list with them as well, under some other category
  ## Hopefully this is then recognized? Or should it be under one of the main ones?
  other_goterms <- goterms[!(goterms %in% names(termlevel))]
  other_list <- as.list(rep(1, length(other_goterms)))
  names(other_list) <- other_goterms
  
  # Adds other go terms to the GO_Levels object
  GO_Levels[["other"]][[1]] <- other_goterms
  
  ## Also need to append termlevel
  termlevel <- c(termlevel, other_list)
  
  return(list(MF = GO_Levels[["MF"]], CC = GO_Levels[["CC"]], 
              BP = GO_Levels[["BP"]], other = GO_Levels[["other"]], termlevel = termlevel))
}

environment(ComputeLevelsNEW) <- asNamespace('SNP2GO')
assignInNamespace("ComputeLevels", ComputeLevelsNEW, ns = "SNP2GO")


GetStatisticNEW <- function (goterm, gotermsByLevel, snpsByLevel, cc, cn) 
{
  ## This determines the level of the goterm
  ## I could add an if statement
  if(goterm %in% gotermsByLevel[["other"]][[1]]) {
    o <- "other" # Or whatever I call it above
  } else {
    o <- Ontology(goterm)
  }
  
  level <- FindLevel(goterm, gotermsByLevel[[o]])
  allCandidates <- snpsByLevel[[o]][["LengthCandidatesByLevel"]][level]
  allNoncandidates <- snpsByLevel[[o]][["LengthNoncandidatesByLevel"]][level]
  nc <- allCandidates - cc
  nn <- allNoncandidates - cn
  list(goterm, fisher.test(matrix(c(cc, cn, nc, nn), byrow = T, 
                                  nrow = 2), alternative = "g")$p.value)
}

environment(GetStatisticNEW) <- asNamespace('SNP2GO')
assignInNamespace("GetStatistic", GetStatisticNEW, ns = "SNP2GO")


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

bkgrd_snps <- GRanges(seqnames=merged$CHROM,ranges=IRanges(merged$POS,merged$POS))

########################################################################
#######   Load significant SNP lists  ##################################
########################################################################
filter_snps <- readRDS(here("sra_reads_nobackup/dump/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

if (args[2] == "pc1") {
  snps <- filter_snps[[4]]
  bkgrd_snps <- bkgrd_snps[!(merged$location %in% filter_snps[[4]]$location)]
}
if (args[2] == "pc2") {
  snps <- filter_snps[[5]]
  bkgrd_snps <- bkgrd_snps[!(merged$location %in% filter_snps[[5]]$location)]
}
if (args[2] == "tr") {
  snps <- filter_snps[[6]]
  bkgrd_snps <- bkgrd_snps[!(merged$location %in% filter_snps[[6]]$location)]
}

snps <- GRanges(seqnames = snps$CHROM, ranges = IRanges(snps$POS, snps$POS))

rm(dfs_all)

########################################################################
#######   Run the analysis per candidate SNP list  E####################
########################################################################

## Run for pc1
y <- SNP2GO:::snp2go(gtf = "genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf",
                    goFile = args[1],
                    candidateSNPs = snps,
                    noncandidateSNPs = bkgrd_snps,
                    FDR = 0.05,
                    runs = 10000,
                    extension = 50)

saveRDS(y, file = paste("snp2go_out", args[2], args[1], ".rds", sep = "_"))
saveRDS(y$enriched, file = paste("snp2go_out", args[2], args[1], "enriched.rds", sep = "_"))



