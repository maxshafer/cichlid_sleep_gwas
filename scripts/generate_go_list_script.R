library(RDAVIDWebService)
library(org.Hs.eg.db)

ortholog_df <- readRDS("GTF_human_orthologs_all.rds")

## OK, for each of these data.frames, run DAVID analysis on the gene lists
ortholog_df <- ortholog_df[ortholog_df$human_ortholog != "",]
# ortholog_df <- ortholog_df[!duplicated(ortholog_df$human_ortholog),]

dbi_out <- select(org.Hs.eg.db, ortholog_df$human_ortholog, "ENTREZID", "SYMBOL")

ortholog_df$human_ENTREZID <- dbi_out$ENTREZID[match(ortholog_df$human_ortholog, dbi_out$SYMBOL)]

categories <- c("KEGG_PATHWAY", "REACTOME_PATHWAY", "DISGENET", "GAD_DISEASE", "UP_TISSUE") # c("UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY", "DISGENET", "GAD_DISEASE", "UP_TISSUE")

markers <- ortholog_df$human_ENTREZID
david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
setAnnotationCategories(david, categories)
result <- addList(david, markers[1:3500], idType = "ENTREZ_GENE_ID", listName = "permutated_snps", listType = "Gene")
out1 <- getFunctionalAnnotationTable(david)
out1_membership <- lapply(seq_along(out1@Membership), function(i) {
  rownames(out1@Membership[[i]]) <- out1@Genes$ID
  colnames(out1@Membership[[i]]) <- paste(names(out1@Membership)[i], colnames(out1@Membership[[i]]), sep = ":")
  return(out1@Membership[[i]])
})
out1_membership <- Reduce(cbind, out1_membership)

david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
setAnnotationCategories(david, categories)
result <- addList(david, markers[3501:7000], idType = "ENTREZ_GENE_ID", listName = "permutated_snps", listType = "Gene")
out2 <- getFunctionalAnnotationTable(david)
out2_membership <- lapply(seq_along(out2@Membership), function(i) {
  rownames(out2@Membership[[i]]) <- out2@Genes$ID
  colnames(out2@Membership[[i]]) <- paste(names(out2@Membership)[i], colnames(out2@Membership[[i]]), sep = ":")
  return(out2@Membership[[i]])
})
out2_membership <- Reduce(cbind, out2_membership)

david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
setAnnotationCategories(david, categories)
result <- addList(david, markers[7001:10500], idType = "ENTREZ_GENE_ID", listName = "permutated_snps", listType = "Gene")
out3 <- getFunctionalAnnotationTable(david)
out3_membership <- lapply(seq_along(out3@Membership), function(i) {
  rownames(out3@Membership[[i]]) <- out3@Genes$ID
  colnames(out3@Membership[[i]]) <- paste(names(out3@Membership)[i], colnames(out3@Membership[[i]]), sep = ":")
  return(out3@Membership[[i]])
})
out3_membership <- Reduce(cbind, out3_membership)

david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
setAnnotationCategories(david, categories)
result <- addList(david, markers[10500:length(markers)], idType = "ENTREZ_GENE_ID", listName = "permutated_snps", listType = "Gene")
out4 <- getFunctionalAnnotationTable(david)
out4_membership <- lapply(seq_along(out4@Membership), function(i) {
  rownames(out4@Membership[[i]]) <- out4@Genes$ID
  colnames(out4@Membership[[i]]) <- paste(names(out4@Membership)[i], colnames(out4@Membership[[i]]), sep = ":")
  return(out4@Membership[[i]])
})
out4_membership <- Reduce(cbind, out4_membership)

## Combine them together?
out1_membership <- as.data.frame(out1_membership)
out2_membership <- as.data.frame(out2_membership)
out3_membership <- as.data.frame(out3_membership)
out4_membership <- as.data.frame(out4_membership)

out_membership <- bind_rows(out1_membership, out2_membership)
out_membership <- bind_rows(out_membership, out3_membership)
out_membership <- bind_rows(out_membership, out4_membership)

rm(out1_membership, out2_membership, out3_membership, out4_membership)

## OK, I now need to somehow convert this to a dataframe, with rows for each GO term, and a column featuring all of the genes that are a part of it
## I should then sanity check to make sure there is a correct ## of GO terms (which will be diff depending on background etc)

## This works, but is returning some weird row_names
## Probably because of overlaps in some go analyses, that then don't match across rows
## Best is to maybe just add all the genes, and then make those lists unique

go_description <- vapply(strsplit(colnames(out_membership),"~"), `[`, 2, FUN.VALUE=character(1))

## Gene names need to be oreochromis gene names
genes <- rownames(out_membership)
genes <- vapply(strsplit(genes,"\\.{2,}"), `[`, 1, FUN.VALUE=character(1))

## Map back to oreochromis gene names? Needs to be done in the lapply below
concatenated_genes <- lapply(c(1:ncol(out_membership)), function(i) {
  names <- unique(genes[which(out_membership[,i])])
  names <- ortholog_df$oreochromis_gene[match(names, ortholog_df$human_ENTREZID)]
  # names <- paste(names, collapse = " ")
  return(names)
})

df <- data.frame(term = colnames(out_membership), description = go_description, genes = unlist(lapply(concatenated_genes, function(x) paste(x, collapse = " "))))

write.table(df, file = "gowinda/custom_gene_set_noGO.txt", quote = F, col.names = F, row.names = F, sep = "\t")


## Make total and candidate SNP files

filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

filter_snps <- lapply(filter_snps[4:6], function(x) {
  x <- x[,c(1,2)]
  x <- x[!duplicated(x[,c(1:2)]),]
})

write.table(filter_snps[[1]], file = "gowinda/filter_snps_pc1_gowinda.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(filter_snps[[2]], file = "gowinda/filter_snps_pc2_gowinda.txt", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(filter_snps[[3]], file = "gowinda/filter_snps_tr_gowinda.txt", quote = F, col.names = F, row.names = F, sep = "\t")

## Make gene by GO term df
## I have a list, which is for each term, with the associated genes (oreochromis ones)
names(concatenated_genes) <- colnames(out_membership)

gene_by_go <- lapply(seq_along(concatenated_genes), function(i) {
  
  df <- data.frame(ids = concatenated_genes[[i]], term = names(concatenated_genes)[i])
  
})

gene_by_go <- Reduce(rbind, gene_by_go)

write.table(gene_by_go, file = "custom_gene_set_geneBygo.txt", quote = F, col.names = F, row.names = F, sep = "\t")




