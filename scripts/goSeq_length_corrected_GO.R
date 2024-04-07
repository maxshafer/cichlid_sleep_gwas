library(goseq)
library(Rgb)
library(here)
library(stringr)

## Load in gtf
gtf <- read.gtf("genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf")
gtf <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]
gtf$gene_length <- gtf$end - gtf$start

### OK, let's use goseq to run length corrected GO analysis
### I will also have to supply custom GO term mappings, but I should be able to get these! (minimally by running DAVID with the full gene list?, or even just using DAVID)

## Load and append gene length to the ortholog's df
ortholog_df <- readRDS("GTF_human_orthologs_all.rds")
ortholog_df$gene_length <- gtf$gene_length[match(ortholog_df$oreochromis_gene, gtf$gene_id)]

## Get human orthologs for filter snps

filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

filter_snps <- lapply(filter_snps, function(x) {
  x$human_ortholog <- ortholog_df$human_ortholog[match(x$ANN_GENE, ortholog_df$oreochromis_gene)]
  x$gene_length <- ortholog_df$gene_length[match(x$ANN_GENE, ortholog_df$oreochromis_gene)]
  
  # I want to know how many times each gene name appears, in a list of only those with unique location and gene name combos
  cont_table <- x[!duplicated(x[,11:12]),]
  cont_table <- table(cont_table$human_ortholog)
  x$number_pres <- cont_table[match(x$human_ortholog, names(cont_table))]
  x$norm_length <- x$gene_length / x$number_pres
  return(x)
})

## Plot the number of HAVs per gene
number_table <- filter_snps[[1]][!duplicated(filter_snps[[1]]$human_ortholog),]

number_table <- table(number_table$number_pres[number_table$number_pres < 12])
number_table <- data.frame(HAV_numb = as.numeric(names(number_table)), value = number_table)

ggplot(number_table, aes(x = HAV_numb, y = value.Freq)) + geom_bar(stat = "identity") + theme_classic() + scale_x_continuous(breaks = number_table$HAV_numb)


## Load the results of the permutation

permute_out <- readRDS("Permutation_corrected_GO_analysis.rds")

permute_table <- permute_out[[1]][!duplicated(permute_out[[1]]$Term),]

ggplot(permute_table, aes(x = Pop.Hits, y = permute_result)) + geom_point() + theme_classic() + xlim(c(0,2000))

## Run for pc1
assayed.genes <- unique(ortholog_df$human_ortholog)
de.genes <- unique(filter_snps[[1]]$human_ortholog)

gene.vector <- as.integer(assayed.genes%in%de.genes)
names(gene.vector) <- assayed.genes
head(gene.vector)

## Make vector of gene lengths
gene_lengths <- ortholog_df$gene_length[match(names(gene.vector), ortholog_df$human_ortholog)]


## run the PWF

pwf <- nullp(gene.vector, bias.data = gene_lengths)

## Run go analysis with basic categories
GO.wall <- goseq(pwf, "hg19", "geneSymbol", test.cats = c("GO:CC", "GO:BP", "GO:MF", "KEGG"))
GO.wall$Bonf_pvalue <- GO.wall$over_represented_pvalue*nrow(GO.wall)
GO.wall$cat_term <- paste(GO.wall$category, GO.wall$term, sep = "~")
head(GO.wall)

GO.nobias <- goseq(pwf, "hg19", "geneSymbol", method="Hypergeometric", test.cats = c("GO:CC", "GO:BP", "GO:MF", "KEGG"))
GO.nobias$Bonf_pvalue <- GO.nobias$over_represented_pvalue*nrow(GO.nobias)
GO.nobias$cat_term <- paste(GO.nobias$category, GO.nobias$term, sep = "~")
head(GO.nobias)

GO.wall$nobiaspvalue <- GO.nobias$over_represented_pvalue[match(GO.wall$cat_term, GO.nobias$cat_term)]

# plot(GO.wall$over_represented_pvalue, GO.wall$nobiaspvalue)

GO.wall$permute <- permute_out[[1]]$permute_result[match(GO.wall$cat_term, permute_out[[1]]$Term)]

ggplot(GO.wall, aes(x = log(over_represented_pvalue)*-1, y = log(nobiaspvalue)*-1)) + geom_point()
ggplot(GO.wall, aes(x = over_represented_pvalue, y = permute)) + geom_point()
ggplot(GO.wall, aes(x = log(over_represented_pvalue - nobiaspvalue), y = permute)) + geom_point()

ggplot(GO.wall[!(is.na(GO.wall$permute)),], aes(x = (log(sub_df$nobiaspvalue)*-1) - (log(sub_df$over_represented_pvalue)*-1), y = numInCat)) + geom_point()


ggplot(GO.wall[!(is.na(GO.wall$permute)),], aes(x = log(over_represented_pvalue)*-1, y = log(nobiaspvalue)*-1, colour = permute)) + geom_point() + scale_colour_viridis()

ggplot(GO.wall[!(is.na(GO.wall$permute)),], aes(x = (log(sub_df$nobiaspvalue)*-1) - (log(sub_df$over_represented_pvalue)*-1), y = 1-permute, colour = 1-permute)) + geom_point() + scale_colour_viridis()

sub_df <- GO.wall[!(is.na(GO.wall$permute)),]
sub_df$pvalue_diff <- (log(sub_df$nobiaspvalue)*-1) - (log(sub_df$over_represented_pvalue)*-1)
summary(lm(permute ~ pvalue_diff, sub_df))



# Get the mapping from ENSEMBL 2 Entrez
en2eg <- as.list(org.Hs.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg <- as.list(org.Hs.egPATH)
# Define a function which gets all unique KEGG IDs
# associated with a set of Entrez IDs
grepKEGG <- function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
# Apply this function to every entry in the mapping from
# ENSEMBL 2 Entrez to combine the two maps
kegg <- lapply(en2eg,grepKEGG,eg2kegg)
head(kegg)


# pwf=nullp(genes,"hg19","ensGene")
KEGG <- goseq(pwf,gene2cat=kegg)
head(KEGG)

KEGG <- goseq(pwf,'hg19','geneSymbol',test.cats="KEGG")
enriched.KEGG <- KEGG[p.adjust(KEGG$over_represented_pvalue, method="BH")<.05,]


