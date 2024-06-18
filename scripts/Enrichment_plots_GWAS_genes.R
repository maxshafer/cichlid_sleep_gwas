library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(org.Hs.eg.db)
library(venneuler)
library(here)
library(patchwork)
library(viridis)
library(Rgb)
library(scales)
library(ggrepel)

source("scripts/variants_functions.R")

theme_set(theme_minimal())
theme_update(axis.text = element_text(colour = "black"))
## OK this script is for making plots from GO and pheno data/outputs
## This includes the DAVID analysis, but also the outputs from the pheno website https://evol.nhri.org.tw/phenome2/

ortholog_df <- readRDS("GTF_human_orthologs_all.rds")

filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

########################################################################
#######   GO Plots  ####################################################
########################################################################
## This script is for loading in the snp2go outputs, and redoing the multiple testing correction
## Will need to read in individual and merge by pc1/pc2/tr

# This is the full list of GO terms tested (for doing the multiple testing correction)
# Also the list including term id matching to term names (and description)
total_go_list <- read.csv("mart_export_human_go_terms_oreochromis_orthos_all.txt", sep = "\t", header = F)
go_names_df <- read.csv("mart_export_human_go_terms.txt", sep = "\t")

go_out <- lapply(c("pc1", "pc2", "tr"), function(x) {
  file_list <- list.files("snp2go_outs/", pattern = x, full.names = T)
  go <- Reduce(rbind, lapply(file_list, function(x) readRDS(x)))
  go$Bonf.pvalue <- go$P*nrow(total_go_list)
  go$FDR_BH <- p.adjust(go$P, n = nrow(total_go_list))
  go$pheno <- as.character(x)
  go$Category <- vapply(strsplit(go$GO,":"), `[`, 1, FUN.VALUE=character(1))
  go$Term <- vapply(strsplit(go$GO,":"), `[`, 2, FUN.VALUE=character(1))
  go$Term[grep("HPO:", go$GO)] <- vapply(strsplit(go$GO[grep("HPO:", go$GO)],":"), `[`, 3, FUN.VALUE=character(1))
  go$Term[grep("MGD:", go$GO)] <- vapply(strsplit(go$GO[grep("MGD:", go$GO)],":"), `[`, 3, FUN.VALUE=character(1))
  go$Term[grep("MGD:", go$GO)] <- vapply(strsplit(go$GO[grep("MGD:", go$GO)],":"), `[`, 3, FUN.VALUE=character(1))
  go$Term[grep("GO:", go$GO)] <- go_names_df$GO.term.name[match(go$GO[grep("GO:", go$GO)], go_names_df$GO.term.accession)]
  
  # HPO and MGD need to be the third
  return(go)
})

go_out <- Reduce(rbind, go_out)
go_out <- go_out[go_out$FDR_BH < 0.05,]

# Remove GO terms with children
go_out <- go_out[go_out$child.GOs =="",]

## Make plots

go_data <- go_out[grepl("GO:", go_out$GO) & go_out$FDR_BH < 0.000001,]
go_data_unique <- length(unique(go_data$Term))
go <- ggplot(data = go_data, aes(y = pheno, x = Term, fill = log(Bonf.pvalue)*-1)) + geom_tile(colour = "white") + theme_classic() + scale_fill_viridis(direction = 1)
go <- go + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

disease_data <- go_out[grepl("GENOMICS_ENGLAND:", go_out$GO) | grepl("PSYGENET:", go_out$GO) | grepl("CTD:", go_out$GO) & go_out$FDR_BH < 0.0001,]
disease_data_unique <- length(unique(disease_data$Term))
disease <- ggplot(data = disease_data, aes(y = pheno, x = Term, fill = log(Bonf.pvalue)*-1)) + geom_tile(colour = "white") + theme_classic() + scale_fill_viridis(direction = 1)
disease <- disease + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

design1 <- "AB"

GO_plot <- go + disease + plot_layout(guides = "collect", design = design1, width = unit(c(107.5,65), "mm"), height = unit(c(10,10), "mm"))


pdf("outs/third_figure_test.pdf", width = 10, height = 10)
GO_plot
dev.off()


########################################################################
#######   Make Venn diagrams  ##########################################
########################################################################

#### Make venn diagram for SNPs and genes in Oreochromis, and test for enrichment?

## Load gene table to convert to newest ids
ncbi_gene_table <- read.csv(here("genome","NCBI_gene-table_Oreochromis-niloticus_UMD_NMBU.tsv"), sep = "\t")

## Load in gtf
gtf <- read.gtf("genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf")
gtf <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]

filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

## Extract gene names
gene_ids <- lapply(filter_snps, function(x) {
  
  gof <- x$ANN_GENE
  gof <- unlist(strsplit(gof, "-"))
  gof <- unique(gof[gof %in% gtf$gene])
  
  
  ncbi_gof_id <- gtf$db_xref[match(gof, gtf$gene)]
  ncbi_gof_id <- substr(ncbi_gof_id, 8, 30)
  
  ncbi_gof <- ncbi_gene_table$Symbol[match(ncbi_gof_id, ncbi_gene_table$NCBI.GeneID)]
  
  return(ncbi_gof)
  
})

names(gene_ids) <- c("pc1_old", "pc2_old", "tr_old", "pc1", "pc2", "tr")

venn.data <- c(pc1_genes = length(gene_ids$pc1),
               pc2_genes = length(gene_ids$pc2),
               tr_genes = length(gene_ids$tr),
               'pc1_genes&pc2_genes' = length(intersect(gene_ids$pc1, gene_ids$pc2)),
               'pc1_genes&tr_genes' = length(intersect(gene_ids$pc1, gene_ids$tr)),
               'pc2_genes&tr_genes' = length(intersect(gene_ids$pc2, gene_ids$tr)),
               'pc1_genes&pc2_genes&tr_genes' = length(Reduce(intersect, list(gene_ids$pc1, gene_ids$pc2, gene_ids$tr))))

print(venn.data)

vd1 <- venneuler(venn.data)
plot(vd1)

## Test for enrichment in pairwise overlaps (all are significantly enriched)
doFisherExact(set1 = gene_ids$pc1, set2= gene_ids$pc2)
doFisherExact(set1 = gene_ids$pc1, set2= gene_ids$tr)
doFisherExact(set1 = gene_ids$pc2, set2= gene_ids$tr)





