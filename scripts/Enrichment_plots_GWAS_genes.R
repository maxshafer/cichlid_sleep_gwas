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

source("scripts/variants_functions.R")

theme_set(theme_minimal())
theme_update(axis.text = element_text(colour = "black"))
## OK this script is for making plots from GO and pheno data/outputs
## This includes the DAVID analysis, but also the outputs from the pheno website https://evol.nhri.org.tw/phenome2/

filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

########################################################################
#######   Plot scatter + circadian genes  ##############################
########################################################################

per_chr <- readRDS("circadian-genes_perchr.rds")

theme_set(theme_classic())
theme_update(axis.text = element_text(colour = "black"))

### Solution is to probably make a flattened version for all dots (somehow), and add these on top?
### Or can I make some kind of density diagram? Probably too computationally intense?

per_chr[[1]]$label <- "circadian"
filter_snps[[2]]$label <- "PC1_enriched"
pc1 <- rbind(filter_snps[[2]], per_chr[[1]])

per_chr[[2]]$label <- "circadian"
filter_snps[[3]]$label <- "PC2_enriched"
pc2 <- rbind(filter_snps[[3]], per_chr[[2]])

per_chr[[3]]$label <- "circadian"
filter_snps[[4]]$label <- "TR_enriched"
tr <- rbind(filter_snps[[4]], per_chr[[3]])


pc1$text <- NA
pc1$text[1396] <- "crcp - calcitonin gene-related peptide-receptor component protein"
pc1_scatter <- ggplot(pc1, aes(x = log(piVals)*-1, y = log(PGLSpiVals)*-1, colour = label, label = text)) + geom_point(size = 0.05) + geom_text_repel() + scale_color_manual(values = c("red", "black")) + xlab("GWAS (log(p-value)*1)") + ylab("PGLS-GWAS (log(p-value)*-1)")

pc2_scatter <- ggplot(pc2, aes(x = log(piVals)*-1, y = log(PGLSpiVals)*-1, colour = label)) + geom_point(size = 0.05) + theme_classic() + scale_color_manual(values = c("red", "black")) + xlab("GWAS (log(p-value)*1)") + ylab("PGLS-GWAS (log(p-value)*-1)")
tr_scatter <- ggplot(tr, aes(x = log(piVals)*-1, y = log(PGLSpiVals)*-1, colour = label)) + geom_point(size = 0.05) + theme_classic() + scale_color_manual(values = c("red", "black")) + xlab("GWAS (log(p-value)*1)") + ylab("PGLS-GWAS (log(p-value)*-1)")

plot2 <- pc1_scatter + ylim(c(0,18)) + xlim(c(0,18)) + theme_classic() + ylab("GWAS-GLM p-value") + xlab("GWAS-pGLS p-value")
plot2 <- plot2 + plot_layout(width = unit(c(25), "mm"), height = unit(c(25), "mm"))

tic()
ggsave("myplot_2.pdf", plot = plot2, width = 4, height = 4, units = "in")
toc()


genes <- readRDS(file = "orthos/orthologs_perchr_1e-05_percentile.rds")
genes_new <- genes$gene_names

annoCharts.2 <- readRDS("DAVID-GO_analysis_perchr_1e-05_percentile.rds")

# Reshape for plotting
measure.vars <- matrix(c("PValue", "Bonferroni", "Benjamini", "FDR", "Count", "Fold.Enrichment", "X.", "PValue", "Bonferroni", "Benjamini", "FDR", "Count", "FE", "X.", "PValue.value", "Bonf.value", "Benj.value", "FDR.value", "Counts.value", "FoldE.value", "X.value"), nrow=7, ncol = 3)
melted <- apply(measure.vars, 1, function(x) reshape2::melt(annoCharts.2, id.vars = c("Category", "Term"), measure.vars = x[1], variable.name = x[2], value.name = x[3]))


go_analysis <- melted[[1]]
go_analysis$Bonf.value <- melted[[2]]$Bonf.value
go_analysis$Benj.value <- melted[[3]]$Benj.value
go_analysis$FDR.value <- melted[[4]]$FDR.value
go_analysis$Count <- melted[[5]]$Counts.value
go_analysis$FoldE <- log(melted[[6]]$FoldE.value)
go_analysis$Perc <- melted[[7]]$X.value
# go_analysis$FoldE[go_analysis$L2 == "oreochromis"] <- go_analysis$FoldE[go_analysis$L2 == "oreochromis"]*-1
go_analysis$L3 <- paste(go_analysis$L1, go_analysis$L2, sep = "_")

## make new Term column without ID #s

go_analysis$Term_label <- go_analysis$Term

go_analysis$Term_label <- ifelse(grepl("~",go_analysis$Term_label), str_split_i(go_analysis$Term_label, "~", 2), go_analysis$Term_label)
go_analysis$Term_label <- ifelse(grepl(":",go_analysis$Term_label), str_split_i(go_analysis$Term_label, ":", 2), go_analysis$Term_label)

# Cluster with function, and plot
# GOTERM_BP_DIRECT
# UP_KEYWORDS
# GOTERM_MF_DIRECT
# GOTERM_CC_DIRECT
# KEGG_PATHWAY
# UP_SEQ_FEATURE
# REACTOME_PATHWAY
# BBID
# GAD_DISEASE
# DISGENET


go_analysis$L1 <- factor(go_analysis$L1, levels = names(genes_new))

go_analysis$L1 <- toupper(go_analysis$L1)
go_analysis$l1 <- str_replace(go_analysis$L1, "TR", "Total rest")

# Can I order the terms by the average pvalue across all three behaviours?
# 0.05, 0.025, 0.0005 for 75 and 200 heights

tissue.data <- subOrder(GO = go_analysis, pvalue = 0.01, species = "human", category = "UP_TISSUE", type = "min") 
tissues <- ggplot(data = tissue.data, aes(x = L1, y = Term_label, size = PValue.value, fill = FoldE, shape = L2)) + geom_point(alpha = 1, shape = 21, colour = "black") + scale_size_continuous(range = c(4, 0.5), limits = c(0, 0.05))
tissues <- tissues + scale_y_discrete(labels = label_wrap(50)) + scale_fill_viridis(limits = c(0,3)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 8, lineheight = 0.65), axis.title.x = element_blank(), axis.title.y = element_blank())

disgenet.data <- subOrder(GO = go_analysis, pvalue = 0.01, species = "human", category = "DISGENET", type = "min")
disgenet <- ggplot(data = disgenet.data, aes(L1, Term_label, size = PValue.value, fill = FoldE, shape = L2)) + geom_point(alpha = 1, shape = 21, colour = "black") + scale_size_continuous(range = c(4, 0.5), limits = c(0, 0.05))
disgenet <- disgenet + scale_y_discrete(labels = label_wrap(50)) + scale_fill_viridis(limits = c(0,3)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 8, lineheight = 0.65), axis.title.x = element_blank(), axis.title.y = element_blank())

mgi.data <- subOrder(GO = go_analysis, pvalue = 0.0001, species = "human", category = "MGI_Phenotype", type = "min")
mgi <- ggplot(data = mgi.data, aes(L1, Term_label, size = PValue.value, fill = FoldE, shape = L2)) + geom_point(alpha = 1, shape = 21, colour = "black") + scale_size_continuous(range = c(4, 0.5), limits = c(0, 0.05))
mgi <- mgi + scale_y_discrete(labels = label_wrap(50)) + scale_fill_viridis(limits = c(0,3)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 8, lineheight = 0.65), axis.title.x = element_blank(), axis.title.y = element_blank()) 

design1 <- "
AC
BC"

GO_plot <- tissues + disgenet + mgi + plot_layout(guides = "collect", design = design1, width = unit(c(12,12), "mm"), height = unit(c(30, 60), "mm"))


pdf("outs/third_figure_test.pdf", width = 10, height = 10)
GO_plot
dev.off()


design <- "
AADF
AAEF
BCEF"


pc1_scatter + pc2_scatter + tr_scatter + tissues + disgenet + mgi + plot_layout(guides = "collect", design = design, width = unit(c(56,56,25,25), "mm"), height = unit(c(56, 56,56), "mm")) # & theme(legend.position='bottom')



pdf("outs/second_figure_test.pdf", width = 5, height = 5)
pc1_scatter + plot_layout(width = unit(c(50), "mm"), height = unit(c(50), "mm"))
dev.off()

########################################################################
#######   Make Venn diagrams  ##########################################
########################################################################



## To map back from an enriched category
## write a function to do this?



plot(plotVennOverlap(category = "C0036341~Schizophrenia"))

plot(plotVennOverlap(category = "C0004352~Autistic Disorder"))

plot(plotVennOverlap(category = "Brain"))



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

names(gene_ids) <- c("summary", "pc1", "pc2", "tr")

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



doFisherExact(set1 = gene_ids$pc1, set2= gene_ids$pc2)
doFisherExact(set1 = gene_ids$pc1, set2= gene_ids$tr)
doFisherExact(set1 = gene_ids$pc2, set2= gene_ids$tr)

snps <- lapply(filter_snps, function(x) unique(x$location))
names(snps) <- c("summary", "pc1", "pc2", "tr")

doFisherExact(set1 = snps$pc1, set2= snps$tr, snps = TRUE)

