library(stringr)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(here)
library(tictoc)
library(ggrastr)

########################################################################
#######   Defining  arguments  #########################################
########################################################################

percentile <- 0.1

########################################################################
#######   Defining  arguments  #########################################
########################################################################

files <- list.files("sra_reads_nobackup/combined_ann/", pattern = "_pvals_ann.gz")

comparison <- c("pc1", "pc2", "total_rest")

########################################################################
#######   Load/merge/filter per comparison  ############################
########################################################################

per_chr <- lapply(comparison, function(comp) {
  
  gwas.datasets <- lapply(files, function(x) {
    df <- fread(paste("sra_reads_nobackup/combined_ann/", x, sep = ""), showProgress = T)
    columns <- c(c(2:3, as.numeric(grep(comp, colnames(df)))))
    df <- df[, ..columns] # Change last two numbers to modify which comparison to keep and filter by
    
    colnames(df) <- c("CHROM","POS","piVals","PGLSpiVals")
    
    df$location <- paste(df$CHROM, df$POS, sep = ":")
    df <- df[!duplicated(df$location),]
    
    
    return(df)
  })
  
  merged <- Reduce(rbind, gwas.datasets)
  
  pi <- unique(merged[,c("piVals", "location")]) %>% slice_min(order_by = piVals, prop = percentile) %>% pull(location)
  pgls <- unique(merged[,c("PGLSpiVals", "location")]) %>% slice_min(order_by = PGLSpiVals, prop = percentile) %>% pull(location)
    
  merged2 <- merged[merged$location %in% unique(pi, pgls),]
  
  
  merged2$label <- "non-significant"
  
  
  ## Load and append circadian genes and significant genes
  
  per_chr <- readRDS("circadian-genes_perchr.rds")
  filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))
  
  per_chr[[1]]$label <- "Circadian"
  filter_snps[[2]]$label <- "PC1_enriched"
  pc1 <- rbind(filter_snps[[2]], per_chr[[1]])
  
  pc1 <- pc1[,c("CHROM", "POS", "piVals", "PGLSpiVals", "location", "label")]
  pc1 <- pc1[!duplicated(pc1$location),]
  
  merged2 <- rbind(merged2, pc1)
  merged2$label <- factor(merged2$label, levels = c("PC1_enriched", "Circadian", "non-significant"))
  
  merged2$piVals <- log(merged2$piVals)*-1
  merged2$PGLSpiVals <- log(merged2$PGLSpiVals)*-1
  
  plot <- ggplot(merged2, aes(x = PGLSpiVals, y = piVals, colour = label)) + scale_color_manual(values = c("red", "yellow", "grey55")) + ggrastr::rasterise(geom_point(alpha = 0.55, size = 0.01, shape = 19)) + theme_classic() + ylab("GWAS-GLM p-value") + xlab("GWAS-pGLS p-value") + ylim(c(0,18)) + xlim(c(0,18))
  plot <- plot + theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"))
  plot <- plot + plot_layout(width = unit(c(50), "mm"), height = unit(c(50), "mm"))
  
  ggsave("scatter_plot_circadian.pdf", plot = plot, width = 10, height = 10, units = "in")
  
  
  tic()
  ggsave("scatter_plot_circadian.png", plot = plot, dpi = 500, width = 4, height = 4, units = "in")
  toc()
  
  pc1_line_dot <- ggplot(merged, aes(x = POS, y = log(piVals)*-1)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(), axis.title = element_blank()) + scale_colour_distiller(palette = "Blues", direction = 1) + facet_wrap(~CHROM, nrow = 1, scales = "free_x")
  pc1_line_dot_pgls <- ggplot(merged, aes(x = POS, y = log(PGLSpiVals)*1)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(), axis.title = element_blank()) + scale_colour_distiller(palette = "Reds", direction = -1) + facet_wrap(~CHROM, nrow = 1, scales = "free_x")
  
  plot2 <- pc1_line_dot + pc1_line_dot_pgls + plot_layout(ncol = 1, width = unit(c(300), "mm"), height = unit(c(50), "mm"))
  tic()
  ggsave("man_plot_test.png", plot = plot2, dpi = 500, width = 4, height = 4, units = "in")
  toc()
  
  return(merged)
  
})


saveRDS(per_chr, file = paste("sra_reads_nobackup/dump/combined_ann/filter_SNPs_perchr", percentile, "percentile", sep = "_"))







#### Below is for making manhattan plots


library(data.table)
library(dplyr)
library(readr)
library(patchwork)
library(here)
library(ggnewscale)
library(stringr)
library(ggtree)
library(tidyr)
library(viridis)


# Load the phylogeny

lt_phylo <- ape::read.nexus("scripts/05_BEAST_RAxML.tre")

# Load phenotype file for plotting and annotating with frequencies
phenotypes <- read.csv("pheno_data/GWAS_data_PCs_TR_2024-01-28.csv", row.names = "X")

## Load in the snps, and subset to unique gene/location matches
filter_snps <- readRDS(here("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile"))

snps <- lapply(filter_snps, function(x) {
  x <- x[!duplicated(x$location),]
  return(x)
})


facet_scales_chrom <- list(
  scale_x_continuous(limits = c(0,40673430)),
  scale_x_continuous(limits = c(0,36523203)),
  scale_x_continuous(limits = c(0,35549522)),
  scale_x_continuous(limits = c(0,39714817)),
  scale_x_continuous(limits = c(0,42433576)),
  scale_x_continuous(limits = c(0,64772279)),
  scale_x_continuous(limits = c(0,30527416)),
  scale_x_continuous(limits = c(0,35850837)),
  scale_x_continuous(limits = c(0,34704454)),
  scale_x_continuous(limits = c(0,39275952)),
  scale_x_continuous(limits = c(0,38600464)),
  scale_x_continuous(limits = c(0,34734273)),
  scale_x_continuous(limits = c(0,40509636)),
  scale_x_continuous(limits = c(0,39688505)),
  scale_x_continuous(limits = c(0,36041493)),
  scale_x_continuous(limits = c(0,38839487)),
  scale_x_continuous(limits = c(0,38636442)),
  scale_x_continuous(limits = c(0,30963196)),
  scale_x_continuous(limits = c(0,37140374)),
  scale_x_continuous(limits = c(0,39199643)),
  scale_x_continuous(limits = c(0,45655644)),
  scale_x_continuous(limits = c(0,87567345))
)

pc1_line_tile <- ggplot(snps[[2]], aes(x = POS, y = 1, fill = log(piVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,24)) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)
pc1_line_tile_pgls <- ggplot(snps[[2]], aes(x = POS, y = 1, fill = log(PGLSpiVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,24)) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)
pc2_line_tile <- ggplot(snps[[3]], aes(x = POS, y = 1, fill = log(piVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,24)) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)
pc2_line_tile_pgls <- ggplot(snps[[3]], aes(x = POS, y = 1, fill = log(PGLSpiVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,24)) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)
tr_line_tile <- ggplot(snps[[4]], aes(x = POS, y = 1, fill = log(piVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,24)) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)
tr_line_tile_pgls <- ggplot(snps[[4]], aes(x = POS, y = 1, fill = log(PGLSpiVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,24)) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)


pc1_line_tile + pc1_line_tile_pgls + pc2_line_tile + pc2_line_tile_pgls + tr_line_tile + tr_line_tile_pgls + plot_layout(ncol = 1, guides = "collect")




## Need to extract the CHROM + POS for the intersect genes, and then make another plot on top with

gene_ids <- lapply(filter_snps, function(x) {
  
  gof <- x$ANN_GENE
  gof <- unlist(strsplit(gof, "-"))
  gof <- unique(gof[gof %in% gtf$gene])
  
  return(gof)
  
})

intersect_genes <- Reduce(intersect, list(gene_ids[[2]], gene_ids[[3]], gene_ids[[4]]))

intersect_df <- data.frame(intersect_genes = intersect_genes)

intersect_df$CHROM <- unlist(apply(intersect_df, 1, function(x) {
  return(unique(filter_snps[[2]]$CHROM[grep(x[1], filter_snps[[2]]$ANN_GENE)])[1])
}))
intersect_df$POS <- unlist(apply(intersect_df, 1, function(x) {
  return(unique(filter_snps[[2]]$POS[grep(x[1], filter_snps[[2]]$ANN_GENE)])[1])
}))



label_data <- snps[[2]]
label_data$label <- intersect_df$intersect_genes[match(label_data$location, intersect_df$location)]

label_plot <- ggplot(data = label_data, aes(x=POS, label = label, y = 1, colour = label)) + geom_point() + geom_text_repel(angle = 90) + theme_minimal() + theme(axis.text = element_blank(), axis.title = element_blank()) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)

label_plot + pc1_line_tile + pc1_line_tile_pgls + pc2_line_tile + pc2_line_tile_pgls + tr_line_tile + tr_line_tile_pgls + plot_layout(ncol = 1, guides = "collect")









pc1_line_tile2 <- ggplot(snps[[2]], aes(x = POS, y = 1, fill = log(piVals*PGLSpiVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Blues", direction = 1) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)
pc2_line_tile2 <- ggplot(snps[[3]], aes(x = POS, y = 1, fill = log(piVals*PGLSpiVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Blues", direction = 1) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)
tr_line_tile2 <- ggplot(snps[[4]], aes(x = POS, y = 1, fill = log(piVals*PGLSpiVals)*-1)) + geom_tile(width = 100000) + theme_classic() + theme(axis.text = element_blank(), axis.title = element_blank()) + scale_fill_distiller(palette = "Blues", direction = 1) + facet_wrap(~CHROM, nrow = 1, scales = "free_x") + ggh4x::facetted_pos_scales(x = facet_scales_chrom)


pc1_line_tile2 + pc2_line_tile2 + tr_line_tile2 + plot_layout(ncol = 1, guides = "collect")


pc1_line_dot <- ggplot(snps[[2]], aes(x = POS, y = log(piVals)*-1)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(), axis.title = element_blank()) + scale_colour_distiller(palette = "Blues", direction = 1) + facet_wrap(~CHROM, nrow = 1, scales = "free_x")
pc1_line_dot_pgls <- ggplot(snps[[2]], aes(x = POS, y = log(PGLSpiVals)*1)) + geom_point() + theme_classic() + theme(axis.text.x = element_blank(), axis.title = element_blank()) + scale_colour_distiller(palette = "Reds", direction = -1) + facet_wrap(~CHROM, nrow = 1, scales = "free_x")


pc1_line_dot + pc1_line_dot_pgls + plot_layout(ncol = 1, guides = "collect")




