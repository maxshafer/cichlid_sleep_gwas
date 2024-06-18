library(stringr)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(here)
library(tictoc)
library(ggrastr)
library(readr)
library(stringr)
library(tidyr)
library(viridis)

########################################################################
#######   Defining  arguments  #########################################
########################################################################

files <- list.files("sra_reads_nobackup/combined_ann/", pattern = "_pvals_ann.gz")

########################################################################
#######   Load/merge/filter per comparison  ############################
########################################################################

dfs_all <- lapply(files, function(x) fread(paste("sra_reads_nobackup/combined_ann/", x, sep = ""), showProgress = T))
names(dfs_all) <- files

dfs_all <- lapply(dfs_all, function(x) {
  x$location <- paste(x$CHROM, x$POS, sep = ":")
  x <- x[!duplicated(x$location),]
  x <- x[, c("CHROM", "POS", "location", "piVals_pc1_60-species", "piVals_pc2_60-species", "piVals_total_rest_60-species", "PGLSpiVals_pc1_60-species", "PGLSpiVals_pc2_60-species", "PGLSpiVals_total_rest_60-species")]
  colnames(x) <- c("CHROM", "POS", "location", "piVals_pc1_60species", "piVals_pc2_60species", "piVals_total_rest_60species", "PGLSpiVals_pc1_60species", "PGLSpiVals_pc2_60species", "PGLSpiVals_total_rest_60species")
  return(x)
})

# This is the background
merged <- Reduce(rbind, dfs_all)

rm(dfs_all)

rm(man_plot)
pc1_scatter <- ggplot(merged, aes(x = log(piVals_pc1_60species)*-1, y = log(PGLSpiVals_pc1_60species)*-1)) + ggrastr::rasterise(geom_point(alpha = 1, size = 8, shape = 19), scale = 0.1) 
pc1_scatter <- pc1_scatter + theme_classic() + ylab("pGLS-GWAS p-value") + xlab("GLM-GWAS p-value") + ylim(c(0,22)) + xlim(c(0,22))
pc1_scatter <- pc1_scatter + theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"))
pc1_scatter <- pc1_scatter + plot_layout(width = unit(c(50), "mm"), height = unit(c(50), "mm"))

ggsave("outs/scatter_plot_pc1.png", plot = pc1_scatter, dpi = 500, width = 10, height = 10, units = "in")
rm(pc1_scatter)

pc2_scatter <- ggplot(merged, aes(x = log(piVals_pc2_60species)*-1, y = log(PGLSpiVals_pc2_60species)*-1)) + ggrastr::rasterise(geom_point(alpha = 1, size = 8, shape = 19), scale = 0.1) 
pc2_scatter <- pc2_scatter + theme_classic() + ylab("pGLS-GWAS p-value") + xlab("GLM-GWAS p-value") + ylim(c(0,22)) + xlim(c(0,22))
pc2_scatter <- pc2_scatter + theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"))
pc2_scatter <- pc2_scatter + plot_layout(width = unit(c(50), "mm"), height = unit(c(50), "mm"))

ggsave("outs/scatter_plot_pc2.png", plot = pc2_scatter, dpi = 500, width = 10, height = 10, units = "in")
rm(pc2_scatter)

tr_scatter <- ggplot(merged, aes(x = log(piVals_total_rest_60species)*-1, y = log(PGLSpiVals_total_rest_60species)*-1)) + ggrastr::rasterise(geom_point(alpha = 1, size = 8, shape = 19), scale = 0.1) 
tr_scatter <- tr_scatter + theme_classic() + ylab("pGLS-GWAS p-value") + xlab("GLM-GWAS p-value") + ylim(c(0,22)) + xlim(c(0,22))
tr_scatter <- tr_scatter + theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"))
tr_scatter <- tr_scatter + plot_layout(width = unit(c(50), "mm"), height = unit(c(50), "mm"))

ggsave("outs/scatter_plot_tr.png", plot = tr_scatter, dpi = 500, width = 10, height = 10, units = "in")
rm(tr_scatter)






#### Below is for making manhattan plots

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

########################################################################
########################  PC1  #########################################
########################################################################

man_plot <- ggplot(merged, aes(x = POS, y = log(piVals_pc1_60species)*-1, colour = CHROM)) + theme_classic() + facet_grid(~CHROM, space = "free_x", scales = "free_x") + theme(panel.spacing.x = unit(0.5, "mm"), strip.text.x = element_blank()) + ggrastr::rasterise(geom_point(alpha = 1, size = 6), scale = 0.1) 
man_plot <- man_plot + theme(axis.title = element_text(size = 10), axis.text.y = element_text(size = 8, colour = "black"), axis.text.x = element_blank(), legend.position = "none")  + ylim(c(0,22))
man_plot <- man_plot + ylab("GLM-GWAS p-value") + xlab("Chromosomal position") + plot_layout(width = unit(c(150), "mm"), height = unit(c(25), "mm"))

ggsave("outs/manhattan_pc1_glmGWAS.png", plot = man_plot, dpi = 500, width = 10, height = 3, units = "in")


man_plot <- ggplot(merged, aes(x = POS, y = log(PGLSpiVals_pc1_60species)*-1, colour = CHROM)) + theme_classic() + facet_grid(~CHROM, space = "free_x", scales = "free_x") + theme(panel.spacing.x = unit(0.5, "mm"), strip.text.x = element_blank()) + ggrastr::rasterise(geom_point(alpha = 1, size = 6), scale = 0.1) 
man_plot <- man_plot + theme(axis.title = element_text(size = 10), axis.text.y = element_text(size = 8, colour = "black"), axis.text.x = element_blank(), legend.position = "none")  + ylim(c(0,22))
man_plot <- man_plot + ylab("GLM-GWAS p-value") + xlab("Chromosomal position") + plot_layout(width = unit(c(150), "mm"), height = unit(c(25), "mm"))

ggsave("outs/manhattan_pc1_pglsGWAS.png", plot = man_plot, dpi = 500, width = 10, height = 3, units = "in")


########################################################################
########################  PC2  #########################################
########################################################################

man_plot <- ggplot(merged, aes(x = POS, y = log(piVals_pc2_60species)*-1, colour = CHROM)) + theme_classic() + facet_grid(~CHROM, space = "free_x", scales = "free_x") + theme(panel.spacing.x = unit(0.5, "mm"), strip.text.x = element_blank()) + ggrastr::rasterise(geom_point(alpha = 1, size = 6), scale = 0.1) 
man_plot <- man_plot + theme(axis.title = element_text(size = 10), axis.text.y = element_text(size = 8, colour = "black"), axis.text.x = element_blank(), legend.position = "none")  + ylim(c(0,22))
man_plot <- man_plot + ylab("GLM-GWAS p-value") + xlab("Chromosomal position") + plot_layout(width = unit(c(150), "mm"), height = unit(c(25), "mm"))

ggsave("outs/manhattan_pc2_glmGWAS.png", plot = man_plot, dpi = 500, width = 10, height = 3, units = "in")


man_plot <- ggplot(merged, aes(x = POS, y = log(PGLSpiVals_pc2_60species)*-1, colour = CHROM)) + theme_classic() + facet_grid(~CHROM, space = "free_x", scales = "free_x") + theme(panel.spacing.x = unit(0.5, "mm"), strip.text.x = element_blank()) + ggrastr::rasterise(geom_point(alpha = 1, size = 6), scale = 0.1) 
man_plot <- man_plot + theme(axis.title = element_text(size = 10), axis.text.y = element_text(size = 8, colour = "black"), axis.text.x = element_blank(), legend.position = "none")  + ylim(c(0,22))
man_plot <- man_plot + ylab("GLM-GWAS p-value") + xlab("Chromosomal position") + plot_layout(width = unit(c(150), "mm"), height = unit(c(25), "mm"))

ggsave("outs/manhattan_pc2_pglsGWAS.png", plot = man_plot, dpi = 500, width = 10, height = 3, units = "in")

########################################################################
########################  TR  ##########################################
########################################################################

man_plot <- ggplot(merged, aes(x = POS, y = log(piVals_total_rest_60species)*-1, colour = CHROM)) rm(merged+ theme_classic() + facet_grid(~CHROM, space = "free_x", scales = "free_x") + theme(panel.spacing.x = unit(0.5, "mm"), strip.text.x = element_blank()) + ggrastr::rasterise(geom_point(alpha = 1, size = 6), scale = 0.1) 
man_plot <- man_plot + theme(axis.title = element_text(size = 10), axis.text.y = element_text(size = 8, colour = "black"), axis.text.x = element_blank(), legend.position = "none")  + ylim(c(0,22))
man_plot <- man_plot + ylab("GLM-GWAS p-value") + xlab("Chromosomal position") + plot_layout(width = unit(c(150), "mm"), height = unit(c(25), "mm"))

ggsave("outs/manhattan_total_rest_glmGWAS.png", plot = man_plot, dpi = 500, width = 10, height = 3, units = "in")


man_plot <- ggplot(merged, aes(x = POS, y = log(PGLSpiVals_total_rest_60species)*-1, colour = CHROM)) + theme_classic() + facet_grid(~CHROM, space = "free_x", scales = "free_x") + theme(panel.spacing.x = unit(0.5, "mm"), strip.text.x = element_blank()) + ggrastr::rasterise(geom_point(alpha = 1, size = 6), scale = 0.1) 
man_plot <- man_plot + theme(axis.title = element_text(size = 10), axis.text.y = element_text(size = 8, colour = "black"), axis.text.x = element_blank(), legend.position = "none")  + ylim(c(0,22))
man_plot <- man_plot + ylab("GLM-GWAS p-value") + xlab("Chromosomal position") + plot_layout(width = unit(c(150), "mm"), height = unit(c(25), "mm"))

ggsave("outs/manhattan_total_rest_pglsGWAS.png", plot = man_plot, dpi = 500, width = 10, height = 3, units = "in")