library(data.table)
library(dplyr)
library(ggplot2)
library(Rgb)
library(gsheet)
library(patchwork)

## OK, this script is for taking the output of step15/step16 (merging and filtering after snpEff annotation)

## The annotated and filtered SNPs are first loaded in
## This is a list of 4 (Summary, pc1, pc2, and total rest)

filter_snps <- readRDS("sra_reads_nobackup/combined_ann/filter_SNPs_perchr_1e-05_percentile")

## Load in gtf
gtf <- read.gtf("genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf")
gtf <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]

########################################################################
#######   ID'ing circadian genes  ######################################
########################################################################

## Circ genes

url <- 'https://docs.google.com/spreadsheets/d/1R29FPj5cWffYq0sMud0R7M-eoyPtSaqi1DiN-VezcRo/edit#gid=0'
genes <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
circ_genes <- c(genes$Tilapia_Gene.Name.1, genes$Tilapia_Gene.Name_2)
circ_genes <- unique(circ_genes)
circ_genes_human <- tolower(genes$Human_Gene.Name)


cgg <- circ_genes[circ_genes %in% gtf$gene]
cghg <- circ_genes_human[circ_genes_human %in% gtf$gene]
clock_genes <- unique(c(cgg, cghg))


files <- list.files("sra_reads_nobackup/combined_ann/", pattern = "_pvals_ann.gz")

comparison <- c("pc1_60-species", "pc2_60-species", "total_rest_60-species")

### Pulling data from snps

dfs_all <- lapply(files, function(x) fread(paste("sra_reads_nobackup/combined_ann/", x, sep = ""), showProgress = T))
names(dfs_all) <- files

per_chr <- lapply(comparison, function(comp) {
  
  gwas.datasets <- lapply(files, function(x) {
    df <- dfs_all[[x]]
    # df <- fread(paste("sra_reads_nobackup/combined_ann/", x, sep = ""), showProgress = T)
    df <- df[df$`ANN[*].GENE` %in% clock_genes,]
    columns <- c(c(2:3,16:21, as.numeric(grep(comp, colnames(df)))))
    df <- df[, ..columns] # Change last two numbers to modify which comparison to keep and filter by
    
    if (comp == "summary") {
      colnames(df) <- c("CHROM","POS","REF","ALT","ANN_GENE","ANN_IMPACT","ANN_EFFECT","ANN_DISTANCE","summary_mean","summary_min", "summary_prod")
    } else {
      colnames(df) <- c("CHROM","POS","REF","ALT","ANN_GENE","ANN_IMPACT","ANN_EFFECT","ANN_DISTANCE","piVals","PGLSpiVals")
    }
    
    return(df)
  })
  
  merged <- Reduce(rbind, gwas.datasets)
  merged$location <- paste(merged$CHROM, merged$POS, sep = ":")
  
  pi <- merged %>% group_by(ANN_GENE) %>% slice_min(order_by = piVals, n = 1)
  pgls <- merged %>% group_by(ANN_GENE) %>% slice_min(order_by = PGLSpiVals, n = 1)
  
  merged <- merged[merged$location %in% c(pi$location, pgls$location),]
  
  return(merged)
  
})

saveRDS(per_chr, file = "circadian-genes_perchr_60sp.rds")







# 
# ########################################################################
# #######   IDing Human Orthologs  #######################################
# ########################################################################
# 
# # Load in gtf for removing non-coding genes
# gtf <- read.gtf("genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.gtf")
# gtf2 <- gtf[gtf$feature == "gene" & gtf$gene_biotype == "protein_coding",]
# 
# # gtf_other <- read.gtf("~/Downloads/ncbi_dataset_gtf/ncbi_dataset/data/GCF_001858045.2/genomic.gtf")
# 
# # Load in ensembl indices for IDing orthologs
# ensembl_btw <- read.csv("ensembl_data/ensembl-O_niloticus_UMD_NMBU-Human.txt", sep = "\t")
# ensembl_ncbi <- read.csv("ensembl_data/O_niloticus_UMD_NMBU-ensembl_ncbi-id.txt")
# 
# ## OK, each list of snps includes ncRNAs/lcnRNAs (and others), but also dashed pairs of gene names (because snp is in between genes)
# 
# gof <- filter_snps[[2]]$ANN_GENE
# gof <- unlist(strsplit(gof, "-"))
# gof <- unique(gof[gof %in% gtf2$gene])
# 
# gof_id <- gtf2$db_xref[match(gof, gtf2$gene)]
# gof_id <- substr(gof_id, 8, 30)
# 
# gof <- unique(c(gof, ncbi_gene_table$Symbol[match(gof_id, ncbi_gene_table$NCBI.GeneID)]))
# 
# gof_id <- paste("LOC", gof_id, sep = "")
# 
# ## Test which ones can be found
# 
# gof_df <- data.frame(gof = gof, 
#                      human_gene_name = ifelse(tolower(gof) %in% tolower(ensembl_btw$Human.gene.name), "yes", "no"),
#                      ensembl_gene_name = ifelse(gof %in% ensembl_ncbi$Gene.name,"yes", "no"),
#                      ensembl_gene_id = ifelse(gof %in% ensembl_ncbi$NCBI.gene..formerly.Entrezgene..accession, "yes", "no"))
# 
# xtabs(~ensembl_gene_name+ensembl_gene_id+human_gene_name, gof_df)
# 
# ## From this, there are 228 genes that can't be found in the ensembl lists, but they are all in the gtf
# ## In the gtf is the ncbi numeric id (under dbkey), which is another way to find orthologs and get ensembl IDs
# 
# ## I think I might have to go gof -> gtf$dbxref -> NCBI_gene_table$gene_name -> Ensembl
# 
# 
# ### OK, so I have the snpEFF annotated gene name/ID from ncbi
# ### I need to convert this to ensembl gene ID
# ### then use this 
# 
# 
# ncbi_gene_table <- read.csv("genome/NCBI_gene-table_Oreochromis-niloticus_UMD_NMBU.tsv", sep = "\t")
# 
# ncbi_gof_id <- gtf2$db_xref[match(gof, gtf2$gene)]
# ncbi_gof_id <- substr(ncbi_gof_id, 8, 30)
# 
# ncbi_gof <- ncbi_gene_table$Symbol[match(ncbi_gof_id, ncbi_gene_table$NCBI.GeneID)]
# 
# ncbi_gof_df <- data.frame(ncbi_gof = ncbi_gof, 
#                           ncbi_gof_id = ncbi_gof_id,
#                           human_gene_name = ifelse(tolower(ncbi_gof) %in% tolower(ensembl_btw$Human.gene.name), "yes", "no"),
#                           ensembl_gene_name = ifelse(ncbi_gof %in% ensembl_btw$Gene.name,"yes", "no"),
#                           ensembl_gene_id = ifelse(ncbi_gof %in% ensembl_ncbi$NCBI.gene..formerly.Entrezgene..accession, "yes", "no"),
#                           ensembl_gene_id_id = ifelse(paste("LOC", ncbi_gof_id, sep = "") %in% ensembl_ncbi$NCBI.gene..formerly.Entrezgene..accession, "yes", "no"))
# 
# xtabs(~ensembl_gene_name+ensembl_gene_id+human_gene_name+ensembl_gene_id_id, ncbi_gof_df)
# 
# ## OK, now only 187 that can't be found, including things like 'myo9aa' which in ENSEMBL is MYO9A (paralog not in zebrafish?), 
# ## that is orthologous to human Myo9a, but is not found in my cross-searches
# 
# ## I think the major problem is that both NCIB and Ensembl orthologies with humans are better, BUT, there is litte
# ## cross over between databases - eg poor quality genes have orthologs with human, but are not ID'd across platforms
# ## as being the same gene
# 
# ## OK, can I revisit the scripts for orthology from ncbi??
# 

# 
# 
# 
# 
# ## OK, query the gtf2 with my gene list to find GeneIDs and see if they are in the ensembl db
# 
# gof <- unique(filter_snps[[2]]$ANN_GENE)
# gof <- unlist(strsplit(gof, "-"))
# idof <- gtf2$GeneID[match(gof, gtf2$gene)]
# ensembl_idof <- ensembl_ncbi$Gene.stable.ID[match(gof, ensembl_ncbi$NCBI.gene..formerly.Entrezgene..ID)]
# 
# # table(idof %in% ensembl_ncbi$NCBI.gene..formerly.Entrezgene..ID)
# 
# hidof <- ensembl_btw$Human.gene.name[match(ensembl_idof, ensembl_btw$Gene.stable.ID)]
# 
# ## Do it by gene name, not id?
# 
# ensembl_ncbi <- read.csv("ensembl_data/O_niloticus_UMD_NMBU-ensembl_ncbi-accession.txt")
# 
# 
# gof <- unique(filter_snps[[2]]$ANN_GENE)
# gof <- unlist(strsplit(gof, "-"))
# 
# ensembl_idof <- ensembl_ncbi$Gene.name[match(gof, ensembl_ncbi$NCBI.gene..formerly.Entrezgene..accession)]
# 
# hidof2 <- ensembl_btw$Human.gene.name[match(ensembl_idof, ensembl_btw$Gene.name)]
# 
# ## Combining the above gets me to 390 genes (out of 700-900), and a bit more than just giving the niloticus gene names (325)
# length(unique(c(hidof, hidof2)))
# 
# write_clip(unique(c(hidof, hidof2)))
# 
# 
# ## Ugh, when there is a dash, it is because snpeff put that there - the snp is between two genes, I can split and use both for enrichment
# ## I could also just use the closest one?
# 
# 
# pc1_genes <- unlist(strsplit(unique(filter_snps[[2]]$ANN_GENE), "-"))
# pc2_genes <- unlist(strsplit(unique(filter_snps[[3]]$ANN_GENE), "-"))
# tr_genes <- unlist(strsplit(unique(filter_snps[[4]]$ANN_GENE), "-"))
# 
# pc1_loc <- unlist(strsplit(unique(filter_snps[[2]]$location), "-"))
# pc2_loc <- unlist(strsplit(unique(filter_snps[[3]]$location), "-"))
# tr_loc <- unlist(strsplit(unique(filter_snps[[4]]$location), "-"))
# 
# 
# table(pc1_genes %in% pc2_genes)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
