

# Make command for Gathering scattered vcf files (output from interval'd HaplotypeCaller)

# intervals <- formatC(c(0:39), width = 4, flag = "0")

# intervals <- noquote(paste("-I SRR96${SLURM_ARRAY_TASK_ID}_", intervals, "-scattered.interval_list_variants.g.vcf.gz", sep = ""))

# write_clip(noquote(paste(intervals, collapse = " ")))


# java -jar $EBROOTPICARD/picard.jar GatherVcfs -I SRR96${SLURM_ARRAY_TASK_ID}_0000-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0001-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0002-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0003-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0004-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0005-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0006-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0007-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0008-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0009-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0010-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0011-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0012-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0013-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0014-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0015-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0016-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0017-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0018-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0019-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0020-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0021-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0022-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0023-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0024-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0025-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0026-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0027-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0028-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0029-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0030-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0031-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0032-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0033-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0034-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0035-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0036-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0037-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0038-scattered.interval_list_variants.g.vcf.gz -I SRR96${SLURM_ARRAY_TASK_ID}_0039-scattered.interval_list_variants.g.vcf.gz -O SRR96${SLURM_ARRAY_TASK_ID}_variants.g.vcf.gz

# Make command for Combining vcfs from the same sample after scatter/gather during Haplotypecalling (step4)

intervals <- formatC(c(0:79), width = 4, flag = "0")
intervals <- paste("-I ${SAMPLE}_", intervals, "-scattered.interval_list_variants.g.vcf.gz", sep = "")
intervals  <- c("java -jar $EBROOTPICARD/picard.jar GatherVcfs -O ../variant_vcf_files/${SAMPLE}_variants.g.vcf.gz", intervals)

print("Command for step4_Haplotype")
print(noquote(paste(intervals, collapse = " ")))

# Make command for Combining vcfs from different samples (for running Genotyping)

# I'm only doing 5 for now, so can do it manually, but can use the list.files() option on the folder like  before

samples <- read.csv("~/projects/cichlids_2/scripts/index_samples.csv", head = F, row.names = "V1")
samples <- substr(samples$V2, start = 1, stop = 10) 
samples <- noquote(paste("-V ", samples, "_variants.g.vcf.gz", sep = ""))
samples <- c("~/gatk-4.2.4.0/gatk GenomicsDBImport --genomicsdb-workspace-path ../$INTERVAL --tmp-dir $TMPDIR -L /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/intervals_ASM185804v2_80x/${INTERVAL}", samples)

print("Command for step6b_Genotyping")
print(noquote(paste(samples, collapse = " ")))

# ~/gatk-4.2.4.1/gatk CombineGVCFs -R /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -O cohort.g.vcf.gz -V SRR9657489_variants.g.vcf.gz -V SRR9657490_variants.g.vcf.gz -V SRR9657495_variants.g.vcf.gz -V SRR9657497_variants.g.vcf.gz -V SRR9657500_variants.g.vcf.gz

# Make command for Combining vcfs after Genotyping (step 7b)
# Files are in ~/projects/cichlids_2/sra_reads_nobackup/intervaled_genotyping_80x/

samples <- list.files("~/projects/cichlids_2/sra_reads_nobackup/intervaled_genotyping_80x/")
samples <- samples[!(grepl(".tbi", samples))]
samples <- noquote(paste("-I ", samples, sep = ""))
samples <- c("java -jar $EBROOTPICARD/picard.jar GatherVcfs -O ../cohort_db_geno.g.vcf.gz", samples)

print("command for step8b_GatherVcfs.sh")
print(noquote(paste(samples, collapse = " ")))

