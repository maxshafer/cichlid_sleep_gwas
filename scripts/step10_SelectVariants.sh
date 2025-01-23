#!/bin/bash

#SBATCH --job-name=VarFilt                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/step10_VarFiltstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/step10_VarFiltstderr.txt

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ayasha.abdallawyse@mail.utoronto.ca        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load java
module load bedtools/2.30.0
module load StdEnv/2020  
module load gcc/9.3.0
module load bcftools/1.16
module load gatk

#export your required environment variables below
#################################################


#add your command lines below
#############################

# First, concatonate the bed files from all the masks and index it
# cd /scicore/home/schiera/gizevo30/projects/cichlids_2/masks/

# cat cohort_db_geno.DPMask.bed cohort_db_geno.LQmask.bed mapability_mask_100_90.bed | sortBed -faidx /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.chrs -i stdin | mergeBed -i stdin > cohort_db_geno.allMasksMerged.bed

# ~/gatk-4.2.4.0/gatk IndexFeatureFile -I cohort_db_geno.allMasksMerged.bed

# Then make a new Filtered g.vcf file, before hard filtereing using the mask
# cd /home/ayasha/projects/def-mshafer/ayasha/sra/

# gatk VariantFiltration -R ~/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa -V cohort_db_geno_filtered.g.vcf.gz --mask ../masks/cohort_db_geno.allMasksMerged.bed --mask-name "combined_filter" -O cohort_db_geno_SoftFilter.g.vcf.gz

# gatk SelectVariants -R ~/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa -V cohort_db_geno_SoftFilter.g.vcf.gz --exclude-filtered true -O cohort_db_geno_hardfiltered.g.vcf.gz


# This is subsetting the VCF to keep only balletic SNPs and excluding singletons and indels
# These numbers correspond to the total haplotypes present in the dataset (119*2 = 238) 
# bcftools view --max-ac 236 --min-ac 5 -m2 -M2 --exclude-types indels -O z -o /home/ayasha/projects/def-mshafer/gwas_output/maf5/NMBU_cohort_genotyped.hardfiltered_SNPS.biallelic.NoSingletons.g.vcf.gz /home/ayasha/projects/def-mshafer/gwas_output/NMBU_cohort_genotyped_hardfiltered.g.vcf.gz
bcftools view --max-ac 236 --min-ac 5 -m2 -M2 --exclude-types indels -O z -o ~/scratch/temp_data/output.vcf.gz NMBU_cohort_genotyped_hardfiltered.g.vcf.gz

gatk IndexFeatureFile -I ~/scratch/temp_data/output.vcf.gz

