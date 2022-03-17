#!/bin/bash

#SBATCH --job-name=VarFilt                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step10_VarFiltstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step10_VarFiltstderr.txt

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load Java
module load BEDTools/2.30.0-GCC-10.3.0

#export your required environment variables below
#################################################


#add your command lines below
#############################

# First, concatonate the bed files from all the masks and index it
# cd /scicore/home/schiera/gizevo30/projects/cichlids_2/masks/

# cat cohort_geno_gathered.DPmask.bed cohort_geno_gathered.g.vcf.gz.LQ.bed mapability_mask_100_90.bed | sortBed -faidx /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.chrs -i stdin | mergeBed -i stdin > cohort_geno_gathered.allMasksMerged.bed

# ~/gatk-4.2.4.1/gatk IndexFeatureFile -I cohort_geno_gathered.allMasksMerged.bed

# # Then make a new Filtered g.vcf file, before hard filtereing using the mask
cd /scicore/home/schiera/gizevo30/projects/cichlids_2/sra_reads_nobackup/

~/gatk-4.2.4.1/gatk VariantFiltration -R ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -V cohort_genotyped_whole.g.vcf.gz --mask ../masks/cohort_geno_gathered.allMasksMerged.bed --mask-name "combined_filter" -O temp_cohort_genotyped_whole_filtered.g.vcf.gz

~/gatk-4.2.4.1/gatk SelectVariants -R ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -V temp_cohort_genotyped_whole_filtered.g.vcf.gz --exclude-filtered true -O cohort_genotyped_whole_hardfiltered.g.vcf.gz


