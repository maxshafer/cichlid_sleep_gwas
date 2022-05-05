#!/bin/bash

#SBATCH --job-name=vcf_filt                  #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/vcf_filtstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/vcf_filtstderr.txt

#SBATCH --mail-type=COMPLETE,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load VCFtools/0.1.16-GCC-10.3.0

#export your required environment variables below
#################################################


#add your command lines below
#############################

vcftools --gzvcf cohort_genotyped_whole_hardfiltered.g.vcf.gz --max-missing 0.1 --recode --out cohort_genotyped_whole_hardfiltered_0.1.g.vcf.gz
vcftools --gzvcf cohort_genotyped_whole_hardfiltered.g.vcf.gz --max-missing 0.25 --recode --out cohort_genotyped_whole_hardfiltered_0.25.g.vcf.gz
vcftools --gzvcf cohort_genotyped_whole_hardfiltered.g.vcf.gz --max-missing 0.5 --recode --out cohort_genotyped_whole_hardfiltered_0.5.g.vcf.gz
