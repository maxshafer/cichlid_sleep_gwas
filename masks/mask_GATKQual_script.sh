#!/bin/bash

#SBATCH --job-name=Qualmask                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/Qualmaskstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/Qualmaskstderr.txt

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ayasha.abdallawyse@mail.utoronto.ca         #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load Java
module load StdEnv/2020
module load bedtools/2.30.0
module load gatk

#export your required environment variables below
#################################################


#add your command lines below
#############################

gatk VariantFiltration -R ~/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa -V ~/scratch/temp_data/NMBU_cohort_genotyped_whole.g.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "LQ_filter" -O ~/projects/def-mshafer/ayasha/cichlids_gwas/NMBU_cohort_genotyped_lq.filt.g.vcf.gz

gunzip -c ~/projects/def-mshafer/ayasha/cichlids_gwas/NMBU_cohort_genotyped_lq.filt.g.vcf.gz | awk '{if ($7=="LQ_filter") { print $1"\t"$2-1"\t"$2; }}' > ~/projects/def-mshafer/genome/masks/NMBU_cohort_genotyped.LQ.tmp.bed

mergeBed -i ~/projects/def-mshafer/genome/masks/NMBU_cohort_genotyped.LQ.tmp.bed > ~/projects/def-mshafer/genome/masks/NMBU_cohort_genotyped.LQ.bed

