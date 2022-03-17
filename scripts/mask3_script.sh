#!/bin/bash

#SBATCH --job-name=Qualmask                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/Qualmaskstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/Qualmaskstderr.txt

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load BEDTools/2.30.0-GCC-10.3.0

#export your required environment variables below
#################################################


#add your command lines below
#############################

# gunzip -c cohort_genotyped_whole_filtered.g.vcf.gz | awk '{if ($7=="my_filter") { print $1"\t"$2-1"\t"$2; }}' > ../masks/cohort_geno_gathered.g.vcf.gz.LQ.tmp.bed

mergeBed -i ../masks/cohort_geno_gathered.g.vcf.gz.LQ.tmp.bed > ../masks/cohort_geno_gathered.g.vcf.gz.LQ.bed

