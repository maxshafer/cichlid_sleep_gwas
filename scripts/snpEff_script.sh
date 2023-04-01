#!/bin/bash

#SBATCH --job-name=Qualmask                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=120:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/snpeff_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/snpeff_stderr.txt

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

#export your required environment variables below
#################################################


#add your command lines below
#############################

java -jar ~/snpEff/snpEff.jar -v Orenil1.0.75 cohort_genotyped_whole_hardfiltered_SNPS.biallelic.NoSingletons.g.vcf.gz > cohort_genotyped_whole_hardfiltered_SNPS.biallelic.NoSingletons.anno.g.vcf.gz



