#!/bin/bash

#SBATCH --job-name=GenoGVCF                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=168:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/GenoGVCFstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/GenoGVCFstderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
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
############################

# # comma separated df with rows, samples, interval
# file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_80x.csv"

# # this is the second column of index_array_40x.csv
# INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 3 -d ','`

~/gatk-4.2.4.1/gatk GenotypeGVCFs -R /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -V cohort.g.vcf.gz -O cohort_genotyped_whole.g.vcf.gz

