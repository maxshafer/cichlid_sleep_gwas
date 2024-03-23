#!/bin/bash

#SBATCH --job-name=7aGenoGVCF                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=168:00:00        #This is the time that your task will run
#SBATCH --qos=7days           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/GenoGVCFstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/GenoGVCFstderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
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
module load gatk

#export your required environment variables below
#################################################


#add your command lines below
############################

# # comma separated df with rows, samples, interval
# file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_80x.csv"

# # this is the second column of index_array_40x.csv
# INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 3 -d ','`

cd ~/scratch/temp_data/
gatk GenotypeGVCFs -R ~/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa -V NMBU_cohort_final.g.vcf.gz -O NMBU_cohort_genotyped_whole.g.vcf.gz

