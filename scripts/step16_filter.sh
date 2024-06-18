#!/bin/bash

#SBATCH --job-name=filterRun                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=128G              #This is the memory reserved per core.
#Total memory reserved: 128GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step16_filterRun_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step16_filterRun_stderr.txt

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load R/4.0.3-foss-2018b

#export your required environment variables below
#################################################


#add your command lines below
#############################

## Final steps!! This filters the SNPs based on p-values from gwas and pgls-gwas
## It also annotates them and separates by comparison

# Rscript /scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/filterRun.R 0.00001
# Rscript /scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/filterRun.R 0.00005

# Run the permutations
Rscript /scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/filterRun_permute.R 800 1000

