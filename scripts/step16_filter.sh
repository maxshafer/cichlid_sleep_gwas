#!/bin/bash

#SBATCH --job-name=filterRun                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=128G              #This is the memory reserved per core.
#Total memory reserved: 128GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/step16_filterRun_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/step16_filterRun_stderr.txt

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ayasha.abdallawyse@mail.utoronto.ca        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load StdEnv/2020
module load r/4.0.2

#export your required environment variables below
#################################################


#add your command lines below
#############################

## Final steps!! Run the GWAS and the PGLS-GWAS!
## These run per chromosome, on the allele frequencies files generated above
## They take awhile to run (especially the PGLS), and require lots of memory
## Run in the folder where you want the output to end up

Rscript /home/ayasha/projects/def-mshafer/cichlid_sleep_gwas/scripts/filterRun.R 0.00001
Rscript /home/ayasha/projects/def-mshafer/cichlid_sleep_gwas/scripts/filterRun.R 0.00005

