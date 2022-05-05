#!/bin/bash

#SBATCH --job-name=GWAS                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=128G              #This is the memory reserved per core.
#Total memory reserved: 128GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step13a_GWASstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step13a_GWASstderr.txt

#You selected an array of jobs from 1 to 25 with 25 simultaneous jobs
#SBATCH --array=1-7,9-25%25
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

## Final steps!! Run the GWAS and the PGLS-GWAS!
## These run per chromosome, on the allele frequencies files generated above
## They take awhile to run (especially the PGLS), and require lots of memory
## Run in the folder where you want the output to end up

# comma separated df with rows, samples, interval
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.chrs"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

Rscript ../../scripts/GWASrun.R ../alleleFreqs/speciesCodesAndGenomeIDs_for_GWAS_2_cohort_db_geno_${INTERVAL}_fromProbabilities_AF.txt ../../scripts/combined_diel_patterns_2022-04-26_dp_updated_names.csv peak_amplitude 

