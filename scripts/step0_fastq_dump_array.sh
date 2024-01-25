#!/bin/bash

#SBATCH --job-name=fastq_dump                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#SBATCH --tmp=128G
#Total memory reserved: 16GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=$HOME/scratch/logs/fastq_dump_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=$HOME/scratch/logs/fastq_dump_stderr.txt

#You selected an array of jobs from 1 to n with n simultaneous jobs
#SBATCH --array=1-119%119
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ayasha.abdallawyse@mail.utoronto.ca        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load StdEnv
module load gcc
module load sra-toolkit

#export your required environment variables below
#################################################


#add your command lines below
#############################

# comma separated df with rows and samples
file_list="$HOME/cichlid_sleep_gwas/scripts/index_samples.csv"

# this is the second column of the index
SAMPLE=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 2 -d ','`
SAMPLE=${SAMPLE%.sra}

fastq-dump --gzip --outdir $TMPDIR --split-files $HOME/projects/def-mshafer/SRA_files/${SAMPLE}.sra

cp $TMPDIR/${SAMPLE}_1.fastq.gz $HOME/projects/def-mshafer/SRA_files/SRA_reads_nobackup/dump/${SAMPLE}_1.fastq.gz
cp $TMPDIR/${SAMPLE}_2.fastq.gz $HOME/projects/def-mshafer/SRA_files/SRA_reads_nobackup/dump/${SAMPLE}_2.fastq.gz


