#!/bin/bash

#SBATCH --job-name=fastq_dump                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#SBATCH --tmp=128G
#Total memory reserved: 16GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/fastq_dump_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/fastq_dump_stderr.txt

#You selected an array of jobs from 1 to n with n simultaneous jobs
#SBATCH --array=1-123%123
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

# module load SRA-Toolkit/2.11.2-gompi-2021a
module load SRA-Toolkit/2.8.1-3-centos_linux64

#export your required environment variables below
#################################################


#add your command lines below
#############################

# comma separated df with rows and samples
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_samples.csv"

# this is the second column of the index
SAMPLE=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 2 -d ','`
SAMPLE=${SAMPLE%.sra}

fastq-dump --gzip --outdir $TMPDIR --split-files /scicore/home/schiera/gizevo30/ncbi/public/sra/${SAMPLE}.sra

cp $TMPDIR/${SAMPLE}_1.fastq.gz /scicore/home/schiera/gizevo30/projects/cichlids_2/sra_reads_nobackup/${SAMPLE}_1.fastq.gz
cp $TMPDIR/${SAMPLE}_2.fastq.gz /scicore/home/schiera/gizevo30/projects/cichlids_2/sra_reads_nobackup/${SAMPLE}_2.fastq.gz


