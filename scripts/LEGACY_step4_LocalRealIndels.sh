#!/bin/bash

#SBATCH --job-name=LocalRealIndels_184                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#Total memory reserved: 16GB

#SBATCH --time=08:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/LocalRealIndels_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/LocalRealIndels_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=57489,57490,57495,57497,57500%5
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

~/gatk-4.2.4.1/gatk RealignerTargetCreator -R /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna -I SRR96${SLURM_ARRAY_TASK_ID}_5_dedup.bam -o SRR96${SLURM_ARRAY_TASK_ID}_6_realign.intervals

~/gatk-4.2.4.1/gatk IndelRealigner -R /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna -targetIntervals SRR96${SLURM_ARRAY_TASK_ID}_6_realign.intervals -I SRR96${SLURM_ARRAY_TASK_ID}_5_dedup.bam -o SRR96${SLURM_ARRAY_TASK_ID}_6_realign.bam
