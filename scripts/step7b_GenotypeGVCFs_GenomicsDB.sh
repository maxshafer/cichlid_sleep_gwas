#!/bin/bash

#SBATCH --job-name=GenoGVCF                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=6G              #This is the memory reserved per core.
#Total memory reserved: 6GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step7b_GenoGVCFstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step7b_GenoGVCFstderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-25%25
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

# comma separated df with rows, samples, interval
# file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_40x.csv"
#file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_80x.csv"
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_chr.csv"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 3 -d ','`

~/gatk-4.2.4.0/gatk GenotypeGVCFs -R ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -V gendb://$INTERVAL -O intervaled_genotyping_chr/cohort_${INTERVAL}.g.vcf.gz --tmp-dir $TMPDIR


# ~/gatk-4.2.4.1/gatk GenotypeGVCFs -R /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -V cohort.g.vcf.gz -L /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/intervals_ASM185804v2_chr/${INTERVAL} -O cohort_interval_vcfs_chr/cohort_geno_${INTERVAL}.g.vcf.gz

