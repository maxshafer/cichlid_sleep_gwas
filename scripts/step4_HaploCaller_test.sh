#!/bin/bash

#SBATCH --job-name=HaploCallerInf                   #This is the name of your job
#SBATCH --cpus-per-task=1                 #This is the number of cores reserved
#SBATCH --mem-per-cpu=1G              #This is the memory reserved per core.
#SBATCH --tmp=6G
#Total memory reserved:1GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step4/step4_HaploCallerG_%a_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step4/step4_HaploCallerG_%a_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-1676%1676
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

# module load GATK/3.7-0-Java-1.8.0_92
module load Java

#export your required environment variables below
#################################################

#  #SBATCH --tmp=128G

#add your command lines below
#############################

# list of jobs that failed
failed_file="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/failed_ids.csv"
NUMB=`sed -n "$SLURM_ARRAY_TASK_ID"p "${failed_file}" | cut -f 2 -d ','`

# comma separated df with rows, samples, interval
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_80x.csv"

# This is the second column of index_array_**x.csv
SAMPLE=`sed -n "$NUMB"p "${file_list}" | cut -f 2 -d ','`
SAMPLE=${SAMPLE%.sra}
# this is the third column of index_array_**x.csv
INTERVAL=`sed -n "$NUMB"p "${file_list}" | cut -f 3 -d ','`

~/gatk-4.2.4.0/gatk HaplotypeCaller -R /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -L /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/intervals_ASM185804v2_80x/${INTERVAL} --max-alternate-alleles 5 -I bam_files/${SAMPLE}_5_dedup.bam -ERC GVCF -O intervaled_vcfs/${SAMPLE}_${INTERVAL}_variants.g.vcf.gz

