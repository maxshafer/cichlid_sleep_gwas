#!/bin/bash

#SBATCH --job-name=VarFilt                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=168:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/step9_VarFiltstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/step9_VarFiltstderr.txt

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

#export your required environment variables below
#################################################


#add your command lines below
#############################

# comma separated df with rows, samples, interval
# file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_1000x.csv"

# this is the second column of index_array_40x.csv
# INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 2 -d ','`

# ~/gatk-4.2.4.1/gatk VariantFiltration -R ~/pro-ects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna -V ~/projects/cichlids_2/sra_reads_nobackup/cohort_geno_gathered.g.vcf.gz -O ~/projects/cichlids_2/sra_reads_nobackup/cohort_geno_gathered_filtered.g.vcf.gz

gatk VariantFiltration -R ~/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa -V ~/projects/cichlids_2/sra_reads_nobackup/cohort_db_geno.g.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "LQ_filter" -O ~/projects/cichlids_2/sra_reads_nobackup/cohort_db_geno_filtered.g.vcf.gz


