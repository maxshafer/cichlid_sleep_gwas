#!/bin/bash

#SBATCH --job-name=chrSub                   #This is the name of your job
#SBATCH --cpus-per-task=1                 #This is the number of cores reserved
#SBATCH --mem-per-cpu=6G              #This is the memory reserved per core.
#Total memory reserved: 6GB

#SBATCH --time=6:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/step11_%a_chrSub_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/step11_%a_chrSub_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-25%25
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ayasha.abdallawyse@mail.utoronto.ca         #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

# module load GATK/3.7-0-Java-1.8.0_92
module load java
module load StdEnv/2020  
module load gcc/9.3.0
module load bcftools/1.13
module load gatk

#export your required environment variables below
#################################################


#add your command lines below
#############################

# The below should be done as an array job
# This is subsetting a whole-genome VCF to get per-chromosome VCFs
# This queries the list of chromosomes, with the unplaced one last

# comma separated df with rows, samples, interval
file_list="/home/ayasha/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.chromosomes"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

bcftools view -r $INTERVAL --threads 8 -O z -o /home/ayasha/projects/def-mshafer/gwas_output/NMBU_cohort_genotyped_${INTERVAL}.hardfiltered_SNPS.biallelic.NoSingletons.g.vcf.gz /home/ayasha/projects/def-mshafer/gwas_output/NMBU_cohort_genotyped_SNPS.biallelic.NoSingletons.g.vcf.gz
# bcftools view -r $INTERVAL --threads 8 -O z -o cohort_db_geno_${INTERVAL}.g.vcf.gz cohort_db_geno.g.vcf.gz

gatk IndexFeatureFile -I /home/ayasha/projects/def-mshafer/gwas_output/NMBU_cohort_genotyped_${INTERVAL}.hardfiltered_SNPS.biallelic.NoSingletons.g.vcf.gz
# ~/gatk-4.2.4.0/gatk IndexFeatureFile -I cohort_db_geno_${INTERVAL}.g.vcf.gz

