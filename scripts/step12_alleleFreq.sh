#!/bin/bash

#SBATCH --job-name=alleleFreq                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step12_alleleFreqstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step12_alleleFreqstderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-23%23
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load Boost/1.76.0-GCC-10.3.0

#export your required environment variables below
#################################################


#add your command lines below
#############################

## OK next steps are to estimate allel frequencies using a Hardy-weinberg equilibrium equation
## taking into account the quality of the genome call and the depth (I think). This uses Milan's
## custom scripts/programs in the 'evo' package on github

# comma separated df with rows, samples, interval
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.chromosomes"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

~/evo/Build/evo alleleFreq -n NMBU_cohort_genotyped_whole_${INTERVAL}_fromProbabilities --use-genotype-probabilities /scicore/home/schiera/gizevo30/projects/cichlids_2/sra_reads_nobackup/dump/final_snps_chr/NMBU_cohort_genotyped_whole_${INTERVAL}_HardFilter_BiAllelic_NoSingletons.g.vcf.gz /scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/speciesCodesAndGenomeIDs_for_GWAS_2.txt


