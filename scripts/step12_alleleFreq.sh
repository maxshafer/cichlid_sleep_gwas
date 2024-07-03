#!/bin/bash

#SBATCH --job-name=alleleFreq                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/step12_alleleFreqstdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/step12_alleleFreqstderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-25%25
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
module load gcc/10.2.0 
module load boost/1.76

#export your required environment variables below
#################################################


#add your command lines below
#############################

## OK next steps are to estimate allel frequencies using a Hardy-weinberg equilibrium equation
## taking into account the quality of the genome call and the depth (I think). This uses Milan's
## custom scripts/programs in the 'evo' package on github

# comma separated df with rows, samples, interval
file_list="/home/ayasha/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.chromosomes"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

~/evo/Build/evo alleleFreq -n NMBU_cohort_genotyped_${INTERVAL}_fromProbabilities --use-genotype-probabilities /home/ayasha/projects/def-mshafer/gwas_output/NMBU_cohort_genotyped_${INTERVAL}.hardfiltered_SNPS.biallelic.NoSingletons.g.vcf.gz /home/ayasha/projects/def-mshafer/cichlid_sleep_gwas/scripts/speciesCodesAndGenomeIDs_for_GWAS_2.txt


# ~/evo/Build/evo alleleFreq -n test_fromProbabilities --use-genotype-probabilities test.g.vcf.gz medians_for_GWAS_withoutfilter_speciesCodesAndGenomeIDs_wOretan.txt


