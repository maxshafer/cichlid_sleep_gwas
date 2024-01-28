#!/bin/bash

#SBATCH --job-name=mask2_DP                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/masks2DPtdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/masks2DPtderr.txt

#SBATCH --mail-type=END,FAIL,TIME_LIMIT,COMPLETE
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load BEDTools/2.30.0-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0

#add your command lines below
#############################

## Getting (overall) depth from a VCF:
## This takes a random subset of nucleotide positions, and outputs the read depth. Use this to ID a good range of depths to use (look for the normal distribution)
# So the full thing runs, but takes a couple of hours in an srun - but looks good, similar to what he sees I guess, would choose something like 900/1000 - 1800/1900

# bcftools query -f '%INFO/DP\n' /scicore/home/schiera/gizevo30/projects/cichlids_2/sra_reads_nobackup/dump/NMBU_cohort_genotyped_whole_filtered.g.vcf.gz | awk '{if (rand() < 0.1) { if ($1 == ".") {$1 = 0;} print;}}' | gzip -c > NMBU_cohort_DP.temp.gz

# Run the command to make the mask using the cutoffs deteremined using the above

min=900
max=1900

bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' /scicore/home/schiera/gizevo30/projects/cichlids_2/sra_reads_nobackup/dump/NMBU_cohort_genotyped_whole_filtered.g.vcf.gz | awk -v mi="$min" -v ma="$max" '{if ($3 == ".") {$3 = 0;} if (($3 <= mi) || ($3 >= ma)) {print $1"\t"$2-1"\t"$2;}}' > NMBU_cohort_DPmask.tmp.bed

mergeBed -i NMBU_cohort_DPmask.tmp.bed > NMBU_cohort_DPmask.bed
