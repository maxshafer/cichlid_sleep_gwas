#!/bin/bash

#SBATCH --job-name=CombineGVCF                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=128G              #This is the memory reserved per core.
#Total memory reserved: 128GB

#SBATCH --time=168:00:00        #This is the time that your task will run
#SBATCH --qos=168hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/CombineGVCFout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/CombineGVCFerr.txt
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ayasha.abdallawyse@mail.utoronto.ca        #You will be notified via email when your task ends or fails

#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load StdEnv/2023
module load java/21.0.1
module load gatk/4.4.0.0 

#export your required environment variables below
#################################################


#add your command lines below
#############################

# This code combines the vcfs in the old way, and is used with step7a for Genotyping
# This takes just under 6 days to run, but also needs to have the updated sample list

cd /home/ayasha/scratch/temp_data/

gatk --java-options "-Xmx8G" CombineGVCFs -R ~/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa -O NMBU_cohort_final.g.vcf.gz -V genotyped_variant_vcfs/SRR9665659_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665692_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665646_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665642_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665633_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665632_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675393_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674053_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657511_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657512_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674003_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674024_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675373_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675378_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673838_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673839_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673914_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673915_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657567_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657492_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674034_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674027_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665717_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665641_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657557_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657562_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674012_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674015_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657560_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657563_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665672_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665673_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673926_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673927_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665644_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665645_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657507_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657508_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665703_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665708_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657490_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657489_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657495_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657497_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657565_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657566_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665640_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665643_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657502_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657501_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673956_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673963_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673978_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673983_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665662_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665663_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665687_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665684_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673934_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673929_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657524_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657525_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673912_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675318_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675315_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673890_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673893_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675334_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675333_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675320_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675319_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657505_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657500_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675329_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675336_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675332_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675331_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665638_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673859_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673852_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675338_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675337_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665648_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665647_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665649_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673987_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675374_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9675375_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673849_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673848_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674035_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673982_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673845_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673925_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673895_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673894_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657533_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657536_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657564_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657561_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657528_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657479_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657480_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673971_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673972_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673984_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673981_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673968_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9674025_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673959_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673964_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673967_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9673969_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657518_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657509_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665682_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9665685_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657568_variants.g.vcf.gz -V genotyped_variant_vcfs/SRR9657546_variants.g.vcf.gz


