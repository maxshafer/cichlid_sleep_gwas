#!/bin/bash

#SBATCH --job-name=snpEff                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=15G              #This is the memory reserved per core.
#Total memory reserved: 15GB

#SBATCH --time=6:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/step14_snpEff_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/step14_snpEff_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-22%22
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

## OK for each chromosome, I need to annotate all SNPs using snpEff
## Will need to unzip, modify (".1"), annotate (genes), cut (first column), then rezip

# comma separated df with rows, samples, interval
file_list="/home/ayasha/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.chromosomes"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

# Unzip it
gunzip -c /home/ayasha/projects/def-mshafer/gwas_output/NMBU_cohort_genotyped_${INTERVAL}.hardfiltered_SNPS.biallelic.NoSingletons.g.vcf.gz > /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_g.vcf

# Modify it

cat /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_g.vcf | sed "s/^${INTERVAL}/${INTERVAL}.1/" > /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_mod.g.vcf

rm /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_g.vcf

# Annotate it
java -Xmx8g -jar ~/snpEff/snpEff.jar -v GCA_001858045.3 /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_mod.g.vcf > /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_mod.ann.g.vcf

rm unzipped_${INTERVAL}_mod.g.vcf

# Extract it

cat /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_mod.ann.g.vcf | ~/snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx8g -jar ~/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].DISTANCE" > /home/ayasha/projects/def-mshafer/gwas_output/extracted_${INTERVAL}.txt
gzip /home/ayasha/projects/def-mshafer/gwas_output/extracted_${INTERVAL}.txt

# Cut it

cat /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_mod.ann.g.vcf | cut -f1 -d: > /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_final.g.vcf

rm /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_mod.ann.g.vcf

# Zip it

gzip -c /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_final.g.vcf > final_${INTERVAL}.g.vcf.gz

rm /home/ayasha/projects/def-mshafer/gwas_output/unzipped_${INTERVAL}_final.g.vcf




