#!/bin/bash

#SBATCH --job-name=snpEff                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=15G              #This is the memory reserved per core.
#Total memory reserved: 15GB

#SBATCH --time=6:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step14_snpEff_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step14_snpEff_stderr.txt

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
#############################

## OK for each chromosome, I need to annotate all SNPs using snpEff
## Will need to unzip, modify (".1"), annotate (genes), cut (first column), then rezip

# comma separated df with rows, samples, interval
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.chrs"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

# Unzip it
# gunzip -c cohort_db_geno_${INTERVAL}.hardfiltered_SNPS.biallelic.NoSingletons.g.vcf.gz > unzipped_${INTERVAL}_g.vcf

# Modify it

# cat unzipped_${INTERVAL}_g.vcf | sed "s/^${INTERVAL}/${INTERVAL}.1/" > unzipped_${INTERVAL}_mod.g.vcf

# rm unzipped_${INTERVAL}_g.vcf

# Annotate it
# java -Xmx8g -jar ~/bin/snpEff/snpEff.jar -v ASM185804v2 unzipped_${INTERVAL}_mod.g.vcf > unzipped_${INTERVAL}_mod.ann.g.vcf

# rm unzipped_${INTERVAL}_mod.g.vcf

# Extract it

# cat unzipped_${INTERVAL}_mod.ann.g.vcf | ~/bin/snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx8g -jar ~/bin/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].DISTANCE" > extracted_${INTERVAL}.txt
gzip extracted_${INTERVAL}.txt

# Cut it

# cat unzipped_${INTERVAL}_mod.ann.g.vcf | cut -f1 -d: > unzipped_${INTERVAL}_final.g.vcf

rm unzipped_${INTERVAL}_mod.ann.g.vcf

# Zip it

# gzip -c unzipped_${INTERVAL}_final.g.vcf > final_${INTERVAL}.g.vcf.gz

# rm unzipped_${INTERVAL}_final.g.vcf



