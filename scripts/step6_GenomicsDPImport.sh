#!/bin/bash

#SBATCH --job-name=GenomicsDBImport                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=24G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=6:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/GenomicsDBImportout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/GenomicsDBImporterr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-80%80
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
# comma separated df with rows, samples, interval
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_array_80x.csv"

# this is the second column of index_array_40x.csv
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 3 -d ','`


~/gatk-4.2.4.1/gatk GenomicsDBImport --genomicsdb-workspace-path $INTERVAL --tmp-dir $TMPDIR --genomicsdb-max-alternate-alleles 7 -L /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/intervals_ASM185804v2_80x/${INTERVAL}  -V SRR9657475_variants.g.vcf.gz -V SRR9657476_variants.g.vcf.gz -V SRR9657489_variants.g.vcf.gz -V SRR9657490_variants.g.vcf.gz -V SRR9657492_variants.g.vcf.gz -V SRR9657495_variants.g.vcf.gz -V SRR9657497_variants.g.vcf.gz -V SRR9657500_variants.g.vcf.gz -V SRR9657501_variants.g.vcf.gz -V SRR9657502_variants.g.vcf.gz -V SRR9657505_variants.g.vcf.gz -V SRR9657507_variants.g.vcf.gz -V SRR9657508_variants.g.vcf.gz -V SRR9657509_variants.g.vcf.gz -V SRR9657511_variants.g.vcf.gz -V SRR9657512_variants.g.vcf.gz -V SRR9657518_variants.g.vcf.gz -V SRR9657524_variants.g.vcf.gz -V SRR9657525_variants.g.vcf.gz -V SRR9657526_variants.g.vcf.gz -V SRR9657528_variants.g.vcf.gz -V SRR9657533_variants.g.vcf.gz -V SRR9657536_variants.g.vcf.gz -V SRR9657546_variants.g.vcf.gz -V SRR9657557_variants.g.vcf.gz -V SRR9657560_variants.g.vcf.gz -V SRR9657561_variants.g.vcf.gz -V SRR9657562_variants.g.vcf.gz -V SRR9657563_variants.g.vcf.gz -V SRR9657564_variants.g.vcf.gz -V SRR9657565_variants.g.vcf.gz -V SRR9657566_variants.g.vcf.gz -V SRR9657567_variants.g.vcf.gz -V SRR9657568_variants.g.vcf.gz -V SRR9665632_variants.g.vcf.gz -V SRR9665633_variants.g.vcf.gz -V SRR9665634_variants.g.vcf.gz -V SRR9665635_variants.g.vcf.gz -V SRR9665638_variants.g.vcf.gz -V SRR9665640_variants.g.vcf.gz -V SRR9665641_variants.g.vcf.gz -V SRR9665642_variants.g.vcf.gz -V SRR9665643_variants.g.vcf.gz -V SRR9665644_variants.g.vcf.gz -V SRR9665645_variants.g.vcf.gz -V SRR9665646_variants.g.vcf.gz -V SRR9665647_variants.g.vcf.gz -V SRR9665648_variants.g.vcf.gz -V SRR9665649_variants.g.vcf.gz -V SRR9665659_variants.g.vcf.gz -V SRR9665662_variants.g.vcf.gz -V SRR9665663_variants.g.vcf.gz -V SRR9665670_variants.g.vcf.gz -V SRR9665671_variants.g.vcf.gz -V SRR9665672_variants.g.vcf.gz -V SRR9665673_variants.g.vcf.gz -V SRR9665684_variants.g.vcf.gz -V SRR9665687_variants.g.vcf.gz -V SRR9665692_variants.g.vcf.gz -V SRR9665703_variants.g.vcf.gz -V SRR9665708_variants.g.vcf.gz -V SRR9665717_variants.g.vcf.gz -V SRR9673838_variants.g.vcf.gz -V SRR9673839_variants.g.vcf.gz -V SRR9673845_variants.g.vcf.gz -V SRR9673848_variants.g.vcf.gz -V SRR9673849_variants.g.vcf.gz -V SRR9673852_variants.g.vcf.gz -V SRR9673859_variants.g.vcf.gz -V SRR9673890_variants.g.vcf.gz -V SRR9673893_variants.g.vcf.gz -V SRR9673894_variants.g.vcf.gz -V SRR9673895_variants.g.vcf.gz -V SRR9673912_variants.g.vcf.gz -V SRR9673914_variants.g.vcf.gz -V SRR9673915_variants.g.vcf.gz -V SRR9673925_variants.g.vcf.gz -V SRR9673926_variants.g.vcf.gz -V SRR9673927_variants.g.vcf.gz -V SRR9673929_variants.g.vcf.gz -V SRR9673934_variants.g.vcf.gz -V SRR9673956_variants.g.vcf.gz -V SRR9673963_variants.g.vcf.gz -V SRR9673967_variants.g.vcf.gz -V SRR9673968_variants.g.vcf.gz -V SRR9673969_variants.g.vcf.gz -V SRR9673971_variants.g.vcf.gz -V SRR9673972_variants.g.vcf.gz -V SRR9673978_variants.g.vcf.gz -V SRR9673981_variants.g.vcf.gz -V SRR9673982_variants.g.vcf.gz -V SRR9673983_variants.g.vcf.gz -V SRR9673984_variants.g.vcf.gz -V SRR9673987_variants.g.vcf.gz -V SRR9674000_variants.g.vcf.gz -V SRR9674003_variants.g.vcf.gz -V SRR9674012_variants.g.vcf.gz -V SRR9674015_variants.g.vcf.gz -V SRR9674024_variants.g.vcf.gz -V SRR9674025_variants.g.vcf.gz -V SRR9674027_variants.g.vcf.gz -V SRR9674034_variants.g.vcf.gz -V SRR9674035_variants.g.vcf.gz -V SRR9674053_variants.g.vcf.gz -V SRR9675315_variants.g.vcf.gz -V SRR9675318_variants.g.vcf.gz -V SRR9675319_variants.g.vcf.gz -V SRR9675320_variants.g.vcf.gz -V SRR9675329_variants.g.vcf.gz -V SRR9675331_variants.g.vcf.gz -V SRR9675332_variants.g.vcf.gz -V SRR9675333_variants.g.vcf.gz -V SRR9675334_variants.g.vcf.gz -V SRR9675336_variants.g.vcf.gz -V SRR9675337_variants.g.vcf.gz -V SRR9675338_variants.g.vcf.gz -V SRR9675373_variants.g.vcf.gz -V SRR9675374_variants.g.vcf.gz -V SRR9675375_variants.g.vcf.gz -V SRR9675377_variants.g.vcf.gz -V SRR9675378_variants.g.vcf.gz -V SRR9675387_variants.g.vcf.gz -V SRR9675393_variants.g.vcf.gz

