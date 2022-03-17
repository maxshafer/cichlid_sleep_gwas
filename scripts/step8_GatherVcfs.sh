#!/bin/bash

#SBATCH --job-name=GatherVcfs                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=24G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step8_GatherVcfs_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step8_GatherVcfs_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
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
# module load Java
module load picard/2.26.2-Java-1.8
# module load BCFtools/1.12-GCC-10.3.0

#export your required environment variables below
#################################################


#add your command lines below
#############################

java -jar $EBROOTPICARD/picard.jar GatherVcfs --OUTPUT cohort_geno_gathered.g.vcf.gz -I cohort_geno_0000-scattered.interval_list.g.vcf.gz -I cohort_geno_0001-scattered.interval_list.g.vcf.gz -I cohort_geno_0002-scattered.interval_list.g.vcf.gz -I cohort_geno_0003-scattered.interval_list.g.vcf.gz -I cohort_geno_0004-scattered.interval_list.g.vcf.gz -I cohort_geno_0005-scattered.interval_list.g.vcf.gz -I cohort_geno_0006-scattered.interval_list.g.vcf.gz -I cohort_geno_0007-scattered.interval_list.g.vcf.gz -I cohort_geno_0008-scattered.interval_list.g.vcf.gz -I cohort_geno_0009-scattered.interval_list.g.vcf.gz -I cohort_geno_0010-scattered.interval_list.g.vcf.gz -I cohort_geno_0011-scattered.interval_list.g.vcf.gz -I cohort_geno_0012-scattered.interval_list.g.vcf.gz -I cohort_geno_0013-scattered.interval_list.g.vcf.gz -I cohort_geno_0014-scattered.interval_list.g.vcf.gz -I cohort_geno_0015-scattered.interval_list.g.vcf.gz -I cohort_geno_0016-scattered.interval_list.g.vcf.gz -I cohort_geno_0017-scattered.interval_list.g.vcf.gz -I cohort_geno_0018-scattered.interval_list.g.vcf.gz -I cohort_geno_0019-scattered.interval_list.g.vcf.gz -I cohort_geno_0020-scattered.interval_list.g.vcf.gz -I cohort_geno_0021-scattered.interval_list.g.vcf.gz -I cohort_geno_0022-scattered.interval_list.g.vcf.gz -I cohort_geno_0023-scattered.interval_list.g.vcf.gz -I cohort_geno_0024-scattered.interval_list.g.vcf.gz

# bcftools sort cohort_geno_gathered.g.vcf.gz -T $TMPDIR -Oz -o cohort_geno_gathered_sorted.g.vcf.g

# ~/gatk-4.2.4.1/gatk IndexFeatureFile -I cohort_geno_gathered.g.vcf.gz

