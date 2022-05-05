#!/bin/bash

#SBATCH --job-name=step8b_GatherVcfs                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=128G              #This is the memory reserved per core.
#SBATCH --tmp=210G
#Total memory reserved: 128GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step8b_GatherVcfs_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step8b_GatherVcfs_stderr.txt

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

#export your required environment variables below
#################################################


#add your command lines below
#############################

java -jar $EBROOTPICARD/picard.jar GatherVcfs -O ../cohort_db_geno.g.vcf.gz -I cohort_0000-scattered.interval_list.g.vcf.gz -I cohort_0001-scattered.interval_list.g.vcf.gz -I cohort_0002-scattered.interval_list.g.vcf.gz -I cohort_0003-scattered.interval_list.g.vcf.gz -I cohort_0004-scattered.interval_list.g.vcf.gz -I cohort_0005-scattered.interval_list.g.vcf.gz -I cohort_0006-scattered.interval_list.g.vcf.gz -I cohort_0007-scattered.interval_list.g.vcf.gz -I cohort_0008-scattered.interval_list.g.vcf.gz -I cohort_0009-scattered.interval_list.g.vcf.gz -I cohort_0010-scattered.interval_list.g.vcf.gz -I cohort_0011-scattered.interval_list.g.vcf.gz -I cohort_0012-scattered.interval_list.g.vcf.gz -I cohort_0013-scattered.interval_list.g.vcf.gz -I cohort_0014-scattered.interval_list.g.vcf.gz -I cohort_0015-scattered.interval_list.g.vcf.gz -I cohort_0016-scattered.interval_list.g.vcf.gz -I cohort_0017-scattered.interval_list.g.vcf.gz -I cohort_0018-scattered.interval_list.g.vcf.gz -I cohort_0019-scattered.interval_list.g.vcf.gz -I cohort_0020-scattered.interval_list.g.vcf.gz -I cohort_0021-scattered.interval_list.g.vcf.gz -I cohort_0022-scattered.interval_list.g.vcf.gz -I cohort_0023-scattered.interval_list.g.vcf.gz -I cohort_0024-scattered.interval_list.g.vcf.gz

~/gatk-4.2.4.0/gatk IndexFeatureFile -I ../cohort_db_geno.g.vcf.gz
