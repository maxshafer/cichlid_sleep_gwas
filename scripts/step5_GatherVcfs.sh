#!/bin/bash

#SBATCH --job-name=GatherVcfs                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=24G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=6:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step5_GatherVcfs_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step5_GatherVcfs_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-123%123
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

# comma separated df with rows and samples
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_samples.csv"

# this is the second column of the index
SAMPLE=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 2 -d ','`
SAMPLE=${SAMPLE%.sra}

java -jar $EBROOTPICARD/picard.jar GatherVcfs -O ${SAMPLE}_variants.g.vcf.gz -I ${SAMPLE}_0000-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0001-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0002-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0003-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0004-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0005-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0006-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0007-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0008-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0009-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0010-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0011-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0012-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0013-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0014-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0015-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0016-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0017-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0018-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0019-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0020-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0021-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0022-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0023-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0024-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0025-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0026-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0027-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0028-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0029-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0030-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0031-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0032-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0033-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0034-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0035-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0036-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0037-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0038-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0039-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0040-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0041-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0042-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0043-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0044-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0045-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0046-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0047-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0048-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0049-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0050-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0051-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0052-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0053-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0054-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0055-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0056-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0057-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0058-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0059-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0060-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0061-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0062-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0063-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0064-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0065-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0066-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0067-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0068-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0069-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0070-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0071-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0072-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0073-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0074-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0075-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0076-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0077-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0078-scattered.interval_list_variants.g.vcf.gz -I ${SAMPLE}_0079-scattered.interval_list_variants.g.vcf.gz


~/gatk-4.2.4.1/gatk IndexFeatureFile -I ${SAMPLE}_variants.g.vcf.gz

