#!/bin/bash

#SBATCH --job-name=uBAM_MarkAdapt                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=24G              #This is the memory reserved per core.
#SBATCH --tmp=128G
#Total memory reserved: 24GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/uBAM_MarkAdapt_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/uBAM_MarkAdapt_stderr.txt

#You selected an array of jobs with 102 simultaneous jobs
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

# Step 1 

java -jar $EBROOTPICARD/picard.jar FastqToSam --FASTQ ${SAMPLE}_1.fastq.gz --FASTQ2 ${SAMPLE}_2.fastq.gz --OUTPUT ${SAMPLE}_1_temp.bam --SAMPLE_NAME ${SAMPLE} --LIBRARY_NAME ${SAMPLE} --PLATFORM illumina --TMP_DIR $TMPDIR


# rm ${SAMPLE}_1.fastq.gz
# rm ${SAMPLE}_2.fastq.gz

# Step 2

java -jar $EBROOTPICARD/picard.jar MarkIlluminaAdapters --INPUT ${SAMPLE}_1_temp.bam --OUTPUT ${SAMPLE}_2_marked.bam --METRICS ${SAMPLE}_2_markilluminaadapters_metrics.txt --TMP_DIR $TMPDIR

# rm ${SAMPLE}_1_temp.bam

