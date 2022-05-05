#!/bin/bash

#SBATCH --job-name=StoF_BWAMEM_MergeBam                   #This is the name of your job
#SBATCH --cpus-per-task=7                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=8G              #This is the memory reserved per core.
#Total memory reserved: 56GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/StoF_BWAMEM_MergeBam_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/StoF_BWAMEM_MergeBam_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=1-119%119
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
# module load bwakit/0.7.15_x64-linux
module load BWA/0.7.17-goolf-1.7.20

#export your required environment variables below
#################################################


#add your command lines below
#############################

# comma separated df with rows and samples
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/index_samples.csv"

# this is the second column of the index
SAMPLE=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 2 -d ','`
SAMPLE=${SAMPLE%.sra}

set -o pipefail

java -jar $EBROOTPICARD/picard.jar SamToFastq --INPUT ${SAMPLE}_2_marked.bam --FASTQ /dev/stdout --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true --TMP_DIR /scicore/home/schiera/gizevo30/mstemp | bwa mem -M -t 7 -p /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna /dev/stdin | java -jar $EBROOTPICARD/picard.jar MergeBamAlignment --ALIGNED_BAM /dev/stdin --UNMAPPED_BAM ${SAMPLE}_2_marked.bam --OUTPUT ${SAMPLE}_3_piped.bam --REFERENCE_SEQUENCE /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna --CREATE_INDEX true --ADD_MATE_CIGAR true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --TMP_DIR /scicore/home/schiera/gizevo30/mstemp

rm ${SAMPLE}_2_marked.bam
