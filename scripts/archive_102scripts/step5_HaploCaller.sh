#!/bin/bash

#SBATCH --job-name=HaploCallerInf                   #This is the name of your job
#SBATCH --cpus-per-task=4                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=6G              #This is the memory reserved per core.
#Total memory reserved: 24GB

#SBATCH --time=167:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step5_HaploCallerG_%a_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/step5_HaploCallerG_%a_stderr.txt

#You selected an array of jobs from 1 to 9 with 9 simultaneous jobs
#SBATCH --array=57489,57490,57495,57497,57500,57501,57502,57505,57507,57508,57509,57511,57512,57518,57524,57525,57546,57557,57560,57561,57562,57563,57564,57565,57566,57568,65632,65633,65634,65635,65638,65640,65641,65643,65644,65645,65647,65648,65649,65659,65662,65663,65670,65671,65672,65673,65684,65687,65692,65703,65708,65717,73848,73849,73852,73859,73890,73893,73894,73895,73912,73929,73934,73956,73963,73968,73971,73972,73978,73981,73982,73983,73984,73987,74003,74012,74015,74024,74025,74027,74034,74035,74053,75315,75318,75319,75320,75329,75331,75332,75333,75334,75336,75337,75338,75361,75362,75373,75374,75375,75378,75393%102
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load GATK/3.7-0-Java-1.8.0_92

#export your required environment variables below
#################################################


#add your command lines below
#############################

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna -I SRR96${SLURM_ARRAY_TASK_ID}_6_realign.bam -nct 4 -ERC GVCF -minPruning 1 -minDanglingBranchLength 4 --max_alternate_alleles 10 --genotyping_mode DISCOVERY -o SRR96${SLURM_ARRAY_TASK_ID}_variants.g.vcf.gz

