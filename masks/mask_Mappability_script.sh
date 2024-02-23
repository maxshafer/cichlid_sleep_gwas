#!/bin/bash

#SBATCH --job-name=Mapmask                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/home/ayasha/scratch/logs/Mappabilitymask_new_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/home/ayasha/scratch/logs/Mappabilitymask_new_stderr.txt

#SBATCH --array=1-23%23
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=ayasha.abdallawyse@mail.utoronto.ca        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

#export your required environment variables below
#################################################


#add your command lines below
#############################


# This is an alternative to step6a, and uses the newer GenomicsDBImport function to generate a DB for each interval
# Each interval is then Genotyped (with GenotypeGVCF), then combined with GatherGVCF before filtering

# comma separated df with rows, samples, interval
file_list="/home/ayasha/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.chromosomes"

# this is the second column of index_array
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

# Subset the whole genome for only the interval (e.g. chromosome) and place it in the temp drive for manipulation
# Will need to load in samtools, subset, then purge modules before loading bwa
module load samtools

samtools faidx /home/ayasha/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa $INTERVAL > /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta.fna

module purge

module load bwa

bwa index /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta.fna


# Run
# Use the github compiled program (or the old binary from blogpost '~/seqbility-20091110/splitfa')
/home/ayasha/seqbility-20091110/splitfa.c /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta.fna 100 > /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1

bwa aln -R 1000000 -O 3 -E 3 /home/ayasha/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1 > /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1.sai

# samse needs also the index and the input reads, as well as the output of aln
bwa samse /home/ayasha/projects/def-mshafer/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1.sai /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1 > /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1_aln-se.sam

# generate the raw and non raw mask
/home/ayasha/seqbility-20091110/gen_raw_mask.pl /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1_aln-se.sam > /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_rawMask_l100_s1.fa

/home/ayasha/seqbility-20091110/gen_mask.c -l 100 -r 0.9 /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_rawMask_l100_s1.fa > /home/ayasha/projects/def-mshafer/genome/intervals/NMBU_${INTERVAL}_mask_l100_s1_r0.9.fa

rm /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta.fna
rm /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1
rm /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1.sai
rm /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_split_l100_s1_aln-se.sam
rm /home/ayasha/projects/def-mshafer/genome/intervals/${INTERVAL}_fasta_rawMask_l100_s1.fa


## After running this, cat all output files together into one final mask file, then run the python script from msmc to generate bed files per chromosome, which can also be cat'd together into the final mappability mask



