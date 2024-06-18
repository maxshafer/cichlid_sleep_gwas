#!/bin/bash

#SBATCH --job-name=Mapmask                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=64G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/Mappabilitymask_new_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/cichlids_2/scripts/logs/Mappabilitymask_new_stderr.txt

#SBATCH --array=1-23%23
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

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
file_list="/scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.chromosomes"

# this is the second column of index_array
INTERVAL=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list}" | cut -f 1 -d ','`

# Subset the whole genome for only the interval (e.g. chromosome) and place it in the temp drive for manipulation
# Will need to load in samtools, subset, then purge modules before loading bwa
module load SAMtools/1.16.1-GCC-10.3.0

samtools faidx /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna $INTERVAL > ./${INTERVAL}_fasta.fna

module purge

module load BWA/0.7.17-goolf-1.7.20

bwa index ${INTERVAL}_fasta.fna


# Run
# Use the github compiled program (or the old binary from blogpost '~/seqbility-20091110/splitfa')
/scicore/home/schiera/gizevo30/seqbility-20091110/splitfa ${INTERVAL}_fasta.fna 100 > ${INTERVAL}_fasta_split_l100_s1

bwa aln -R 1000000 -O 3 -E 3 /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna ${INTERVAL}_fasta_split_l100_s1 > ${INTERVAL}_fasta_split_l100_s1.sai

# samse needs also the index and the input reads, as well as the output of aln
bwa samse /scicore/home/schiera/gizevo30/projects/cichlids_2/genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna ${INTERVAL}_fasta_split_l100_s1.sai ${INTERVAL}_fasta_split_l100_s1 > ${INTERVAL}_fasta_split_l100_s1_aln-se.sam

# generate the raw and non raw mask
/scicore/home/schiera/gizevo30/seqbility-20091110/gen_raw_mask.pl ${INTERVAL}_fasta_split_l100_s1_aln-se.sam > ${INTERVAL}_fasta_rawMask_l100_s1.fa

/scicore/home/schiera/gizevo30/seqbility-20091110/gen_mask -l 100 -r 0.9 ${INTERVAL}_fasta_rawMask_l100_s1.fa > /scicore/home/schiera/gizevo30/projects/cichlids_2/masks/NMBU_${INTERVAL}_mask_l100_s1_r0.9.fa

rm ${INTERVAL}_fasta.fna
rm ${INTERVAL}_fasta_split_l100_s1
rm ${INTERVAL}_fasta_split_l100_s1.sai
rm ${INTERVAL}_fasta_split_l100_s1_aln-se.sam
rm ${INTERVAL}_fasta_rawMask_l100_s1.fa


## After running this, cat all output files together into one final mask file, then run the python script from msmc to generate bed files per chromosome, which can also be cat'd together into the final mappability mask



