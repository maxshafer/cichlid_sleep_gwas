module load SAMtools
module load Java

# Create Dictionary file for reference fasta
~/gatk-4.2.4.1/gatk CreateSequenceDictionary -R GCF_001858045.1_ASM185804v2_genomic_edit.fna

# Create index file for reference fasta
samtools faidx GCF_001858045.1_ASM185804v2_genomic_edit.fna

# bwa index also
bwa index GCF_001858045.1_ASM185804v2_genomic_edit.fna

# Create interval files for 40 intervals or 1000 intervals, or 80 (to fit in 6hour que)

~/gatk-4.2.4.1/gatk SplitIntervals -R GCF_001858045.1_ASM185804v2_genomic_edit.fna -scatter 80 -O intervals_ASM185804v2_80x

~/gatk-4.2.4.1/gatk SplitIntervals -R GCF_001858045.1_ASM185804v2_genomic_edit.fna -scatter 1000 -O intervals_ASM185804v2_1000x


