module load SAMtools
module load Java

# Create Dictionary file for reference fasta
# ~/gatk-4.2.4.1/gatk CreateSequenceDictionary -R GCF_001858045.1_ASM185804v2_genomic_edit.fna
~/gatk-4.2.4.1/gatk CreateSequenceDictionary -R GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna

# Create index file for reference fasta
# samtools faidx GCF_001858045.1_ASM185804v2_genomic_edit.fna
samtools faidx GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna

# bwa index also
# bwa index GCF_001858045.1_ASM185804v2_genomic_edit.fna
bwa index GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna

# Create interval files for 40 intervals or 1000 intervals, or 80 (to fit in 6hour que)

# ~/gatk-4.2.4.1/gatk SplitIntervals -R GCF_001858045.1_ASM185804v2_genomic_edit.fna -scatter 80 -O intervals_ASM185804v2_80x

# ~/gatk-4.2.4.1/gatk SplitIntervals -R GCF_001858045.1_ASM185804v2_genomic_edit.fna -scatter 1000 -O intervals_ASM185804v2_1000x

~/gatk-4.2.4.1/gatk SplitIntervals -R GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna -scatter 80 -O intervals_UMD_NMBU_80x

~/gatk-4.2.4.1/gatk SplitIntervals -R GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna -scatter 1000 -O intervals_UMD_NMBU_1000x

