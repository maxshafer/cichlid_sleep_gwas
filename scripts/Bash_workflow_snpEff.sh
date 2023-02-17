## Feb 14th, 2023
## MERS

### OK I want to annotate all of the SNPs in the vcf files
### My script for annotating is quite slow, even with some improvements
### If I want to do GSEA, I need the full list of genes, so I would take the lowest p-value per gene from the annotation

## I want to use snpeff, which is now on the cluster (I installed it, so to keep the same version etc)
## Problem is, the version of the genome I used is old, and doesn't have a snpEff database
## I can create my own database, or alternatively, see if I can just convert the vcf files to be compatible (change chr identifiers etc)

## Start an srun

srun --nodes=1 --cpus-per-task=8 --mem=64G --pty bash

## snpeff is installed in my home drive, in the /bin directory
module load Java

## Check which databases exist with snpeff

java -jar ~/bin/snpEff/snpEff.jar databases > snp_databases.txt

## Search for specific ones (O. niloticus)

grep Oreochromis snp_databases.txt

## downloading the databases for Oreochromis niloticus from snpeff

java -jar snpEff.jar download O_niloticus_UMD_NMBU.99
java -jar ~/bin/snpEff/snpEff.jar download -v Orenil1.0.75
java -jar ~/bin/snpEff/snpEff.jar download -v O_niloticus_UMD_NMBU.105

## snpeff runs on gz files as well
## This doesn't match anything though, because the chr identifiers are wrong

java -Xmx8g -jar ~/bin/snpEff/snpEff.jar -v O_niloticus_UMD_NMBU.105 cohort_0060-scattered.interval_list.g.vcf.gz > test.output.ann.105.vcf

## This runs sed to replace a chr identifier in a gunzipped vcf file

gunzip cohort_0060-scattered.interval_list.g.vcf.gz

# Replaces 'NC_031985' with 'LG22', which might not be the correct match, I should check
cat cohort_0060-scattered.interval_list.g.vcf | sed "s/^NC_031985/LG22/" > 0063_test.vcf

# This runs the annotation on the gunzipped test vcf (its actually 0060, but named wrongly)
# Doesn't take so long, and does produce results, but the mappings are wrong for both the 99 or 105 genome (from Ensembl)
# For example, it will annotate a SNP, but that position is near a gene in the ensembl database, but not in my genome

java -Xmx8g -jar ~/bin/snpEff/snpEff.jar -v O_niloticus_UMD_NMBU.99 0063_test.vcf > test.output.ann.99.vcf
java -Xmx8g -jar ~/bin/snpEff/snpEff.jar -v O_niloticus_UMD_NMBU.105 0063_test.vcf > test.output.ann.105.vcf


## I really think matching my chrs to some database isn't going to work, so I will make my own
## Here are the steps from https://pcingola.github.io/SnpEff/se_buildingdb/#building-a-database

# Step 1: Configure a new genome in SnpEff's config file snpEff.config.
# 
# Add genome entry to snpEff's configuration
# If the genome uses a non-standard codon table: Add codon table parameter
# Step 2: Build using gene annotations and reference sequences
# 
# Option 1: Building a database from GTF files (recommended for large genomes)
# Option 2: Building a database from GenBank files (recommended for small genomes)
# Option 3: Building a database from GFF files
# Option 4: Building a database from RefSeq table from UCSC
# Step 3: Checking the database: SnpEff will check the database by comparing predicted protein sequences and CDS sequences with ones provided by the user.
# 
# Checking CDS sequences
# Checking Protein sequences

## I will have to use option 3, because that's the only annotation file available for my genome from Refseq (newer versions have gtf)

# Edit the config file to create the new genome:
# use grep to find where the other oreochromis genomes are put

grep Oreochromis ~/bin/snpEff/snpEff.config -n

vi +1185 ~/bin/snpEff/snpEff.config 

## I'll call it ASM185804v2, like my gff files etc
# Add: 'ASM185804v2.genome : Oreochromis_niloticus'

# For now, I will not configure codon tables, the ones for vertebrates should suffice?

## Step 2, build database from GFF file
## First to get the gff file, it's already on cluster

mkdir ~/bin/snpEff/data/ASM185804v2
cd ~/bin/snpEff/data/ASM185804v2

cp ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic.gff genes.gff

## Note: GFF3 files can include the reference sequence in the same file. This is done by dumping the fasta file after a '##FASTA' line. You can also add the sequence fasta file to the 'data/genomes/' directory, like it is done in when using GTF format.

cp ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic_edit.fna sequences.fa

cd ~/bin/snpEff
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -v ASM185804v2

## The above runs but

# There might be differences in the chromosome names used in the genes file ('/scicore/home/schiera/gizevo30/bin/snpEff/./data/ASM185804v2/genes.gff')
# and the chromosme names used in the 'reference sequence' file.
# Please check that chromosome names in both files match.
# 	Chromosome names missing in 'reference sequence' file:	'NC_013663.1', 'NC_031965.1',

# I hate the '.1's, ugh, ok maybe I can just use a sed command, I know in my vcf they don't have .1s, so let's make this consistant with that
# Seems like the fasta matches the vcfs, but the gff does not, so lets sed the gff
# Might be safe to just sed out the '.1', or can do a for loop - maybe check how many lines match teh .1?

# (base) [gizevo30@login01 ASM185804v2]$ head genes.gff 
# ##gff-version 3
# #!gff-spec-version 1.21
# #!processor NCBI annotwriter
# #!genome-build ASM185804v2
# #!genome-build-accession NCBI_Assembly:GCF_001858045.1
# #!annotation-source NCBI Oreochromis niloticus Annotation Release 103
# ##sequence-region NC_031965.1 1 38372991
# ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=8128
# NC_031965.1	Gnomon	gene	11049	13118	.	-	.	ID=gene0;Dbxref=GeneID:102075536;Name=LOC102075536;gbkey=Gene;gene=LOC102075536;gene_biotype=lncRNA
# NC_031965.1	Gnomon	lnc_RNA	11049	13118	.	-	.	ID=rna0;Parent=gene0;Dbxref=GeneID:102075536,Genbank:XR_266817.3;Name=XR_266817.3;gbkey=ncRNA;gene=LOC102075536;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 13 samples with support for all annotated introns;product=uncharacterized LOC102075536;transcript_id=XR_266817.3
# 
# (base) [gizevo30@login01 ASM185804v2]$ head sequences.fa 
# >NC_031965
# AAAATAgctgaacatgtttaaaatataaattagaaaaagatgagtaatgacaaaagaaaatttagtCAGAAAACCTCTGA
# ATGAATAACAAGTTCATATTCTCTAATCACTTAAAGAGTGGAATTATGTATTATAAATCtctaattcataaaaaaataat
# aaaaataaatgtaaaacaaatttgttgaattgaactgaaaagatGAACAAACATTCTGGAGAAAACACAAGCTTTAATAT
# ACAGCATAATGTTTTCAGACTATACATGTGTAATCATGTTAATGGTATTACATACATCAATTTGTTTTGTATGTCAtaga
# aataacataaaaatcttAAGTCATATTGTACAGCTTCTTTCTGAAAAGGACTTAAATTAAttcaacaaatgttttttttt
# ttgaaatagaataaaaatttttaaaaagagacgtttgcAATGTTTAAGGTGTTAATAATCTGATCCTGTCATGTGTCTGT
# TCAACTTTAGTTTGGCACAGACTCATCTTAAAGAACAAATTaccaaataataaacaaataaattgaaattaaagtagctt
# ttttgttaaacttcacaggagaataaataaaaataactaaattaaataaccaaataaaataataagcaacataatacata
# aaattcAACTTTTAAAACCACAATCCTGAACCTTTAAATTGTGGGTTGGTTTCTTAGAACATGGCTGAACAGCTGTTTCT
# (base) [gizevo30@login01 ASM185804v2]$ 

## This matches the chr identifiers, but also many gene ids, there is a lot of stuff in the gff file (less would be in the fasta)

grep '\.1' genes.gff

grep 'NC_031967\.1' genes.gff

## However, they are only ever in the first column! I can use awk and sed together!
## This doesn't work correctly

awk '{gsub("\.1","",$0)}1' genes.gff > test.gff

## OK this is super close, but removes '.1' from gene ids, so the '$0' doesn't seem to work with this file
## I don't know why the above also removes some of the numbers ('03') from the chr ID...

awk 'BEGIN{FS=OFS="NC_......"} {gsub(/\.1/, "", $0)} 1' genes.gff > test.gff

## Aha! this seems to work! replaces the first occurance of .1 in each line!
sed 's/\.1//1' genes.gff > test.gff

# it does change this line '#!genome-build-accession NCBI_Assembly:GCF_001858045.1' and removes the '.1', but lets see if it builds now


cd ~/bin/snpEff
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -v ASM185804v2


## OK this runs, but fails the check for CDS and protein sequences - can I get these from Refseq?
## Seems so https://www.ncbi.nlm.nih.gov/assembly/GCF_001858045.1/
## These have been downloaded, now rsync into the cluster

rsync ~/Downloads/genome_assemblies_cds_fasta/ncbi-genomes-2023-02-14/GCF_001858045.1_ASM185804v2_cds_from_genomic.fna.gz  gizevo30@login.scicore.unibas.ch:/scicore/home/schiera/gizevo30/projects/cichlids_2/genome/
rsync ~/Downloads/genome_assemblies_prot_fasta/ncbi-genomes-2023-02-14/GCF_001858045.1_ASM185804v2_protein.faa.gz gizevo30@login.scicore.unibas.ch:/scicore/home/schiera/gizevo30/projects/cichlids_2/genome/


## Copy to the folder, then try the build again

cp ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_cds_from_genomic.fna ~/bin/snpEff/data/ASM185804v2/cds.fa
cp ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_protein.faa ~/bin/snpEff/data/ASM185804v2/protein.fa

# Try build again

cd ~/bin/snpEff
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -v ASM185804v2

## Ok different error this time, but I think it's because the new fasta files have the '.1'
## This makes me think it might be easier to leave the '.1', and change my vcf files?
## The matching will work better, because I will write a script that uses the whole chr names, and I will just ignore anything on the contigs (like Milan does)
## Yeah, that seems right

# Re copy the gff file (with '.1')

cp ~/projects/cichlids_2/genome/GCF_001858045.1_ASM185804v2_genomic.gff genes.gff

## But now, need to change the sequence fasta, or maybe re-download it? Maybe it will have the '.1' AND/OR try the build again (doesn't work, non-matching chr IDs)


rsync ~/Downloads/genome_assemblies_genome_fasta/ncbi-genomes-2023-02-14/GCF_001858045.1_ASM185804v2_genomic.fna.gz  gizevo30@login.scicore.unibas.ch:/scicore/home/schiera/gizevo30/bin/snpEff/data/ASM185804v2/

## This now matches, and has ".1" attached to chr IDs, now re try the build

cd ~/bin/snpEff
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -v ASM185804v2


### Hmm I still get this:

FATAL ERROR: No CDS checked. This is might be caused by differences in FASTA file transcript IDs respect to database's transcript's IDs.
Transcript IDs from database (sample):
	'rna25981'
	'rna70978'
	'rna6538'
	'rna6537'
	'rna41728'
	'rna63592'
	'rna6536'
	'rna24040'
	'rna41729'
	'rna67'
	'rna24041'
	'rna67042'
	'rna66'
	'rna3987'
	'rna67041'
	'rna43768'
	'rna3985'
	'rna3986'
	'rna6534'
	'rna6535'
	'rna24042'
	'rna29190'
Transcript IDs from database (fasta file):
	'1_cds_XP_005471655'
	'lcl|NC_031976.1_cds_XP_013131958.1_26186'
	'lcl|NC_031966.1_cds_XP_003443474.2_2345'
	'lcl|NC_031969.1_cds_XP_019213678.1_8851'
	'lcl|NC_031981.1_cds_XP_019201999.1_41098'
	'lcl|NC_031983.1_cds_XP_003442956.1_45682'
	'lcl|NC_031972.1_cds_XP_019216573.1_16306'
	'1_cds_XP_005471650'
	'1_cds_XP_005471651'
	'lcl|NC_031972.1_cds_XP_005470592.1_17290'
	'1_cds_XP_005471643'
	'1_cds_XP_005471644'
	'1_cds_XP_005471647'
	'1_cds_XP_005471648'
	'1_cds_XP_005471649'
	'lcl|NW_017613886.1_cds_XP_019209092.1_54439'
	'1_cds_XP_005471640'
	'1_cds_XP_005471642'
	'1_cds_XP_005471676'
	'1_cds_XP_005471677'
	'lcl|NC_031985.1_cds_XP_013123567.2_51323'
	'1_cds_XP_005471679'



## I don't know why this is happening, everything comes from the same place? I've download the cds, protein, and genome sequences today - maybe I should also get the gff freshly?


rsync ~/Downloads/genome_assemblies_genome_gff/ncbi-genomes-2023-02-14/GCF_001858045.1_ASM185804v2_genomic.gff  gizevo30@login.scicore.unibas.ch:/scicore/home/schiera/gizevo30/bin/snpEff/data/ASM185804v2/


cd ~/bin/snpEff
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -v ASM185804v2

## Still doesn't take, don't know what is going on, but the IDs don't match between the gff and the cds
## I could try option #2, and build from genbank? I think the instructions are a bit defunk

# Try build again with option -gff2

java -jar ~/bin/snpEff/snpEff.jar build -gff2 -v ASM185804v2

# nope :(

## Lets try the same but with -gtf22

mv genes.gff genes.gtf
cd ~/bin/snpEff
java -jar ~/bin/snpEff/snpEff.jar build -gtf22 -v ASM185804v2

# This may have borken the cluster

## Alternatively, I can run it without checking, then try the annotation to see if it works
## Maybe the database will match my vcf files, but not the cds? I will have to add '.1' to chr IDs in my vcf though


cd ~/bin/snpEff
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -noCheckCds -v ASM185804v2

# and/or no check protein
java -jar ~/bin/snpEff/snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v ASM185804v2


## OK this works, and finishes without errors, lets try to run it on a vcf

# Replaces 'NC_031985' with 'NC_031985.1', which should match my new database
cat cohort_0060-scattered.interval_list.g.vcf | sed "s/^NC_031985/NC_031985.1/" > 0060_test.vcf

java -Xmx8g -jar ~/bin/snpEff/snpEff.jar -v ASM185804v2 0060_test.vcf > test.output.ann.ASM.vcf


## Hey, this seems to actually work, and produces something that is correct
## This is a line fromm the vcf file 'test.out.ann.ASM.vcf', and it matches the ncbi genome browser
## But, the warning "WARNING_TRANSCRIPT_NO_START_CODON", might actually be informative, maybe I haven't assigned a codon table to the new database!


ANN=T|intron_variant|MODIFIER|id695969|GENE_id695969|transcript|rna58924|protein_coding|8/8|c.1140+603C>T||||||,T|non_coding_transcript_variant|MODIFIER|LOC100694011|gene30699|transcript|XM_019350307.1|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON

NC_031985.1     13912761        .       C       T       10980.57        .       AC=152;AF=0.639;AN=238;DP=305;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=0.4154;MLEAC=230;MLEAF=0.966;MQ=37.86;QD=30.28;SOR=1.580;ANN=T|intron_variant|MODIFIER|id695969|GENE_id695969|transcript|rna58924|protein_coding|8/8|c.1140+603C>T||||||,T|non_coding_transcript_variant|MODIFIER|LOC100694011|gene30699|transcript|XM_019350307.1|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON      GT:AD:DP:GQ:PGT:PID:PL:PS       1/1:0,2:2:6:.:.:51,6,0  1/1:0,3:3:9:.:.:121,9,0 0/0:1,0:1:0:.:.:0,0,0   1/1:0,2:2:6:.:.:57,6,0  1|1:0,5:5:15:1|1:13912759_T_C:225,15,0:13912759 0/0:1,0:1:0:.:.:0,0,0   1|1:0,1:1:3:1|1:13912741_GC_G:26,3,0:13912741   1/1:0,3:3:9:.:.:93,9,0  1/1:0,2:2:6:.:.:63,6,0  0/0:0,0:0:0:.:.:0,0,0   0/0:0,0:0:0:.:.:0,0,0   1/1:0,2:2:6:.:.:41,6,0  0/0:0,0:0:0:.:.:0,0,0   1/1:0,2:2:6:.:.:74,6,0  0/0:0,0:0:0:.:.:0,0,0   1/1:0,3:3:9:.:.:125,9,0 1/1:0,2:2:6:.:.:80,6,0  1/1:0,8:8:24:.:.:314,24,0       1/1:0,2:2:6:.:.:57,6,0  0/0:1,0:1:0:.:.:0,0,0   1/1:0,4:4:12:.:.:137,12,0       0/0:2,0:2:0:.:.:0,0,0   0/0:4,0:4:0:.:.:0,0,0   0/0:4,0:4:0:.:.:0,0,0   1/1:0,4:4:12:.:.:157,12,0       0/0:1,0:1:0:.:.:0,0,0   0/0:1,0:1:0:.:.:0,0,0   1/1:0,3:3:9:.:.:115,9,0 1/1:0,5:5:15:.:.:175,15,0       0/0:1,0:1:0:.:.:0,0,0   1/1:0,2:2:6:.:.:74,6,0  1|1:0,3:3:9:1|1:13912759_T_C:131,9,0:13912759   1/1:0,4:4:12:.:.:121,12,0       1/1:0,2:2:6:.:.:87,6,0  1/1:0,5:5:15:.:.:167,15,0       0/0:1,0:1:0:.:.:0,0,0   1/1:0,3:3:9:.:.:95,9,0  1/1:0,4:4:12:.:.:143,12,0       1/1:0,4:4:12:.:.:135,12,0       1/1:0,3:3:9:.:.:99,9,0  0/0:2,0:2:0:.:.:0,0,0   1/1:0,2:2:6:.:.:74,6,0  1/1:0,5:5:15:.:.:173,15,0       0/0:0,0:0:0:.:.:0,0,0   1/1:0,3:3:9:.:.:105,9,0 1/1:0,2:2:6:.:.:57,6,0  1|1:0,2:2:6:1|1:13912741_GC_G:77,6,0:13912741   0/0:1,0:1:0:.:.:0,0,0   1/1:0,4:4:12:.:.:143,12,0       1/1:0,5:5:14:.:.:114,14,0       1/1:0,3:3:9:.:.:88,9,0  0/0:1,0:1:0:.:.:0,0,0   1/1:0,2:2:6:.:.:57,6,0  1/1:0,3:3:9:.:.:67,9,0  0/0:0,0:0:0:.:.:0,0,0   1/1:0,3:3:9:.:.:95,9,0  1/1:0,4:4:12:.:.:109,12,0       0/0:2,0:2:0:.:.:0,0,0   1/1:0,3:3:9:.:.:100,9,0 0/0:0,0:0:0:.:.:0,0,0   0/0:1,0:1:0:.:.:0,0,0   0/0:0,0:0:0:.:.:0,0,0   0/0:0,0:0:0:.:.:0,0,0   0/0:1,0:1:0:.:.:0,0,0   0/0:1,0:1:0:.:.:0,0,0   1/1:0,3:3:9:.:.:93,9,0  1/1:0,2:2:6:.:.:78,6,0  0/0:0,0:0:0:.:.:0,0,0   1/1:0,5:5:15:.:.:205,15,0       0/0:1,0:1:0:.:.:0,0,0   0/0:1,0:1:0:.:.:0,0,0   0/0:1,0:1:0:.:.:0,0,0   1/1:0,4:4:12:.:.:137,12,0       1/1:0,2:2:6:.:.:77,6,0  1/1:0,2:2:6:.:.:64,6,0  0/0:0,0:0:0:.:.:0,0,0   1/1:0,5:5:15:.:.:136,15,0       1/1:0,2:2:6:.:.:64,6,0  1/1:0,2:2:6:.:.:68,6,0  0/0:2,0:2:0:.:.:0,0,0   1/1:0,2:2:6:.:.:61,6,0  1/1:0,2:2:6:.:.:68,6,0  0/0:0,0:0:0:.:.:0,0,0   0/0:0,0:0:0:.:.:0,0,0   0/0:0,0:0:0:.:.:0,0,0   0/0:2,0:2:0:.:.:0,0,0   1/1:0,2:2:6:.:.:68,6,0  1/1:0,4:4:12:.:.:131,12,0       0/0:1,0:1:0:.:.:0,0,0   1|1:0,1:1:3:1|1:13912741_GC_G:32,3,0:13912741   1/1:0,2:2:6:.:.:51,6,0  1/1:0,2:2:6:.:.:61,6,0  0/0:1,0:1:0:.:.:0,0,0   1/1:0,8:8:24:.:.:284,24,0       1/1:0,5:5:15:.:.:175,15,0       1/1:0,4:4:12:.:.:141,12,0       1/1:0,4:4:12:.:.:157,12,0       1/1:0,2:2:6:.:.:57,6,0  1/1:0,3:3:9:.:.:83,9,0  1/1:0,7:7:21:.:.:252,21,0       1/1:0,2:2:6:.:.:68,6,0  1/1:0,1:1:3:.:.:42,3,0  1/1:0,4:4:12:.:.:120,12,0       1/1:0,4:4:12:.:.:121,12,0       1/1:0,3:3:9:.:.:121,9,0 0/0:1,0:1:0:.:.:0,0,0   1/1:0,3:3:9:.:.:95,9,0  0/0:2,0:2:0:.:.:0,0,0   1/1:0,4:4:12:.:.:125,12,0       1/1:0,3:3:9:.:.:122,9,0 1/1:0,3:3:9:.:.:106,9,0 0/0:0,0:0:0:.:.:0,0,0   1/1:0,4:4:12:.:.:121,12,0       1/1:0,3:3:9:.:.:109,9,0 1/1:0,2:2:6:.:.:74,6,0  1/1:0,3:3:9:.:.:119,9,0 0/0:1,0:1:0:.:.:0,0,0   1/1:0,6:6:18:.:.:208,18,0       0/0:1,0:1:0:.:.:0,0,0

### OK, I either can't assign a codonTable, or what I did doesn't work, or it doesn't need to be added
### I think I might go ahead with the unchecked version, it at least has the gene names


### Ahh ok, I think this is working, but there are issues with how it calls variant effects
### I can still retrieve genes and relative positions to SNPs









