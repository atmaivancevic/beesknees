# Download bee reference genome(s)

# Amel_4.5

srun --pty bash

cd /Users/CL_Shared/db/genomes
mkdir Amel_4.5 
cd Amel_4.5

mkdir assembled_chromosomes
cd assembled_chromosomes
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/195/GCA_000002195.1_Amel_4.5/GCA_000002195.1_Amel_4.5_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/*.gz'

mkdir ../placed_scaffolds
cd ../placed_scaffolds
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/195/GCA_000002195.1_Amel_4.5/GCA_000002195.1_Amel_4.5_assembly_structure/Primary_Assembly/placed_scaffolds/FASTA/*.gz'

mkdir ../unplaced_scaffolds
cd ../unplaced_scaffolds
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/195/GCA_000002195.1_Amel_4.5/GCA_000002195.1_Amel_4.5_assembly_structure/Primary_Assembly/unplaced_scaffolds/FASTA/*.gz'

cd ../assembled_chromosomes
gunzip *

# rename headers to match UCSC/BeeBase (e.g. the header for chrLG1.fa should be >Group1, not >CM000054.5 Apis mellifera strain DH4 linkage group 1)
for i in {1..16}; do cat chrLG"$i".fna | grep -v ">" | sed 1i">Group$i" > chrLG"$i".renamed.fa; done

# concatenate renamed chr files into one main.fa file
chromList=$(ls *.renamed.fa | sort -V)
cat $chromList > Amel_4.5.main.fa

# move to new dir
mkdir ../main
mv Amel_4.5.main.fa  ../main

# make a chrom sizes file 
cd ../main
module load samtools
samtools faidx Amel_4.5.main.fa
cut -f1,2 Amel_4.5.main.fa.fai > Amel_4.5.chrom.sizes

# make hisat and bwa indexes
module load hisat2
hisat2-build Amel_4.5.main.fa Amel_4.5.main
# created idx files: 
# Amel_4.5.main.1.ht2 
# Amel_4.5.main.2.ht2 
# Amel_4.5.main.3.ht2 
# Amel_4.5.main.4.ht2 
# Amel_4.5.main.5.ht2 
# Amel_4.5.main.6.ht2 
# Amel_4.5.main.7.ht2
# Amel_4.5.main.8.ht2

module load bwa
bwa index Amel_4.5.main.fa Amel_4.5.main
# created idx files:
# Amel_4.5.main.fa.amb
# Amel_4.5.main.fa.ann
# Amel_4.5.main.fa.bwt
# Amel_4.5.main.fa.fai
# Amel_4.5.main.fa.pac
# Amel_4.5.main.fa.sa

# rearrange dirs and files to match hg38 layout
# tree:

# /Users/CL_Shared/db/genomes/Amel_4.5
# .
# ├── fa
# │   ├── Amel_4.5.chrom.sizes
# │   └── Amel_4.5.main.fa
# ├── index
# │   ├── bwa
# │   │   ├── Amel_4.5.main.fa.amb
# │   │   ├── Amel_4.5.main.fa.ann
# │   │   ├── Amel_4.5.main.fa.bwt
# │   │   ├── Amel_4.5.main.fa.fai
# │   │   ├── Amel_4.5.main.fa.pac
# │   │   └── Amel_4.5.main.fa.sa
# │   └── hisat2
# │       ├── Amel_4.5.main.1.ht2
# │       ├── Amel_4.5.main.2.ht2
# │       ├── Amel_4.5.main.3.ht2
# │       ├── Amel_4.5.main.4.ht2
# │       ├── Amel_4.5.main.5.ht2
# │       ├── Amel_4.5.main.6.ht2
# │       ├── Amel_4.5.main.7.ht2
# │       └── Amel_4.5.main.8.ht2
# └── raw
#     └── full_assembly
#         ├── assembled_chromosomes
#         │   ├── chrLG10.fna.gz
#         │   ├── chrLG11.fna.gz
#         │   ├── chrLG12.fna.gz
#         │   ├── chrLG13.fna.gz
#         │   ├── chrLG14.fna.gz
#         │   ├── chrLG15.fna.gz
#         │   ├── chrLG16.fna.gz
#         │   ├── chrLG1.fna.gz
#         │   ├── chrLG2.fna.gz
#         │   ├── chrLG3.fna.gz
#         │   ├── chrLG4.fna.gz
#         │   ├── chrLG5.fna.gz
#         │   ├── chrLG6.fna.gz
#         │   ├── chrLG7.fna.gz
#         │   ├── chrLG8.fna.gz
#         │   └── chrLG9.fna.gz
#         ├── placed_scaffolds
#         │   ├── chrLG10.placed.scaf.fna.gz
#         │   ├── chrLG11.placed.scaf.fna.gz
#         │   ├── chrLG12.placed.scaf.fna.gz
#         │   ├── chrLG13.placed.scaf.fna.gz
#         │   ├── chrLG14.placed.scaf.fna.gz
#         │   ├── chrLG15.placed.scaf.fna.gz
#         │   ├── chrLG16.placed.scaf.fna.gz
#         │   ├── chrLG1.placed.scaf.fna.gz
#         │   ├── chrLG2.placed.scaf.fna.gz
#         │   ├── chrLG3.placed.scaf.fna.gz
#         │   ├── chrLG4.placed.scaf.fna.gz
#         │   ├── chrLG5.placed.scaf.fna.gz
#         │   ├── chrLG6.placed.scaf.fna.gz
#         │   ├── chrLG7.placed.scaf.fna.gz
#         │   ├── chrLG8.placed.scaf.fna.gz
#         │   └── chrLG9.placed.scaf.fna.gz
#         └── unplaced_scaffolds
#             └── unplaced.scaf.fna.gz
# 9 directories, 49 files

###############################################################################

# Amel_HAv3.1

# Do the same thing for this newer bee assembly (ftp link on NCBI: ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_assembly_structure)

# final tree:

# /Users/CL_Shared/db/genomes/Amel_HAv3.1
# .
# ├── fa
# │   ├── Amel_HAv3.1.chrom.sizes
# │   └── Amel_HAv3.1.main.fa
# ├── index
# │   ├── bwa
# │   │   ├── Amel_HAv3.1.main.fa.amb
# │   │   ├── Amel_HAv3.1.main.fa.ann
# │   │   ├── Amel_HAv3.1.main.fa.bwt
# │   │   ├── Amel_HAv3.1.main.fa.fai
# │   │   ├── Amel_HAv3.1.main.fa.pac
# │   │   └── Amel_HAv3.1.main.fa.sa
# │   └── hisat2
# │       ├── Amel_HAv3.1.main.1.ht2
# │       ├── Amel_HAv3.1.main.2.ht2
# │       ├── Amel_HAv3.1.main.3.ht2
# │       ├── Amel_HAv3.1.main.4.ht2
# │       ├── Amel_HAv3.1.main.5.ht2
# │       ├── Amel_HAv3.1.main.6.ht2
# │       ├── Amel_HAv3.1.main.7.ht2
# │       └── Amel_HAv3.1.main.8.ht2
# └── raw
#     └── full_assembly
#         ├── assembled_chromosomes
#         │   ├── chrLG10.fna.gz
#         │   ├── chrLG11.fna.gz
#         │   ├── chrLG12.fna.gz
#         │   ├── chrLG13.fna.gz
#         │   ├── chrLG14.fna.gz
#         │   ├── chrLG15.fna.gz
#         │   ├── chrLG16.fna.gz
#         │   ├── chrLG1.fna.gz
#         │   ├── chrLG2.fna.gz
#         │   ├── chrLG3.fna.gz
#         │   ├── chrLG4.fna.gz
#         │   ├── chrLG5.fna.gz
#         │   ├── chrLG6.fna.gz
#         │   ├── chrLG7.fna.gz
#         │   ├── chrLG8.fna.gz
#         │   └── chrLG9.fna.gz
#         ├── placed_scaffolds
#         │   ├── chrLG10.placed.scaf.fna.gz
#         │   ├── chrLG11.placed.scaf.fna.gz
#         │   ├── chrLG12.placed.scaf.fna.gz
#         │   ├── chrLG13.placed.scaf.fna.gz
#         │   ├── chrLG14.placed.scaf.fna.gz
#         │   ├── chrLG15.placed.scaf.fna.gz
#         │   ├── chrLG16.placed.scaf.fna.gz
#         │   ├── chrLG1.placed.scaf.fna.gz
#         │   ├── chrLG2.placed.scaf.fna.gz
#         │   ├── chrLG3.placed.scaf.fna.gz
#         │   ├── chrLG4.placed.scaf.fna.gz
#         │   ├── chrLG5.placed.scaf.fna.gz
#         │   ├── chrLG6.placed.scaf.fna.gz
#         │   ├── chrLG7.placed.scaf.fna.gz
#         │   ├── chrLG8.placed.scaf.fna.gz
#         │   └── chrLG9.placed.scaf.fna.gz
#         ├── unlocalized_scaffolds
#         │   ├── chrLG10.unlocalized.scaf.fna.gz
#         │   ├── chrLG11.unlocalized.scaf.fna.gz
#         │   ├── chrLG12.unlocalized.scaf.fna.gz
#         │   ├── chrLG14.unlocalized.scaf.fna.gz
#         │   ├── chrLG15.unlocalized.scaf.fna.gz
#         │   ├── chrLG2.unlocalized.scaf.fna.gz
#         │   ├── chrLG5.unlocalized.scaf.fna.gz
#         │   ├── chrLG6.unlocalized.scaf.fna.gz
#         │   ├── chrLG7.unlocalized.scaf.fna.gz
#         │   ├── chrLG8.unlocalized.scaf.fna.gz
#         │   └── chrLG9.unlocalized.scaf.fna.gz
#         └── unplaced_scaffolds
#             └── unplaced.scaf.fna.gz

# 10 directories, 60 files

