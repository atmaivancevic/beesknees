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

# download gene annotations
# gtf format
wget 'ftp://ftp.ensemblgenomes.org/pub/metazoa/current/gtf/apis_mellifera/Apis_mellifera.Amel_4.5.42.chr.gtf.gz'

# gff3 format
wget 'ftp://ftp.ensemblgenomes.org/pub/metazoa/release-42/gff3/apis_mellifera/Apis_mellifera.Amel_4.5.42.chr.gff3.gz'

# check the chromosome sizes (listed in the gff3) match what I have in the Amel4.5.chrom.sizes file

# fiji-1:/Users/CL_Shared/db/genomes/Amel_4.5/genes$ head Apis_mellifera.Amel_4.5.42.chr.gff3 -n 20
##gff-version 3
##sequence-region   1 1 29893408
##sequence-region   10 1 12965953
##sequence-region   11 1 14726556
##sequence-region   12 1 11902654
##sequence-region   13 1 10288499
##sequence-region   14 1 10253655
##sequence-region   15 1 10167229
##sequence-region   16 1 7207165
##sequence-region   2 1 15549267
##sequence-region   3 1 13234341
##sequence-region   4 1 12718334
##sequence-region   5 1 14363272
##sequence-region   6 1 18472937
##sequence-region   7 1 13219345
##sequence-region   8 1 13546544
##sequence-region   9 1 11120453

# fiji-1:/Users/CL_Shared/db/genomes/Amel_4.5/genes$ cat ../fa/Amel_4.5.chrom.sizes 
# Group1	29893408
# Group2	15549267
# Group3	13234341
# Group4	12718334
# Group5	14363272
# Group6	18472937
# Group7	13219345
# Group8	13546544
# Group9	11120453
# Group10	12965953
# Group11	14726556
# Group12	11902654
# Group13	10288499
# Group14	10253655
# Group15	10167229
# Group16	7207165

# Yay they're the same! 

# need to prefix the gene gtf with "Group" to make the names match the fasta chr files
cat Apis_mellifera.Amel_4.5.42.chr.gtf | awk '{print "Group" $0}' | sed 's/Group#/#/g' | head > Amel_4.5.annotation.gtf 
# use Amel_4.5.annotation.gtf  for featureCounts etc.

# download repeats

wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/195/GCA_000002195.1_Amel_4.5/GCA_000002195.1_Amel_4.5_rm.out.gz'

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

# repeats sourced from:
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_rm.out.gz'

# genes sourced from:
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Apis_mellifera/GFF/ref_Amel_HAv3.1_top_level.gff3.gz'

# as above, check that the chromosome sizes match in the gene files and full chr sizes file

# fiji-1:/Users/CL_Shared/db/genomes/Amel_HAv3.1/genes$ cat ref_Amel_HAv3.1_top_level.gff3 | grep Name=LG
# NC_037638.1	RefSeq	region	1	27754200	.	+	.	ID=id0;Dbxref=taxon:7460;Name=LG1;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG1;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037639.1	RefSeq	region	1	16089512	.	+	.	ID=id33275;Dbxref=taxon:7460;Name=LG2;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG2;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037640.1	RefSeq	region	1	13619445	.	+	.	ID=id57466;Dbxref=taxon:7460;Name=LG3;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG3;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037641.1	RefSeq	region	1	13404451	.	+	.	ID=id78108;Dbxref=taxon:7460;Name=LG4;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG4;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037642.1	RefSeq	region	1	13896941	.	+	.	ID=id97092;Dbxref=taxon:7460;Name=LG5;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG5;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037643.1	RefSeq	region	1	17789102	.	+	.	ID=id112983;Dbxref=taxon:7460;Name=LG6;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG6;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037644.1	RefSeq	region	1	14198698	.	+	.	ID=id129022;Dbxref=taxon:7460;Name=LG7;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG7;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037645.1	RefSeq	region	1	12717210	.	+	.	ID=id144059;Dbxref=taxon:7460;Name=LG8;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG8;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037646.1	RefSeq	region	1	12354651	.	+	.	ID=id164757;Dbxref=taxon:7460;Name=LG9;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG9;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037647.1	RefSeq	region	1	12360052	.	+	.	ID=id178189;Dbxref=taxon:7460;Name=LG10;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG10;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037648.1	RefSeq	region	1	16352600	.	+	.	ID=id191710;Dbxref=taxon:7460;Name=LG11;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG11;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037649.1	RefSeq	region	1	11514234	.	+	.	ID=id213107;Dbxref=taxon:7460;Name=LG12;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG12;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037650.1	RefSeq	region	1	11279722	.	+	.	ID=id224682;Dbxref=taxon:7460;Name=LG13;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG13;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037651.1	RefSeq	region	1	10670842	.	+	.	ID=id233955;Dbxref=taxon:7460;Name=LG14;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG14;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037652.1	RefSeq	region	1	9534514	.	+	.	ID=id246717;Dbxref=taxon:7460;Name=LG15;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG15;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body
# NC_037653.1	RefSeq	region	1	7238532	.	+	.	ID=id259304;Dbxref=taxon:7460;Name=LG16;collection-date=2003;country=USA;gbkey=Src;genome=chromosome;isolation-source=bee hive;linkage-group=LG16;mol_type=genomic DNA;sex=male;strain=DH4;tissue-type=whole body

# fiji-1:/Users/CL_Shared/db/genomes/Amel_HAv3.1/genes$ cat ../fa/Amel_HAv3.1.chrom.sizes 
# Group1	27754200
# Group2	16089512
# Group3	13619445
# Group4	13404451
# Group5	13896941
# Group6	17789102
# Group7	14198698
# Group8	12717210
# Group9	12354651
# Group10	12360052
# Group11	16352600
# Group12	11514234
# Group13	11279722
# Group14	10670842
# Group15	9534514
# Group16	7238532

###################

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

