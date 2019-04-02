#!/bin/bash

## Script for running hisat2
## Date: 13 Feb 2019 

## Example usage:
## inDir=/Users/CL_Shared/data/atma/beenomics_rna/trimmed_fastq outDir=/Users/CL_Shared/data/atma/beenomics_rna/bam hisatIdxDir=/Users/CL_Shared/db/genomes/Amel_4.5/index/hisat2 hisatIdx=Amel_4.5.main sbatch --array 0-7 hisat2_PE.q

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=4:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J hisat2
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@colorado.edu

# Set constant variables
rnaStrandness=RF # specify whether the library is forward (FR) or reverse stranded (RF)
numThreads=8

# Set query files
# note: remove R1/R2 to generate a unique identifier for each pair of files
queries=($(ls ${inDir}/*_paired.fastq.gz | xargs -n 1 basename | sed 's/_R1_paired.fastq.gz//g' | sed 's/_R2_paired.fastq.gz//g' | uniq))

# Load modules
module load hisat2 samtools

# Run hisat2
pwd; hostname; date

echo "Starting HISAT2 alignment..."
echo "HISAT2 version: "$(hisat2 --version)
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

hisat2 --no-softclip \
--rna-strandness ${rnaStrandness} \
-p ${numThreads} \
-x ${hisatIdxDir}/${hisatIdx} \
-1 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_paired.fastq.gz \
-2 ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_paired.fastq.gz \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_paired.fastq.gz}.uniq.bam

# Create a bai index file for upload to UCSC browser
samtools index ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%_paired.fastq.gz}.uniq.bam

echo $(date +"[%b %d %H:%M:%S] Done!")

## Explanation of arguments:
# '-p <int>' - number of threads
# '--rna-strandness <rf/fr>' - specify strand-specific information; library kit typically has this information
# '--no-softclip' - disallows soft clipping
# '-x <file.basename>' - hisat2 index basename; include everything up to .main.1.ht2
# '-1 <file_R1_paired.fastq.gz' - mate file 1
# '-2 <file_R2_paired.fastq.gz' - mate file 2
# '-q <int>' - filter out mapped reads by quality threshold
# '-Sb' - convert sam to bam
# '-@' - number of threads
