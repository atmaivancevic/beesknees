#!/bin/bash

## Script to trim adapters using trimmomatic (paired end reads)
## Date: 13 Feb 2019
##
## Example usage:
## adapterFile=/opt/trimmomatic/0.36/adapters/TruSeq3-SE.fa inDir=/Users/CL_Shared/data/atma/beenomics_rna/fastq outDir=/Users/CL_Shared/data/atma/beenomics_rna/trimmed_fastq sbatch --array 0-7 trimAdapters_PE.q 

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=6:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J trimAdapters
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=atma.ivancevic@colorado.edu

# set constant variables
trimExe=/opt/trimmomatic/0.36/trimmomatic-0.36.jar
numCpus=8

# set query files
queries=($(ls $inDir/*.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//g' | sed 's/_R2_001.fastq.gz//g' | uniq))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo "Trimmomatic version: "$(java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar --version)
echo $(date +"[%b %d %H:%M:%S] Starting trimmomatic...")

java -jar $trimExe \
PE -threads $numCpus \
-trimlog ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_trimlog.txt \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1.fastq.gz ${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2.fastq.gz \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_paired.fastq.gz ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R1_unpaired.fastq.gz \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_paired.fastq.gz ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]}_R2_unpaired.fastq.gz \
ILLUMINACLIP:${adapterFile}:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo $(date +"[%b %d %H:%M:%S] Done!")

## Explanation of arguments:
# 'PE' - paired-end input
# '-threads <int>' - number of threads to use
# '-trimlog <file.txt>' - creates a log of all read trimmings
# 'ILLUMINACLIP <file.fa>' - cut adapter and other illumina-specific sequences from the read
# 'SLIDINGWINDOW <int:int>' - scans 5' end of reads and clips once average quality within the specified window falls below a threshold
# 'LEADING <int>' - cut bases off the start of a read, if below a threshold quality
# 'TRAILING <int>' - cut bases off the end of a read, if below a threshold quality
# 'MINLEN <int>' - retain only reads of some length (in bp)