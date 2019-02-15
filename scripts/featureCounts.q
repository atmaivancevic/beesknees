#!/bin/bash

## Script for running feature counts
## Date: 15 Feb 2019 

## Example usage:
## inDir=/Users/CL_Shared/data/atma/beenomics_rna/bam outDir=/Users/CL_Shared/data/atma/beenomics_rna/featureCounts gtfFile=/Users/CL_Shared/db/genomes/Amel_4.5/annotations/Amel_4.5.annotation.gtf feature=exon strandOption=2 sbatch featureCounts.q

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2:00:00
#SBATCH --mem=16GB

# Job name and output
#SBATCH -J featureCounts
#SBATCH -o /Users/%u/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/slurmErr/slurm-%A_%a.err

# Email notifications 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user="%u"@colorado.edu

# load modules
module load subread

# Run featureCounts, outputting a single table containing count information for all samples
pwd; hostname; date

echo "Starting featureCounts..."
echo "featureCounts version: "$(featureCounts -v)

inputBams=${inDir}/*.uniq.bam

featureCounts -T 1 -O -p -s $strandOption -t $feature -g gene_id -a ${gtfFile} -o ${outDir}/featureCounts.txt $inputBams

echo $(date +"[%b %d %H:%M:%S] Done!")

## Explanation of arguments:
# '-O' - assign reads to all their overlapping meta-features
# '-T <int>' - number of threads; 1 by default
# '-s <int>' - perform strand-specific read counting; options: 0 (unstranded), 1 (stranded), and 2 (reversely stranded); 0 by default
# '-g <option>' - specify attribute type in GTF annotation; 'gene_id' by default
# '-a <file.gtf' - name of annotation file; gtf format by default
