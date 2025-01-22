#!/bin/bash

##### SLURM #####
#SBATCH --job-name=sort_bam
#SBATCH --partition=express
#SBATCH --time=01:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-28
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load SAMtools/1.9-foss-2018b

##### VARIABLES #####
bam_dir="/data/scratch/flanary/atac-seq/untreated/GSE138293/bam"
replicate_list="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/GSE138293/06_format_bam/sample_replicates.txt"
replicate=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$replicate_list")

##### COMMANDS #####
# Sort by coordinates
samtools sort -o $bam_dir/"$replicate"_sort.bam $bam_dir/"$replicate"_mapped.bam

##### END #####
echo "Finished sorting BAM for: $replicate"