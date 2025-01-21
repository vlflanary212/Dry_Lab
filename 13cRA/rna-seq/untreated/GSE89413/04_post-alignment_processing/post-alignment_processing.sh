#!/bin/bash

##### SLURM #####
#SBATCH --job-name=samtools
#SBATCH --partition=medium
#SBATCH --time=23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-40
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES #####
module load SAMtools/1.9-foss-2018b

###### ARRAY #######
samples="/home/flanary/Dry_Lab/13cRA/rna-seq/untreated/GSE89413/SRR_Acc_List.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$samples")

##### VARIABLES #####
bam_dir="/data/scratch/flanary/rna-seq/untreated/GSE89413/bam"
output=$bam_dir/$sample

##### COMMANDS #####
echo "Processing sample: $sample"

# index
samtools index $output/"$sample".Aligned.sortedByCoord.out.bam

# filter bam files
samtools view -b -F 4 -q 20 $output/"$sample".Aligned.sortedByCoord.out.bam > $output/$sample.bam

# remove duplicates
samtools rmdup -S $output/$sample.bam $output/$sample.rmdup.bam

# index
samtools index $output/$sample.rmdup.bam

##### END #####
echo "Finish post-alignment processing with SAMtools for $sample"
