#!/bin/bash

##### SLURM #####
#SBATCH --job-name=merge_bam
#SBATCH --partition=short
#SBATCH --time=11:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load SAMtools/1.9-foss-2018b

##### VARIABLES #####
bam_dir="/data/scratch/flanary/atac-seq/untreated/GSE108517/bam"
sample_list="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/GSE108517/cell_lines.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$sample_list")

##### COMMANDS #####
# Merge sorted bam files for replicates
samtools merge $bam_dir/"$sample"_merge.bam $bam_dir/"$sample"_rep1_sort.bam $bam_dir/"$sample"_rep2_sort.bam

# remove duplicates
samtools rmdup -S $bam_dir/"$sample"_merge.bam $bam_dir/"$sample"_rmdup.bam

# index 
samtools index $bam_dir/"$sample"_rmdup.bam

# remove mitochondrial reads
samtools idxstats $bam_dir/"$sample"_rmdup.bam | cut -f 1 | grep -v chrM | xargs samtools view -b $bam_dir/"$sample"_rmdup.bam > $bam_dir/"$sample"_final.bam

# index again
samtools index $bam_dir/"$sample"_final.bam

##### END #####
echo "Post-alignment processing complete for $sample"
