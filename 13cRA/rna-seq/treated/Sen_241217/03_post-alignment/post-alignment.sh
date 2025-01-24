#!/bin/bash

##### SLURM #####
#SBATCH --job-name=samtools
#SBATCH --partition=medium
#SBATCH --time=24:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES #####
module load SAMtools/1.9-foss-2018b

###### ARRAY #######
samples="/home/flanary/Dry_Lab/13cRA/rna-seq/treated/Sen_241217/cell_lines.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$samples")

##### VARIABLES #####
bam_dir="/data/scratch/flanary/rna-seq/treated/Sen_241217/bam"
output_dir="/data/scratch/flanary/rna-seq/treated/Sen_241217/final_bam"

##### COMMANDS #####
# merge bam for technical replicates
samtools merge $bam_dir/"$sample"_merged.bam $bam_dir/"$sample"_rep1.Aligned.sortedByCoord.out.bam $bam_dir/"$sample"_rep2.Aligned.sortedByCoord.out.bam

# index
samtools index $bam_dir/"$sample"_merged.bam

# filter bam files
samtools view -b -F 4 -q 20 $bam_dir/"$sample"_merged.bam > $bam_dir/"$sample"_filt.bam

# remove duplicates
samtools rmdup -S $bam_dir/"$sample"_filt.bam $bam_dir/"$sample"_rmdup.bam

# index
samtools index $bam_dir/"$sample"_rmdup.bam

# move the final bam files to a separate directory
cp $bam_dir/"$sample"_rmdup.bam $output_dir/"$sample"_final.bam
cp $bam_dir/"$sample"_rmdup.bam.bai $output_dir/"$sample"_final.bam.bai

##### END #####
echo "Finish post-alignment processing with SAMtools"
