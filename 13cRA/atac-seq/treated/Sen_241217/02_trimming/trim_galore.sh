#!/bin/bash

##### SLURM #####
#SBATCH --job-name=trim
#SBATCH --partition=short
#SBATCH --time=11:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-8
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES #####
module load Anaconda3
conda activate atac-seq

##### VARIABLES #####
wd="/data/scratch/flanary/atac-seq/treated/Sen_241217"
fastq_dir=$wd/"fastq"
output_dir=$wd/"trimmed_fastq"
samples="/home/flanary/Dry_Lab/13cRA/atac-seq/treated/Sen_241217/replicates.txt"

# Get sample name for this array task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samples")

# Input FASTQ files
r1="$fastq_dir/${sample}_R1.fastq.gz"
r2="$fastq_dir/${sample}_R2.fastq.gz"

##### COMMANDS #####
echo "Processing sample: $sample"

# Trim adapters and low quality reads (Phred score < 20 or length < 30 bp)
trim_galore --paired --cores 8 --output_dir $output_dir $r1 $r2

##### END #####
conda deactivate
