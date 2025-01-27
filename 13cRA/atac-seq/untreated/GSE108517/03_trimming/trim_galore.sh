#!/bin/bash

##### SLURM #####
#SBATCH --job-name=trim
#SBATCH --partition=short
#SBATCH --time=11:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-2
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES #####
module load Anaconda3
conda activate atac-seq

##### VARIABLES #####
wd="/data/scratch/flanary/atac-seq/untreated/GSE108517"
fastq_dir=$wd/"fastq"
output_dir=$wd/"trimmed_fastq"
samples="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/GSE108517/SRR_Acc_List.txt"

# Get sample name for this array task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samples")

# Input FASTQ files
r1="$fastq_dir/${sample}_1.fastq"
r2="$fastq_dir/${sample}_2.fastq"

##### COMMANDS #####
echo "Processing sample: $sample"

# Trim adapters and low quality reads (Phred score < 20 or length < 30 bp)
trim_galore --paired --cores 8 --output_dir "$output_dir" \
--quality 20 --length 30 "$r1" "$r2"


if [[ $? -eq 0 ]]; then
  echo "Trimming completed successfully for sample: $sample"
else
  echo "Error during trimming for sample: $sample"
  exit 1
fi

##### END #####
conda deactivate
