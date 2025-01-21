#!/bin/bash

##### SLURM #####
#SBATCH --job-name=align
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

##### PACKAGES #####
module load STAR/2.7.3a-GCC-6.4.0-2.28

###### ARRAY #######
samples="/home/flanary/Dry_Lab/13cRA/rna-seq/untreated/GSE28875/SRR_Acc_List.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$samples")

##### VARIABLES #####
wd="/data/scratch/flanary/rna-seq/untreated/GSE28875"
fastq_dir=$wd/"fastq"
output_dir=$wd/"bam"

##### COMMANDS #####
echo "Aligning sample: $sample"

# make output directory
output=$output_dir/$sample
mkdir -p $output

# STAR indices
star_index="/data/project/sen-lab/genome/hg38/star_index"
echo "Using STAR index: $star_index"

# align
echo "FASTQ files: $fastq_dir/$sample.fastq"

STAR --runThreadN 8 --genomeDir $star_index \
--readFilesIn $fastq_dir/"$sample".fastq \
--outFileNamePrefix $output/$sample. \
--outSAMtype BAM SortedByCoordinate

echo "Aligned files for $sample are in: $output"

##### END #####
echo "STAR alignment complete"