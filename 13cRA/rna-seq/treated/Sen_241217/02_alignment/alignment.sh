#!/bin/bash

##### SLURM #####
#SBATCH --job-name=star
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
module load STAR/2.7.3a-GCC-6.4.0-2.28

###### ARRAY #######
samples="/home/flanary/Dry_Lab/13cRA/rna-seq/treated/Sen_241217/replicates.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$samples")

##### VARIABLES #####
wd="/data/scratch/flanary/rna-seq/treated/Sen_241217"
fastq_dir=$wd/"fastq"
output_dir=$wd/"bam"

##### COMMANDS #####
echo "Processing sample $sample"

# STAR indices
star_index="/data/project/sen-lab/genome/hg38/star_index"
echo "Using STAR index: $star_index"

# align
echo "FASTQ files: $fastq_dir/"$sample"_R1.fastq.gz, $fastq_dir/"$sample"_R2.fastq.gz"

STAR --runThreadN 8 --genomeDir $star_index \
     --readFilesIn $fastq_dir/"$sample"_R1.fastq.gz $fastq_dir/"$sample"_R2.fastq.gz \
     --outFileNamePrefix $output_dir/$sample. \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat

echo "Aligned files for $sample are in: $output_dir"

##### END #####
echo "STAR alignment complete"
