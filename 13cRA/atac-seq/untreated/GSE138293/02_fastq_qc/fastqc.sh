#!/bin/bash

##### SLURM #####
#SBATCH --job-name=fastqc
#SBATCH --partition=short
#SBATCH --time=11:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-55
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### VARIABLES #####
fastq_dir="/data/scratch/flanary/atac-seq/untreated/GSE138293/fastq"
fastqc_output="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/GSE138293/02_fastq_qc/fastqc_output"

##### PACKAGES #####
module load FastQC/0.11.9-Java-11

##### COMMANDS #####
# Get a list of FASTQ files
fastq_files=(${fastq_dir}/*.fastq)

# Get a single file for the current task
fastq_file=${fastq_files[$SLURM_ARRAY_TASK_ID]}
fastqc -o $fastqc_output "$fastq_file"

##### END #####
echo "done"