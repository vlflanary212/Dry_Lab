#!/bin/bash

##### SLURM #####
#SBATCH --job-name=fastq
#SBATCH --partition=express
#SBATCH --time=01:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load SRA-Toolkit/3.0.0-centos_linux64

##### VARIABLES #####
wd="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/GSE108517"
output_dir="/data/scratch/flanary/atac-seq/untreated/GSE108517/fastq"

##### ARRAY #####
sample_list=$wd/"SRR_Acc_List.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$sample_list")
echo "sample:$sample"

##### COMMANDS #####
# download fastq files for specified samples from the SRR_Acc_List.txt
fasterq-dump $sample --outdir $output_dir

##### END #####
echo "done"
