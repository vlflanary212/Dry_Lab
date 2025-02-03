#!/bin/bash

##### SLURM #####
#SBATCH --job-name=multiqc
#SBATCH --partition=express
#SBATCH --time=01:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### VARIABLES #####
wd="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/GSE138293/04_post-trim_qc"
fastqc_output=$wd/"fastqc_output"
multiqc_output=$wd/"multiqc_output"

##### PACKAGES #####
module load Anaconda3
conda activate atac-seq

##### COMMANDS #####
multiqc -v $fastqc_output -o $multiqc_output

##### END #####
echo "done"
