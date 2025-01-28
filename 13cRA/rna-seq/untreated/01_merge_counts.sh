#!/bin/bash

##### SLURM #####
#SBATCH --job-name=untx_counts
#SBATCH --partition=express
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load Singularity

##### VARIABLES #####
sif="/home/flanary/sif_files/scrna-seq_1.0.sif"
script="/home/flanary/Dry_Lab/13cRA/rna-seq/untreated/01_merge_counts.R"

##### COMMANDS #####
singularity exec $sif Rscript $script

##### END #####
echo "done"