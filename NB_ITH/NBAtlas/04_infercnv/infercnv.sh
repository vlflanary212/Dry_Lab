#!/bin/bash

##### SLURM #####
#SBATCH --job-name=infercnv
#SBATCH --partition=long, largemem-long
#SBATCH --time=149:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load Singularity

##### VARIABLES #####
sif="/home/flanary/sif_files/scrna-seq_1.2.sif"
script="/home/flanary/Projects/Dry_Lab/NB_ITH/NBAtlas/04_infercnv/infercnv.R"

##### COMMANDS #####
singularity exec $sif Rscript $script

##### END #####
echo "done"
