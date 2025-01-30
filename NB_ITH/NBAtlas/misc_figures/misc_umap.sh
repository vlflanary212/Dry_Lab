#!/bin/bash

##### SLURM #####
#SBATCH --job-name=alldata
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load Singularity

##### VARIABLES #####
sif="/home/flanary/sif_files/scrna-seq_1.0.sif"
script="/home/flanary/Dry_Lab/NB_ITH/NBAtlas/misc_figures/misc_umap.R"

##### COMMANDS #####
singularity exec $sif Rscript $script

##### END #####
echo "done"