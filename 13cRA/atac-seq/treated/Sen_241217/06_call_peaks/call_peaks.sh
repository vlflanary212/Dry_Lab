#!/bin/bash

##### SLURM #####
#SBATCH --job-name=call_peaks
#SBATCH --partition=short
#SBATCH --time=11:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES #####
module load Anaconda3
conda activate atac-seq

##### VARIABLES #####
wd="/data/scratch/flanary/atac-seq/treated/Sen_241217"
bam_dir=$wd/"final_bam"
bigwig_dir=$wd/"bw"
peak_dir=$wd/"peaks"
samples="/home/flanary/Dry_Lab/13cRA/atac-seq/treated/Sen_241217/cell_lines.txt"

# Get sample name for this array task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samples")

##### COMMANDS #####
# Create bigwig file
bamCoverage -b $bam_dir/"$sample"_final.bam -o $bigwig_dir/"$sample".bw --binSize 200 --normalizeUsing None --effectiveGenomeSize 2913022398 

# Call peaks
macs2 callpeak --treatment $bam_dir/"$sample"_final.bam --name "$sample" --outdir $peak_dir --gsize hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits --pvalue 0.01
 
##### END #####
echo "done"
