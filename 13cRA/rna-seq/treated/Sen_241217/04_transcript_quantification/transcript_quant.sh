#!/bin/bash

##### SLURM #####
#SBATCH --job-name=tx_quant
#SBATCH --partition=express
#SBATCH --time=01:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES #####
module load Anaconda3
conda activate rna-seq

##### VARIABLES #####
wd="/data/scratch/flanary/rna-seq/treated/Sen_241217"
bam_dir=$wd/"final_bam"
output_dir=$wd/"counts"

###### ARRAY #######
samples="/home/flanary/Dry_Lab/13cRA/rna-seq/treated/Sen_241217/cell_lines.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$samples")

##### COMMANDS #####
# make output directory
sample_output=$output_dir/$sample
mkdir -p $sample_output

# gtf file
gtf="/data/project/sen-lab/genome/hg38/gencode.v22.annotation.gtf"

# get counts
featureCounts -p -T 8 -t exon -g gene_id -a $gtf -o $sample_output/raw_counts.txt $bam_dir/"$sample"_final.bam

# edit counts
cat $sample_output/raw_counts.txt| cut -f1,7- | sed 1d > $sample_output/counts.txt

##### END #####
echo "Transcript quantification complete for $sample"
conda deactivate