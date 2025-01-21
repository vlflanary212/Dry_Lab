#!/bin/bash

##### SLURM #####
#SBATCH --job-name=tx_quant
#SBATCH --partition=express
#SBATCH --time=01:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-40
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES #####
module load Anaconda3
conda activate featurecounts

##### VARIABLES #####
bam_dir="/data/scratch/flanary/rna-seq/untreated/GSE89413/bam"
output_dir="/data/scratch/flanary/rna-seq/untreated/GSE89413/counts"

###### ARRAY #######
samples="/home/flanary/Dry_Lab/13cRA/rna-seq/untreated/GSE89413/SRR_Acc_List.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$samples")

##### COMMANDS #####
echo "Processing sample: $sample"

# make output directory
sample_output=$output_dir/$sample
mkdir -p $sample_output

# gtf file
gtf="/data/project/sen-lab/genome/hg38/gencode.v22.annotation.gtf"

# get counts
featureCounts -p -T 8 -t exon -g gene_id -a $gtf -o $sample_output/raw_counts.txt $bam_dir/$sample/$sample.rmdup.bam

# edit counts
cat $sample_output/raw_counts.txt| cut -f1,7- | sed 1d > $sample_output/counts.txt

##### END #####
echo "Transcript quantification complete for sample: $sample"
conda deactivate
