#!/bin/bash

##### SLURM #####
#SBATCH --job-name=align
#SBATCH --partition=medium
#SBATCH --time=49:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-8
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load BWA/0.7.17-foss-2018b
module load SAMtools/1.9-foss-2018b

##### ARRAY #####
sample_list="/home/flanary/Dry_Lab/13cRA/atac-seq/treated/Sen_241217/replicates.txt"
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$sample_list")

##### VARIABLES #####
# directories
file_dir="/data/scratch/flanary/atac-seq/treated/Sen_241217"
fastq_dir=$file_dir/"trimmed_fastq"  #fastq storage directory
bam_dir=$file_dir/"bam"  #bam storage directory

# reference genome
ref_genome="/data/project/sen-lab/genome/hg38/bwa/Homo_sapiens_assembly38.fasta"

# fastq files (2/sample due to paired-end sequencing)
fastq1="$fastq_dir/$sample"_R1_val_1.fq.gz
fastq2="$fastq_dir/$sample"_R2_val_2.fq.gz

##### COMMANDS #####
# map reads and make bam files
bwa mem -M -t 8 $ref_genome $fastq1 $fastq2 | samtools view -S -b -h -F 4 -q 20 > $bam_dir/"$sample"_mapped.bam

# sort bam files by coordinate
samtools sort -o $bam_dir/"$sample"_sort.bam $bam_dir/"$sample"_mapped.bam

##### END #####
echo "Alignment complete for sample: $sample"
