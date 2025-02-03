#!/bin/bash

##### SLURM #####
#SBATCH --job-name=bwa
#SBATCH --partition=medium
#SBATCH --time=49:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load BWA/0.7.17-foss-2018b
module load SAMtools/1.9-foss-2018b

##### VARIABLES #####
# directories
file_dir="/data/scratch/flanary/atac-seq/untreated/Sen_240503"
fastq_dir=$file_dir/"fastq"  #fastq storage directory
bam_dir=$file_dir/"bam"  #bam storage directory

# reference genome
ref_genome="/data/project/sen-lab/genome/hg38/bwa/Homo_sapiens_assembly38.fasta"

# sample - re-running GIMEN separately as single-ended since R2 file is a duplicate of R1
sample="GIMEN"

# fastq file
fastq="$fastq_dir/$sample"_R1.fastq.gz

##### COMMANDS #####
# map reads and make bam files
bwa mem -M -t 8 $ref_genome $fastq | samtools view -S -b -h -F 4 -q 20 > $bam_dir/"$sample"_mapped.bam

##### END #####
echo "Alignment complete for sample: $sample"
