#!/bin/bash

##### VARIABLES #####
fastq_dir="/data/scratch/flanary/atac-seq/untreated/Sen_240503/fastq"
mapping="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/Sen_240503/00_rename_fastq/mapping.txt"

##### COMMANDS #####
# Enter the fastq directory
cd "$fastq_dir"

# Read the mapping file line by line
while IFS=$'\t' read -r old_filename new_filename; do
    # Define the old and new filenames
    old_filename="${old_filename}"
    new_filename="${new_filename}"

    # Rename the file
    mv "$old_filename" "$new_filename"
    echo "Renamed $old_filename to $new_filename"
done < "$mapping"

##### END #####
echo "Renaming complete!"
