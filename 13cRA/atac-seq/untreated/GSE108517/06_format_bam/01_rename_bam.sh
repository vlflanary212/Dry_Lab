#!/bin/bash

##### VARIABLES #####
bam_dir="/data/scratch/flanary/atac-seq/untreated/GSE108517/bam"
mapping="/home/flanary/Dry_Lab/13cRA/atac-seq/untreated/GSE108517/06_format_bam/bam_filename_mapping.txt"

# Enter the bam directory
cd "$bam_dir"

# Read the mapping file line by line
while IFS=$'\t' read -r old_prefix new_prefix; do
    # Define the old and new filenames
    old_bam="${old_prefix}_mapped.bam"
    new_bam="${new_prefix}_mapped.bam"

    # Rename the file
    mv "$old_bam" "$new_bam"
    echo "Renamed $old_bam to $new_bam"
done < "$mapping"

echo "Renaming complete!"
