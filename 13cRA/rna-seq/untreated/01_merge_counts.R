# Author: Victoria Flanary
# Date: 240127
# Objective: Generate a combined counts matrix for bulk rna-seq data of 
# untreated neuroblastoma cell lines

## Set-up
library(data.table)
library(tidyverse)
library(here)

# GSE28875
## Only has a single sample - retrieve counts directly
hncc_file <- "/data/scratch/flanary/rna-seq/untreated/GSE28875/SRR191362/counts.txt"
hncc_counts <- fread(hncc_file, header = TRUE) |> as.data.frame()
colnames(hncc_counts) <- c("gene_id", "SRR191362")

# GSE89413
## Define filepaths
nb_dir <- "/data/scratch/flanary/rna-seq/untreated/GSE89413"
sample_list <- "/home/flanary/Dry_Lab/13cRA/rna-seq/untreated/GSE89413/SRR_Acc_List.txt"

## Read in sample list
samples <- readLines(sample_list)

## Retrieve counts for each sample in a list
counts_list <- list()

for (i in samples) {
  counts_file <- file.path(nb_dir, i, "counts.txt")
  counts <- fread(counts_file, header = TRUE)
  colnames(counts) <- c("gene_id", i)
  counts_list[[i]] <- counts
}

## Merge the counts list
nb_counts <- purrr::reduce(counts_list, full_join, by = "gene_id") |>
  as.data.frame()

# Merge the counts 
merged_counts <- full_join(nb_counts, hncc_counts, by = "gene_id")

# Move gene ids to rownames
rownames(merged_counts) <- merged_counts$gene_id
merged_counts <- merged_counts[, colnames(merged_counts) != "gene_id"]

# Save the merged counts
write.table(
  merged_counts,
  here("/home/flanary/Dry_Lab/13cRA/rna-seq/untreated/raw_counts.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# End of script