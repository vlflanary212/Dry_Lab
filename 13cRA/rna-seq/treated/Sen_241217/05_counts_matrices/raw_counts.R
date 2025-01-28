# Author: Victoria Flanary
# Date: 240127
# Objective: Generate a combined counts matrix for bulk rna-seq data of 
# treated neuroblastoma cell lines

## Set-up
library(data.table)
library(tidyverse)
library(here)

# GSE89413
## Define filepaths
file_dir <- "/Users/victoriaflanary/Library/CloudStorage/Box-Box//Data/rna-seq/treated/Sen_241217/counts"
sample_list <- here("13cRA", "rna-seq", "treated", "Sen_241217", "cell_lines.txt")

## Read in sample list
samples <- readLines(sample_list)

## Retrieve counts for each sample in a list
counts_list <- list()

for (i in samples) {
  counts_file <- file.path(file_dir, i, "counts.txt")
  counts <- fread(counts_file, header = TRUE)
  colnames(counts) <- c("gene_id", i)
  counts_list[[i]] <- counts
}

## Merge the counts list
merged_counts <- purrr::reduce(counts_list, full_join, by = "gene_id") |>
  as.data.frame()

# Move gene ids to rownames
rownames(merged_counts) <- merged_counts$gene_id
merged_counts <- merged_counts[, colnames(merged_counts) != "gene_id"]

# Save the merged counts
write.table(
  merged_counts,
  here("13cRA", "rna-seq", "treated", "Sen_241217", "raw_counts.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# End of script