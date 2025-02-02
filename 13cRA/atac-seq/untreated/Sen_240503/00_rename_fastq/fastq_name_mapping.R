# Create a text file mapping cell line replicates to their appropriate fastq

# Load packages
library(here)

# Create mapping
## Current file names 
cell_lines <- c("GIMEN", "SY5Y", "SHEP", "IMR5")
cell_lines <- rep(cell_lines, each = 2)

sample_suffixes <- paste0("S", 1:4)
sample_suffixes <- rep(sample_suffixes, each = 2)

samples <- paste(cell_lines, sample_suffixes, sep = "_")

reads <- c("L001_R1", "L001_R2")
reads <- rep(read_prefixes, 4)

file_suffixes <- rep("001.fastq.gz", 8)

old_filenames <- paste(samples, reads, file_suffixes, sep = "_")

## New file names
cell_lines <- c("GIMEN", "SH-SY5Y", "SHEP-1", "IMR-5")
cell_lines <- rep(cell_lines, each = 2)

file_suffixes <- paste0(c("R1", "R2"), ".fastq.gz")
file_suffixes <- rep(file_suffixes, 4)

new_filenames <- paste(cell_lines, file_suffixes, sep = "_")

## create the mapping dataframe
mapping_df <- data.frame(old_filenames, new_filenames)

# Save the mapping to a text file
write.table(
  mapping_df,
  here("13cRA", "atac-seq", "untreated", "Sen_240503", "00_rename_fastq", "mapping.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Also save renamed cell lines in a separate list
write.table(
  unique(cell_lines),
  here("13cRA", "atac-seq", "untreated", "Sen_240503", "cell_lines.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
