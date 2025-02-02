# Create a text file mapping cell line replicates to their appropriate fastq

# Load packages
library(here)

# Create mapping
## current file names 
cell_line_replicates <- c("IMR5_1", "IMR5_2", "SHEP1", "SHEP2", "SKNAS1",
                          "SKNAS2", "SY5Y_1", "SY5Y_2")
cell_line_replicates <- rep(cell_line_replicates, each = 2)
sample_num <- rep(paste0("S", 29:36), each = 2)
old_prefixes <- paste(cell_line_replicates, sample_num, sep = "_")
read <- rep(c("R1", "R2"), 8)
old_suffixes <- paste0(read, "_001.fastq.gz")
old_filenames <- paste(old_prefixes, old_suffixes, sep = "_")

## new filenames
cell_lines <- c("IMR5", "SHEP", "SKNAS", "SY5Y")
cell_lines <- rep(cell_lines, each = 4)
replicates <- rep(c("rep1", "rep2"), 4, each = 2)
new_prefixes <- paste(cell_lines, replicates, sep = "_")
new_suffixes <- paste0(read, ".fastq.gz")
new_filenames <- paste(new_prefixes, new_suffixes, sep = "_")
  
## create the mapping dataframe
mapping_df <- data.frame(old_filenames, new_filenames)

# Save the mapping to a text file
write.table(
  mapping_df,
  here("13cRA", "atac-seq", "treated", "Sen_241217", "00_rename_fastq", "mapping.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Also save replicates and cell lines in separate lists
write.table(
  unique(new_prefixes),
  here("13cRA", "atac-seq", "treated", "Sen_241217", "replicates.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  unique(cell_lines),
  here("13cRA", "atac-seq", "treated", "Sen_241217", "cell_lines.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
