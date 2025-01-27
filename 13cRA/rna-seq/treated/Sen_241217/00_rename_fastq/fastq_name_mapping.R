# Create a text file mapping cell line replicates to their appropriate fastq

# Load packages
library(here)

# Create mapping
## current file names
prefixes <- paste0(paste0("AS", 1:10), paste0("_S", 19:28))
prefixes <- rep(prefixes, each = 2)
old_suffixes <- paste0(rep(c("R1", "R2"), 10), "_001.fastq.gz")
old_filenames <- paste(prefixes, old_suffixes, sep = "_")

## samples 1-8 are for the 13cRA project
ra_lines <- rep(c("IMR5", "SHEP", "SKNAS", "SY5Y"), each = 2)
replicates <- rep(c("rep1", "rep2"), 4)
ra_replicates <- paste0(ra_lines, "_", replicates)

## samples 9-10 are for the FOXJ3 project
foxj3_samples <- c("FOXJ3_siRNA", "BE2C_neg_control")

## make a sample vector with the desired suffixes
samples <- c(ra_replicates, foxj3_samples)
samples <- rep(samples, each = 2)
new_suffixes <- paste0(rep(c("R1", "R2"), 10), ".fastq.gz")
new_filenames <- paste(samples, new_suffixes, sep = "_")

## create the mapping dataframe
mapping_df <- data.frame(old_filenames, new_filenames)

# Save the mapping to a text file
write.table(
  mapping_df,
  here("13cRA", "rna-seq", "treated", "Sen_241217", "00_rename_fastq", "mapping.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Also save ra-treated samples and cell lines in separate lists
write.table(
  ra_replicates,
  here("13cRA", "rna-seq", "treated", "Sen_241217", "replicates.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  unique(ra_lines),
  here("13cRA", "rna-seq", "treated", "Sen_241217", "cell_lines.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
