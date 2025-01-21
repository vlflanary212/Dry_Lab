# Generate metadata dataframes for the untreated ATAC-seq data

# Load packages
library(tidyverse)
library(here)

# Load SRA Run Tables
wd <- here("13cRA", "atac-seq", "untreated")
hncc <- read_csv(here(wd, "GSE108517", "SraRunTable.csv"))
nb <- read_csv(here(wd, "GSE138293", "SraRunTable.csv"))

# Check available metadata variables
intersect(colnames(hncc), colnames(nb))
setdiff(colnames(hncc), colnames(nb))

# Subset for variables of use of interest
hncc_filt <- dplyr::select(hncc, run = "Run", sample = "Sample Name", 
                      cell_line = "cell_line", assay = "Assay Type", 
                      instrument = "Instrument", layout = "LibraryLayout",
                      platform = "Platform", study = "SRA Study"
                      )
hncc_filt$cell_line <- "hNCC"

nb_filt <- dplyr::select(nb, run = "Run", sample = "Sample Name", 
                    cell_line = "cell_line", assay = "Assay Type", 
                    instrument = "Instrument", layout = "LibraryLayout",
                    platform = "Platform", study = "SRA Study"
                      )
# Add replicate column
hncc_rep <- mutate(hncc_filt, replicate = paste0("rep", c(1,2)))
hncc_rep <- mutate(hncc_rep, rep_name = paste0(cell_line, "_", replicate))
hncc_final <- relocate(hncc_rep, replicate, rep_name, .after = cell_line)

rep_num <- rep(c(1,2), 14)
nb_rep <- mutate(nb_filt, replicate = paste0("rep", rep_num))
nb_rep <- mutate(nb_rep, rep_name = paste0(cell_line, "_", replicate))
nb_final <- relocate(nb_rep, replicate, rep_name, .after = cell_line)

# Save the data frames
write_csv(
  hncc_final,
  here(wd, "GSE108517", "metadata_filt.csv")
)

write_csv(
  nb_final,
  here(wd, "GSE138293", "metadata_filt.csv")
)

# Save run and rep_names separately to rename bam files for merging
hncc_mapping <- data.frame(
  run = hncc_final$run,
  name = hncc_final$rep_name
)

write.table(
  hncc_mapping,
  here(wd, "GSE108517", "06_format_bam", "bam_filename_mapping.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  unique(hncc_final$cell_line),
  here(wd, "GSE108517", "06_format_bam", "cell_lines.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  hncc_mapping$name,
  here(wd, "GSE108517", "06_format_bam", "sample_replicates.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

nb_mapping <- data.frame(
  run = nb_final$run,
  name = nb_final$rep_name
)

write.table(
  nb_mapping,
  here(wd, "GSE138293", "06_format_bam", "bam_filename_mapping.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  unique(nb_final$cell_line),
  here(wd, "GSE138293", "06_format_bam", "cell_lines.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  nb_mapping$name,
  here(wd, "GSE138293", "06_format_bam", "sample_replicates.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
