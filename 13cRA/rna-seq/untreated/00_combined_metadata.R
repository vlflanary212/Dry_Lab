# Perform EDA using the metadata of untreated neuroblastoma lines acquired
# from SRAToolkit

# Load packages
library(tidyverse)
library(here)

# Load metadata
wd <- here("13cRA", "rna-seq", "untreated")
hncc <- read_csv(here(wd, "GSE28875", "SraRunTable.csv"))
nb <- read_csv(here(wd, "GSE89413", "SraRunTable.csv"))

# Select columns of interest
intersect(colnames(hncc), colnames(nb))
setdiff(colnames(hncc), colnames(nb))

# Filter for metadata of interest
metadata <- c(study = "SRA Study", run = "Run", 
              sample = "Sample Name", source = "source_name",
              organism = "Organism", assay = "Assay Type", 
              platform = "Platform", instrument = "Instrument",
              layout = "LibraryLayout", selection = "LibrarySelection",
              data_type = "LibrarySource", center = "Center Name",
              experiment = "Experiment")

hncc_filt <- dplyr::select(hncc, all_of(metadata))
nb_filt <- dplyr::select(nb, all_of(metadata))

# Add cell line info
## Read in samples
nb_lines <- read.table(
  here(wd, "GSE89413", "cell_lines_by_sample.txt"),
  header = TRUE
)

## Ensure both objects are in the same order
identical(nb_lines$sample, nb_filt$sample)  # TRUE

## Add the samples
nb_df <- left_join(nb_filt, nb_lines, by = "sample")

## Add cell_line col to hncc df too
hncc_df <- hncc_filt
hncc_df$cell_line <- "hNCC"

# Combine dataframes
metadata_df <- rbind(nb_df, hncc_df)

# Move the cell_line column
metadata_df <- relocate(metadata_df, cell_line, .after = sample)

# Change the formatting for the tissue source column
metadata_final <- mutate(
  metadata_df, 
  source = if_else(source != "neuroblastoma-derived cell line", "normal_tissue", source)
)

# Save the metadata
write_csv(
  metadata_final,
  here(wd, "combined_metadata.csv")
)
