# Author: Victoria Flanary
# Date: 240122
# Objective: Visualize the number of tumors belonging to each metadata
# variable in the NBAtlas

# Set-up
## Set a seed
set.seed(42)

## Load  packages
library(Seurat)
library(tidyverse)
library(here)

## Load data
seurat_obj <- readRDS(
  here("NB_ITH", "NBAtlas", "data", "alldata", "00_init_obj.rds")
)

## Set results directory
results_dir <- here("NB_ITH", "NBAtlas", "01_dataset_eda", "dataset_overview", "bar_plots")

# Extract metadata
metadata <- seurat_obj@meta.data

#  Note the types of metadata available
colnames(metadata)
# [1] "orig.ident"         "nCount_RNA"         "nFeature_RNA"       "percent_mito"      
# [5] "Study"              "Assay"              "Platform"           "Sample"            
# [9] "Patient_No"         "Timepoint"          "INSS_stage"         "MYCN_amplification"
# [13] "Gender"             "Risk"               "Cell_condition"     "Cell_type" 

# Continuous technical: nCount_RNA, nFeature_RNA, percent_mito
# Discrete technical: Study, Assay, Platform, Cell_condition
# Discrete biological: orig.ident, Sample, Patient_No, Timepoint, INSS_stage,
# MYCN_amplification, Gender, Risk, Cell_type

# orig.ident seems to refer to original cell type annotations given by the 
# NBAtlas authors instead of patient or tumor sample
identical(metadata$orig.ident, metadata$Cell_type)  # TRUE

# Isolate metadata of interest for plotting
## Want to know how many tumors belong to each metadata variable
## For Cell_type, how many tumors contain that cell type
vars <- c("Study", "Assay", "Platform", "Sample", "Patient_No", "Gender",
          "MYCN_amplification", "Risk", "INSS_stage", "Timepoint", "Cell_type")

metadata_filt <- select(metadata, all_of(vars))

# Plot the number of tumors per metadata variable as bar plots
## Plot in a single pdf
pdf(
    here(results_dir, paste0("Num_Tumors_by_Metadata.pdf")),
    height = 4, width = 6
  )

for (i in vars) {
  # Select the relevant columns and summarize data
  metadata_summary <- metadata_filt[, c("Sample", i)] |>  # Select Sample and current metadata column
    distinct() |>
    group_by(.data[[i]]) |>   # Group by current metadata variable
    summarise(n = n(), .groups = "drop")  # Count the number of Samples per metadata value
  
  # Assign consistent column names
  colnames(metadata_summary) <- c("MetadataValue", "n")
  
  # Adjust y-axis limit
  max_y <- max(metadata_summary$n, na.rm = TRUE)
  adjusted_max_y <- max_y * 1.1
  
  # Create and print the plot
  print(
    ggplot(metadata_summary, aes(x = MetadataValue, y = n, fill = MetadataValue)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = n), vjust = -0.5, size = 4) +
      ylim(0, adjusted_max_y) +
      labs(
        title = paste("Number of Samples by", i),
        x = NULL,
        y = "Number of Samples"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        legend.position = "none",  # Remove the plot legend
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text
      )
  )
}

dev.off()

## Plot as individual files
for (i in vars) {
  pdf(
    here(results_dir, paste0("Num_Tumors_by_", i, ".pdf")),
    height = 4, width = 6
  )
  
  # Select the relevant columns and summarize data
  metadata_summary <- metadata_filt[, c("Sample", i)] |>  # Select Sample and current metadata column
    distinct() |>
    group_by(.data[[i]]) |>   # Group by current metadata variable
    summarise(n = n(), .groups = "drop")  # Count the number of Samples per metadata value
  
  # Assign consistent column names
  colnames(metadata_summary) <- c("metadata_value", "n")
  
  # Adjust y-axis limit
  max_y <- max(metadata_summary$n, na.rm = TRUE)
  adjusted_max_y <- max_y * 1.1
  
  # Create and print the plot
  print(
    ggplot(metadata_summary, aes(x = metadata_value, y = n, fill = metadata_value)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = n), vjust = -0.5, size = 4) +
      ylim(0, adjusted_max_y) +
      labs(
        title = paste("Number of Samples by", i),
        x = NULL,
        y = "Number of Samples"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),  # Center the plot title
        legend.position = "none",  # Remove the plot legend
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text
      )
  )
  
  dev.off()
}

# End of script