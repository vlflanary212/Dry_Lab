# Author: Victoria Flanary
# Date: 250131
# Objective: Run infercnv on the NBAtlas

# Set seed
set.seed(42)

# Load libraries
library(Seurat)
library(infercnv)
library(tidyverse)
library(here)
library(future)

# Future parameters
options(future.globals.maxSize = 256 * 1024^3)  # 256 GB maximum RAM
plan("multisession", workers = 16)  # parallelize across 4 cores

# Set filepaths
wd <- here("NB_ITH", "NBAtlas", "04_infercnv")

# Load data
seurat_obj <- readRDS(
  here("NB_ITH", "NBAtlas", "data", "alldata", "01_filt_obj.rds")
)

# Extract the counts matrix
counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

# Prepare the annotations file
cell_annotations <- data.frame(
  cell_id = as.character(colnames(seurat_obj)),
  cell_type = as.character(seurat_obj@active.ident)
)

write.table(
  cell_annotations,
  here(wd, "inputs", "annotations.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Create the InferCNV object
all_cell_types <- unique(as.character(seurat_obj@active.ident))
normal_cell_types <- setdiff(all_cell_types, "Neuroendocrine")

infercnv_object <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = cell_annotations,
  gene_order_file = here(wd, "inputs", "gene_positions.txt"),
  delim = "\t",
  ref_group_names = normal_cell_types
) # Use all non-neuroendocrine cells as reference normal
# NBAtlas authors already showed that the neuroendocrine cells were tumor,
# And I confirmed their annotations through marker gene expression analysis

# Run InferCNV
options(scipen = 100)

infercnv_final <- infercnv::run(
  infercnv_obj = infercnv_object,
  analysis_mode = "subclusters",
  cluster_by_groups = TRUE,
  tumor_subcluster_pval = 0.05,
  tumor_subcluster_partition_method = "leiden", 
  leiden_resolution = 0.1, 
  cutoff = 0.1,  # set to 0.1 for 10X Genomics data
  out_dir = here(wd, "outputs"),
  denoise = TRUE,  # Enable denoising
  HMM = TRUE,  # Enable HMM for calling CNVs
  per_chr_hmm_subclusters = TRUE, # Run suclustering per chromosome over all cells
  HMM_report_by = "cell", # Get the HMM results for each cell, not subcluster
  num_threads = 16,
  resume_mode = FALSE  # Disable reloading from previous runs
)

# End of analysis
print("InferCNV on entire NBAtlas complete")