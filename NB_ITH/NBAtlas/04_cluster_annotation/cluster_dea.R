# Author: Victoria Flanary
# Date: 250129
# Objective: Run DEA to identify cluster-specific marker genes.

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(presto)
library(tidyverse)
library(here)
library(future)

# Future parameters
options(future.globals.maxSize = 256 * 1024^3)  # 256 GB maximum RAM
plan("multisession", workers = 4)  # parallelize across 4 cores

# Set filepaths
data_dir <- here("NB_ITH", "NBAtlas", "data", "alldata")
results_dir <- here("NB_ITH", "NBAtlas", "04_cluster_annotation", "cluster_dea")

# Load data
seurat_obj <- readRDS(here(data_dir, "05_harmony_clust.rds"))

# JoinLayers for DEA
seurat_obj <- JoinLayers(seurat_obj)

# Set Ident to desired cluster size
seurat_obj <- SetIdent(seurat_obj, value = "RNA_snn_res.0.1")
table(seurat_obj@active.ident)

# Save current object
saveRDS(
  seurat_obj,
  here(data_dir, "06_dea_obj.rds")
)

# Clusters drastically range in size from several thousand to less than 100
# Want to balance groups to reduce bias toward larger clusters while 
# Preserving biological signals from these large clusters

# Subsample clusters so that no cluster has more than 1000 cells
seurat_subset <- subset(
  seurat_obj,
  cells = WhichCells(seurat_obj, downsample = 1000)
)
table(seurat_subset@active.ident)

# DEA
## Wilcoxon rank-sum (default)
marker_genes <- FindAllMarkers(
  seurat_subset,
  test.use = "wilcox",
  logfc.threshold = 2,
  return.thresh = 0.05,
  min.pct = 0.20,
  only.pos = TRUE,
  assay = "RNA"
)

## Filter and arrange DEGs
marker_genes_filt <- marker_genes |>
  filter(p_val_adj < 0.05) |>  # only keep degs that are significant after multiple hypothesis correction
  filter(pct.1 > 0.5) |>  # only keep degs expressed by most cells in the cluster
  mutate(specificity = pct.1 - pct.2) |>  # calculate difference in expression between cluster vs rest
  relocate(specificity, .before = p_val_adj) |>
  filter(specificity > 0.5) |> # only keep cluster-specific degs
  arrange(cluster, desc(avg_log2FC)) 

## Save the DEGs
write_csv(
  marker_genes,
  here(results_dir, "marker_genes.csv")
)

## Plot the top 25 cluster-specific DEGs
top25 <- marker_genes_filt |>
  group_by(cluster) |>
  slice_head(n = 25)

pdf(
  here(results_dir, "marker_genes_bar_plots.pdf"),
  height = 7, width = 9
)

par(mfrow = c(2, 4), mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  barplot(
    sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], FALSE),
    horiz = TRUE,
    las = 1,
    main = paste0(i, " vs. rest"),
    border = "white",
    yaxs = "i"
  )
  abline(v = c(0, 0.25), lty = c(1, 2))
}

dev.off()

# End of script
