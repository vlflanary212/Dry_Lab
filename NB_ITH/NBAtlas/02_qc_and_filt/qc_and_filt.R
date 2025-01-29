# Author: Victoria Flanary
# Date: 250129
# Objective: Run QC metrics on the NBAtlas data. Filter if needed.

# Set-up
## Load  packages
library(Seurat)
library(tidyverse)
library(here)

## Load the original object
seurat_obj <- readRDS(
  here("NB_ITH", "NBAtlas", "data", "alldata", "00_init_obj.rds")
)

## Set results directory
results_dir <- here("NB_ITH", "NBAtlas", "02_qc_and_filt")

# Split the object by batch
seurat_obj@meta.data <- mutate(
  seurat_obj@meta.data,
  Batch = paste(Study, Assay, Platform, Cell_condition, sep = "_")
)

## Save batches in a text file
write.table(
  unique(seurat_obj$Batch),
  here("NB_ITH", "NBAtlas", "batches.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Calculate QC metrics
seurat_obj <- seurat_obj |>
  PercentageFeatureSet("^MT-", col.name = "percent_mito") |>
  PercentageFeatureSet("^RP[SL]", col.name = "percent_ribo") |>
  PercentageFeatureSet("^HB[^(P|E|S)]", col.name = "percent_hb")

# Visualize QC metrics across studies
## Violin plots
features <- c("nCount_RNA", "nFeature_RNA", "percent_mito", 
              "percent_ribo", "percent_hb")
pdf(
  here(results_dir, "unfilt", "qc_violin_plots.pdf"),
  height = 8, width = 8
)

for (i in features) {
  print(VlnPlot(seurat_obj, group.by = "Batch", features = i, raster = FALSE) +
          NoLegend())
}

dev.off()

## Ridge plots
pdf(
  here(results_dir, "unfilt", "qc_ridge_plots.pdf"),
  height = 6, width = 12
)

for (i in features) {
  print(
    RidgePlot(seurat_obj, group.by = "Batch", features = i) +
      theme(axis.title.y = element_blank()) +
      NoLegend()
  )
}

dev.off()

## Scatter plots
pdf(
  here(results_dir, "unfilt", "nCount_nFeature_scatter.pdf"),
  height = 4, width = 6
)

print(
  FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                 group.by = "Batch", shuffle = TRUE, raster = FALSE) +
    NoLegend()
)

dev.off()


pdf(
  here(results_dir, "unfilt", "nCount_pct_mito_scatter.pdf"),
  height = 4, width = 6
)

print(
  FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent_mito",
                 group.by = "Batch", shuffle = TRUE, raster = FALSE) +
    NoLegend()
)

dev.off()

## Box plot
### Extract the raw counts matrix from the Seurat object
counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

### Normalize the counts to percent counts per cell
counts@x <- counts@x / rep.int(colSums(counts), diff(counts@p)) * 100

### Get the top 20 most expressed genes
most_expressed <- order(Matrix::rowSums(counts), decreasing = TRUE)[20:1]

### Get the gene names for the most expressed genes
gene_names <- rownames(counts)[most_expressed]

### Convert to matrix while preserving row names (gene names)
counts_matrix <- as.matrix(t(counts[most_expressed, ]))

# Generate the box plot with correct labels
pdf(
  here(results_dir, "unfilt", "top_expressed_genes_box.pdf"),
  height = 6, width = 9
)

par(oma = c(0, 1, 0, 0))

print(
  boxplot(
    counts_matrix,
    cex = 0.1, las = 1, xlab = "Percent counts per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
    names = gene_names # Set the gene names on the y-axis
  )
)

par(oma = c(0, 1, 0, 0))

dev.off()

## Histogram of qc metrics for all studies combined
pdf(
  here(results_dir, "unfilt", "qc_hist.pdf"),
  height = 4, width = 6
)

for (i in features) {
  print(
    ggplot(seurat_obj@meta.data, aes(x = .data[[i]])) +
      geom_histogram() +
      labs(
        title = paste0(i, " Distribution"),
        x = i,
        y = "Frequency"
      )
  ) 
}

dev.off()

# Filtering
## Rm cells with outlier values for nCount_RNA and nFeature_RNA
## Likely doublets missed by doublet detection algorithms
seurat_filt <- subset(seurat_obj, subset = nCount_RNA < 30000 & nFeature_RNA < 6000)

## Rm cells with percent mitochondrial genes > 20%
## All cells with percent_mito above this point also have low nCount_RNA
selected_mito <- WhichCells(seurat_filt, expression = percent_mito < 20)
seurat_filt <- subset(seurat_filt, cells = selected_mito)

## Rm cells with percent_ribo > 60%
## Very few cells, but likely technical artifacts
selected_ribo <- WhichCells(seurat_filt, expression = percent_ribo < 60)
seurat_filt <- subset(seurat_filt, cells = selected_ribo)

## Rm cells with any expression of hemoglobin genes (RBC contamination)
selected_hb <- WhichCells(seurat_filt, expression = percent_hb == 0)
seurat_filt <- subset(seurat_filt, cells = selected_hb)

## Rm MALAT1 - disproportionately expressed compared to other genes
seurat_filt <- seurat_filt[!grepl("MALAT1", rownames(seurat_filt)), ]

## Rm mitochondrial and ribosomal genes
## Abundance of these transcripts masks interesting biological variability
seurat_filt <- seurat_filt[!grepl("^MT-", rownames(seurat_filt)), ]
seurat_filt <- seurat_filt[!grepl("^RP[SL]", rownames(seurat_filt)), ]

## Record difference in data size before and after filtering
dim(seurat_obj)  # 42413 362991
dim(seurat_filt)  # 42412 310165

# Re-plot QC metrics
## Violin plots
pdf(
  here(results_dir, "filt", "qc_violin_plots.pdf"),
  height = 8, width = 8
)

for (i in features) {
  print(VlnPlot(seurat_filt, group.by = "Batch", features = i, raster = FALSE) +
          NoLegend())
}

dev.off()

## Ridge plots
pdf(
  here(results_dir, "filt", "qc_ridge_plots.pdf"),
  height = 6, width = 12
)

for (i in features) {
  print(
    RidgePlot(seurat_filt, group.by = "Batch", features = i) +
      theme(axis.title.y = element_blank()) +
      NoLegend()
  )
}

dev.off()

## Scatter plots
pdf(
  here(results_dir, "filt", "nCount_nFeature_scatter.pdf"),
  height = 4, width = 6
)

print(
  FeatureScatter(seurat_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                 group.by = "Batch", shuffle = TRUE, raster = FALSE) +
    NoLegend()
)

dev.off()


pdf(
  here(results_dir, "filt", "nCount_pct_mito_scatter.pdf"),
  height = 4, width = 6
)

print(
  FeatureScatter(seurat_filt, feature1 = "nCount_RNA", feature2 = "percent_mito",
                 group.by = "Batch", shuffle = TRUE, raster = FALSE) +
    NoLegend()
)

dev.off()

## Box plot
### Extract the raw counts matrix from the Seurat object
counts <- GetAssayData(seurat_filt, assay = "RNA", layer = "counts")

### Normalize the counts to percent counts per cell
counts@x <- counts@x / rep.int(colSums(counts), diff(counts@p)) * 100

### Get the top 20 most expressed genes
most_expressed <- order(Matrix::rowSums(counts), decreasing = TRUE)[20:1]

### Get the gene names for the most expressed genes
gene_names <- rownames(counts)[most_expressed]

### Convert to matrix while preserving row names (gene names)
counts_matrix <- as.matrix(t(counts[most_expressed, ]))

# Generate the box plot with correct labels
pdf(
  here(results_dir, "filt", "top_expressed_genes_box.pdf"),
  height = 6, width = 9
)

par(oma = c(0, 1, 0, 0))

print(
  boxplot(
    counts_matrix,
    cex = 0.1, las = 1, xlab = "Percent counts per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
    names = gene_names # Set the gene names on the y-axis
  )
)

par(oma = c(0, 1, 0, 0))

dev.off()

## Histogram of qc metrics for all studies combined
pdf(
  here(results_dir, "filt", "qc_hist.pdf"),
  height = 4, width = 6
)

for (i in features) {
  print(
    ggplot(seurat_filt@meta.data, aes(x = .data[[i]])) +
      geom_histogram() +
      labs(
        title = paste0(i, " Distribution"),
        x = i,
        y = "Frequency"
      )
  ) 
}

dev.off()

# Save the filtered object
saveRDS(
  seurat_filt,
  here("NB_ITH", "NBAtlas", "data", "alldata", "01_filt_obj.rds")
)
