# Author: Victoria Flanary
# Date: 250130
# Objective: Plot neuroendocrine marker genes

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(tidyverse)
library(here)
library(future)

# Future parameters
options(future.globals.maxSize = 256 * 1024^3)  # 256 GB maximum RAM
plan("multisession", workers = 4)  # parallelize across 4 cores

# Set filepaths
data_dir <- here("NB_ITH", "NBAtlas", "data", "alldata")
results_dir <- here("NB_ITH", "NBAtlas", "04_cluster_annotation")

# Load data
seurat_obj <- readRDS(here(data_dir, "05_harmony_clust.rds"))

# Plot Genes
pdf(
  here(results_dir, "cluster1_feature_plots.pdf"),
  height = 8, width = 20
)  # differentiating neuronal-like nb?

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c("MYCN", "NEFL", "NYAP2", "LY6H", 
                           "SIX3", "PFN2", "NEFM"))
)

dev.off()


pdf(
  here(results_dir, "cluster13_feature_plots.pdf"),
  height = 12, width = 20
) # progenitor population with neuropeptide function?

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c("FOXR2", "ZFP42", "CYP26B1", "POU3F1",
                           "NPFF", "HCRT", "CDH22", "FZD8", 
                           "MST1L"))
)  # CYP26B1 crucial for retinoic acid metabolism

dev.off()


pdf(
  here(results_dir, "cluster14_feature_plots.pdf"),
  height = 12, width = 20
) # developing neurons defined by synaptic formation, axon guidance, and early differentiation

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c("EIF2S3L", "MIR490", "NBAS", "TLL2", 
                           "CORIN", "PCDH11Y", "ST6GALNAC3", "DSCAM",
                           "GREB1L", "SRPX2"))
) 


pdf(
  here(results_dir, "cluster15_feature_plots.pdf"),
  height = 20, width = 20
) # more developing neurons

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c("TSPAN8", "TP53TG3C", "ADCY2", "LGR5",
                           "PDE11A", "NRG3", "SLC44A5", "KCNH5", 
                           "TMC2", "FREM1", "PRR16", "TCAF2", 
                           "TPTE2"))
) # SLC44A5 plays a role in choline transport

dev.off()


pdf(
  here(results_dir, "cluster16_feature_plots.pdf"),
  height = 12, width = 20)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c("RET", "ERBB4", "OPRM1", "IGF2",
                           "LHX9", "KCNT1", "STRA6", "CACNA2D2",
                           "ADARB2", "DSCAML1"))
) # ERBB4 expression, bridge cells?
# STRA6 also involved in retinoic acid transport

dev.off()

pdf(
  here(results_dir, "cluster16_feature_plots.pdf"),
  height = 12, width = 20
)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c("RET", "ERBB4", "OPRM1", "IGF2",
                           "LHX9", "KCNT1", "STRA6", "CACNA2D2",
                           "ADARB2", "DSCAML1"))
) # ERBB4 expression, bridge cells?
# STRA6 also involved in retinoic acid transport

dev.off()

pdf(
  here(results_dir, "cluster17_feature_plots.pdf"),
  height = 12, width = 20
)  # metabolically active neuronal phenotype with metoxidative stress response, DNA repair, and neuroprotection

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c("LAMP5", "GTSF1", "SLC5A12", "XRCC2",
                           "GSTT2B", "NMNAT3", "TDRD12", "SKA1",
                           "PCDHA10", "CRYBA2"))
) 
# SLC5A12 involved in lactate metabolism
# XRCC2 involved in DNA repair
# GSTT2B functions in oxidative stress resistance
dev.off()


pdf(
  here(results_dir, "cluster18_feature_plots.pdf"),
  height = 12, width = 20
) 

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = c())
) 
dev.off()
