# Author: Victoria Flanary
# Date: 250202
# Objective: Annotate peaks for untreated neuroblastoma cell line ATAC-seq data

# Load packages
library(tidyverse)
library(biovizBase)
library(GenomicFeatures)
library(genomation)
library(ChIPseeker)
library(rtracklayer)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(here)
library(ChrAccR)

# Set working directory
wd <- here("13cRA", "atac-seq", "untreated")

# Create a custom function for reading in narrow peaks as GRanges
NarrowToGRanges <- function(data_dir, study, sample_name) {
  peak_dir <- file.path(data_dir, study, "peaks")
  peak_file <- file.path(peak_dir, paste0(sample_name, "_peaks.narrowPeak"))
  narrow_granges <- readNarrowPeak(peak_file)
  mcols(narrow_granges) <- cbind(mcols(narrow_granges), DataFrame(sample = sample_name))
  return(narrow_granges)
}

data_dir <- "~/Library/CloudStorage/Box-Box/Data/atac-seq/untreated"

# hNCC data from GSE108517
## Load the sample list
hncc_samples <- readLines(here(wd, "GSE108517", "cell_lines.txt"))

# Create GRanges object for this study
hncc_peak_list <- NarrowToGRanges(data_dir, study = "GSE108517", sample_name = hncc_samples)
names(hncc_peak_list) <- hncc_samples

# Neuroblastoma data
## Load sample list
nb_samples <- readLines(here(wd, "GSE138293", "cell_lines.txt"))

# Create GRanges object for this study
nb_peak_list <- lapply(nb_samples, function(sample) {
  NarrowToGRanges(data_dir, study = "GSE138293", sample_name = sample)
})
names(nb_peak_list) <- nb_samples

# Supplementary cell line ATAC-seq data generated internally
## Load sample list
sen_samples <- readLines(here(wd, "Sen_240503", "samples.txt"))

## Peak calling failed for GIMEN - rm for now
sen_samples <- setdiff(sen_samples, "GIMEN")

# Create GRanges object for this study
sen_peak_list <- lapply(sen_samples, function(sample) {
  NarrowToGRanges(data_dir, study = "Sen_240503", sample_name = sample)
})
names(sen_peak_list) <- sen_samples

# Merge peak lists
peak_list <- c(sen_peak_list, nb_peak_list, hncc_peak_list)

# Fix names for the peak list
names <- setdiff(names(peak_list), "")
names <- c(names, "hNCC")
names(peak_list) <- names

# Append lists into a single GRanges object
peaks <- peak_list |> GRangesList() |> unlist()

# Format peak names to chr_start_end
peaks$name <- paste0(seqnames(peaks), "_", start(peaks), "_", end(peaks))

# Filter to only include canonical autosomes
chromosomes <- paste0("chr", 1:22)
peaks_filt <- peaks[seqnames(peaks) %in% chromosomes]
seqlevels(peaks_filt) <- intersect(seqlevels(peaks_filt), chromosomes)

# Annotate peaks
## Load annotations
txdb_file <- "~/Library/CloudStorage/Box-Box/Sen_Lab/Computational/Genome/hg38/gencode.v22.annotation.txdb"
txdb <- loadDb(txdb_file)

## Get genomic tiling
### Sort the peaks
peaks_sorted <- sort(peaks_filt)

### Get list of 1kb non-overlapping genomic regions
tiling <- muRtools::getTilingRegions("hg38", width = 1000L, onlyMainChrs = TRUE)

### Filter tiles for those that overlap the ATAC-seq peaks
tiling <- subsetByOverlaps(tiling, peaks_sorted)

### Add fragment names to tiling
tiling$name <- paste0(seqnames(tiling), "_", start(tiling), "_", end(tiling))

## Annotate peaks
peaks_anno <- annotatePeak(
  tiling,
  tssRegion = c(-1000, 500),
  TxDb = txdb, 
  annoDb = "org.Hs.eg.db"
) |>
  as.data.frame()

## Check all regions
unique(peaks_anno$annotation)

## Save the annotated peaks
saveRDS(
  peaks_anno,
  here(wd, "data", "annotated_peaks.rds")
)

# Subset regions into separate GRanges
promoters <- peaks_anno[grep("Promoter", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

exons <- peaks_anno[grep("Exon", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

introns <- peaks_anno[grep("Intron", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

distal <- peaks_anno[grep("Distal Intergenic", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

downstream <- peaks_anno[grep("Downstream (<=300bp)", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

utr_3 <- peaks_anno[grep("3' UTR", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

utr_5 <- peaks_anno[grep("5' UTR", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Create the DsATAC object
## Create a vector for all bam files

## Generate a combined metadata dataframe
metadata <- data.frame(
  sample = names(peak_list)
  age = ,
  sex = ,
  mycn_status = ,
  phenotype = ,
  bam_dir =
)

