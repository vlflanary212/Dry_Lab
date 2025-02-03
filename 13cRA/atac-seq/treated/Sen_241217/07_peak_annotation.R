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

# Set filepaths
wd <- here("13cRA", "atac-seq", "treated", "Sen_241217")
data_dir <- "/data/scratch/flanary/atac-seq/treated/Sen_241217"

# Create a custom function for reading in narrow peaks as GRanges
NarrowToGRanges <- function(data_dir, sample_name) {
  peak_dir <- file.path(data_dir, "peaks")
  peak_file <- file.path(peak_dir, paste0(sample_name, "_peaks.narrowPeak"))
  narrow_granges <- readNarrowPeak(peak_file)
  mcols(narrow_granges) <- cbind(mcols(narrow_granges), DataFrame(sample = sample_name))
  return(narrow_granges)
}

# Internal data: 4 cell lines treated with 5 microM 13cRA for 24 hours
## Load the sample list
samples <- readLines(here(wd, "cell_lines.txt"))

# Create GRanges object for this study
peak_list <- lapply(samples, function(sample) {
  NarrowToGRanges(data_dir, sample_name = sample)
})
names(peak_list) <- samples

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
txdb_file <- "/data/project/sen-lab/genome/hg38/gencode.v22.annotation.txdb"
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

utr_3 <- peaks_anno[grep("3' UTR", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

utr_5 <- peaks_anno[grep("5' UTR", peaks_anno$annotation), ] |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Create the DsATAC object
## Create a vector for all bam files
bam_dir <- file.path(data_dir, "final_bam")
bam_files <- file.path(bam_dir, paste0(samples, "_final.bam"))
  
## Generate a combined metadata dataframe
sample_anno <- data.frame(
  sample = c("IMR-5", "SHEP", "SK-N-AS", "SH-SY5Y"),
  age_months = c("13", "48", "72", "48"),
  sex = c("Male", "Female", "Female", "Female"),
  mycn_status = c("Amplified", "Non-Amplified", "Non-Amplified", "Non-Amplified"),
  phenotype = c("ADR", "MES", "MES", "ADR"),
  bam_files = bam_files
)

## Create the DsATAC object
dsa <- DsATAC.bam(sample_anno, "bam_files", "hg38", 
                  regionSets = list(
                    promoters = promoters,
                    exons = exons,
                    introns = introns,
                    distal = distal,
                    utr_3 = utr_3,
                    utr_5 = utr_5
                  ), 
                  sampleIdCol = "sample", diskDump = FALSE)

## Save the DsATAC object
dest <- file.path(wd, "DsAtacDataset")
saveDsAcc(dsa, dest)

# End of script
