# Author: Victoria Flanary
# Date: 250128
# Objective: Run ATAC-seq QC metrics (TSS enrichment score and fragment size
# distribution) using ChrAccR

# Load packages
library(tidyverse)
library(ggpubr)
library(here)
library(ChrAccR)

# ChrAccR config settings
theme_set(muRtools::theme_nogrid())