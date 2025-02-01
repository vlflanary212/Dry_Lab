# Author: Victoria Flanary
# Date: 250131
# Objective: Create a gene positions file to run infercnv

# Load in the gtf file
gtf_file <- "~/Library/CloudStorage/Box-Box/Sen_Lab/Computational/Genome/hg38/gencode.v22.annotation.gtf"

gtf <- read.table(
  file = gtf_file,
  sep = "\t",
  header = FALSE,
  comment.char = "#"
)

# Format the gtf file to only include gene and chromosome positions
genes <- gtf %>%
  filter(V3 == "gene") %>%
  mutate(gene_name = sapply(strsplit(V9, ";"), function(x) {
    # Find the part that contains the 'gene_name' field
    gene_field <- grep("gene_name", x, value = TRUE)
    # Extract the actual gene name, removing extra spaces and quotes
    if (length(gene_field) > 0) {
      gsub('gene_name |"| ', '', gene_field)
    } else {
      NA
    }
  })) %>%
  dplyr::select(gene_name, chromosome = V1, start = V4, end = V5)

# Remove duplicate gene entries
genes_unique <- genes[!duplicated(genes$gene_name), ]

# Save the positions object
write.table(
  genes_unique,
  file = here("NB_ITH", "NBAtlas", "04_infercnv", "inputs", "gene_positions.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# End of script
