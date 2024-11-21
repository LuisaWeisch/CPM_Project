# Load libraries
library(tidyverse)
library(readr)
library(DESeq2)

# Load TCGA data
tcga_data <- read.delim("_raw/TCGA-PAAD.star_counts.tsv", header = TRUE, sep = "\t", row.names = 1)

# Format data
rownames(tcga_data) <- gsub("\\.\\d+$", "", rownames(tcga_data))

# Load Deeploc data
deeploc <- read.delim("_raw/Membrane_lars_ENSG.txt")

# Match TCGA & Deeploc data
cancer_data_subset <- tcga_data[rownames(tcga_data) %in% deeploc$To,]

#Undo the log2(counts + 1) transformation
cancer_data_subset[,] <- 2^cancer_data_subset[,] - 1

# Load healthy data
healthy_data <- read.delim("_raw/gtex_gene_expected_count", header = TRUE, sep = "\t", row.names = 1)

# Load phenotype data
phenotype_data <- read.delim("_raw/GTEX_phenotype", header = TRUE, sep = "\t", row.names = 1)

# Filter by Pancreas
phenotype_data <- filter(phenotype_data, X_primary_site == "Pancreas") 

# Get sample IDs from the filtered phenotype data
pancreas_samples <- rownames(phenotype_data)

# Filter the healthy data to include only the samples in 'pancreas_samples'
healthy_data_subset <- healthy_data[, colnames(healthy_data) %in% pancreas_samples]

# Format healthy data
rownames(healthy_data_subset) <- gsub("\\.\\d+$", "", rownames(healthy_data_subset))

# Match Healthy & Deeploc data
### ERROR!!!!
healthy_data_subset <- healthy_data[rownames(healthy_data_subset) %in% deeploc$To,]

# Undo the log(counts + 1) transformation
healthy_data_subset[,] <- 2^healthy_data_subset - 1

# Save subset_healthy in data
# write.table(healthy_data_subset, "data/healthy_data_subset.tsv", sep = "\t", quote = FALSE, col.names = NA)

## The datasets have different number of genes apparently, so another subset is made to join them with cbind and 
# perform the DESeq2 analysis

# Identify common genes
common_genes <- intersect(rownames(healthy_data_subset), rownames(cancer_data_subset))

# Subset both datasets to only include common genes
healthy_counts <- healthy_data_subset[common_genes, ]
cancer_counts <- cancer_data_subset[common_genes, ]

# Combine the counts into a single matrix (only the common genes)
counts_matrix <- cbind(healthy_counts, cancer_counts)

counts_matrix <- as.matrix(counts_matrix)  # Ensure counts_matrix is a matrix

# Create a metadata table
metadata <- data.frame(
  Sample = colnames(counts_matrix),                 # Sample names
  Condition = c(rep("Healthy", ncol(healthy_counts)), 
                rep("Cancer", ncol(cancer_counts))) # Group labels
)



# Create a DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ Condition # Use "Condition" as the design variable
)

# Pre-filter low-count genes (optional but recommended)
dds <- dds[rowSums(counts(dds)) > 10, ]
