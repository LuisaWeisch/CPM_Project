# Load libraries
library(dplyr)
library(readr)

# Load TCGA data
tcga_data <- read.delim("_raw/TCGA-PAAD.star_tpm.tsv", header = TRUE, sep = "\t", row.names = 1)

# Format data
rownames(tcga_data) <- gsub("\\.\\d+$", "", rownames(tcga_data))

# Load Deeploc data
deeploc <- read.delim("_raw/Membrane_lars_ENSG.txt")

# Match TCGA & Deeploc data
subset_data <- tcga_data[rownames(tcga_data) %in% deeploc$To,]

# Load healthy data
healthy_data <- read.delim("_raw/gtex_RSEM_gene_tpm", header = TRUE, sep = "\t", row.names = 1)

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
healthy_data_subset <- tcga_data[rownames(healthy_data_subset) %in% deeploc$To,]

# Save subset_healthy in data
write.table(healthy_data_subset, "data/healthy_data_subset.tsv", sep = "\t", quote = FALSE, col.names = NA)



