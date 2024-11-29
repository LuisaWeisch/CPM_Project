# Load libraries
library(dplyr)
library(readr)

# Load data
# Check if the variables already exists, otherwise read the file
if (!exists("cancer_data")) {
  cancer_data <- read.delim("_raw/TCGA-PAAD.star_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("membrane_data")) {
  membrane_data <- read.delim("_raw/Membrane_lars_ENSG.txt")
}

if (!exists("healthy_data")) {
  healthy_data <- read.delim("_raw/gtex_gene_expected_count", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("phenotype_data")) {
  phenotype_data <- read.delim("_raw/GTEX_phenotype", header = TRUE, sep = "\t", row.names = 1)
}


###################
#   CANCER DATA   #     
###################

# Format data
rownames(cancer_data) <- gsub("\\.\\d+$", "", rownames(cancer_data))

# Match TCGA & Deeploc data
cancer_mb <- cancer_data[rownames(cancer_data) %in% membrane_data$To,]

# Save subset_healthy in data
output_folder <- "data"

# Create data folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

write.table(cancer_mb, "data/cancer_data_subset.tsv", sep = "\t", quote = FALSE, col.names = NA)

####################
#   HEALTHY DATA   #     
####################

# Filter by Pancreas
pancreas_data <- filter(phenotype_data, X_primary_site == "Pancreas") 
                        
# Get sample IDs from the filtered phenotype data
pancreas_samples <- rownames(pancreas_data)

## !! The sample IDs in healthy_data use dots (.), while the pancreas samples use hyphens (-). 
# Replace the dots in healthy_data column names with hyphens:
colnames(healthy_data) <- gsub("\\.", "-", colnames(healthy_data))

# Filter the healthy data to include only the samples in 'pancreas_samples'
healthy_data_pancreas <- healthy_data[, colnames(healthy_data) %in% pancreas_samples]

# Format healthy data gene IDs
rownames(healthy_data_pancreas) <- gsub("\\.\\d+$", "", rownames(healthy_data_pancreas))

# Match Healthy & Deeploc data
healthy_data_pancreas_mb <- healthy_data_pancreas[rownames(healthy_data_pancreas) %in% membrane_data$To,]

# Save subset_healthy in data
output_folder <- "data"

# Create data folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

write.table(healthy_data_pancreas_mb, "data/healthy_data_subset.tsv", sep = "\t", quote = FALSE, col.names = NA)



