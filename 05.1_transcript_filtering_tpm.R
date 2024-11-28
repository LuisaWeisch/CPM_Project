# Load libraries
library(dplyr)
library(readr)

# Load data
# Check if the variables already exists, otherwise read the file
if (!exists("cancer_data")) {
  cancer_data <- read.delim("_raw/tcga_rsem_isoform_tpm", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("membrane_data")) {
  membrane_data <- read.delim("_raw/Membrane_lars_ENSG.txt")
}

if (!exists("healthy_data")) {
  healthy_data <- read.delim("_raw/gtex_rsem_isoform_tpm", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("phenotype_data")) {
  phenotype_data <- read.delim("_raw/GTEX_phenotype", header = TRUE, sep = "\t", row.names = 1)
}

list_ids <- read.csv("_raw/list_ids.csv")
###################
#   CANCER DATA   #     
###################

# Format data
rownames(cancer_pancreas) <- gsub("\\.\\d+$", "", rownames(cancer_pancreas))

cancer_pancreas <- cancer_data[,colnames(cancer_data) %in% list_ids$Sample_ID]
write.csv(cancer_pancreas, "data/cancer_pancreas_transcripts.csv")

# Match TCGA & Deeploc data
cancer_mb <- cancer_pancreas[rownames(cancer_pancreas) %in% membrane_data$From,]

# Save subset_healthy in data
output_folder <- "data"

# Create data folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

write.table(cancer_mb, "data/cancer_mb_transcripts.tsv", sep = "\t", quote = FALSE, col.names = NA)

####################
#   HEALTHY DATA   #     
####################

# Filter by Pancreas
pancreas_data <- phenotype_data |> 
  filter(x_primary_site == "Pancreas")

# pancreas_data <- filter(phenotype_data, X_primary_site == "Pancreas") 
                        
# Get sample IDs from the filtered phenotype data
pancreas_samples <- rownames(pancreas_data)

# Filter the healthy data to include only the samples in 'pancreas_samples'
healthy_data_pancreas <- healthy_data[, colnames(healthy_data) %in% pancreas_samples]

# Match Healthy & Deeploc data
healthy_data_pancreas_mb <- healthy_data_pancreas[rownames(healthy_data_pancreas) %in% membrane_data$From,]

# Save subset_healthy in data
output_folder <- "data"

# Create data folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

write.table(healthy_data_pancreas_mb, "data/healthy_transcripts_tpm.tsv", sep = "\t", quote = FALSE, col.names = NA)



