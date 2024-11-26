# Set the output folder
output_folder <- "_raw"

# Create _raw folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Healthy data set - Counts:

url_healthy <- "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_gene_expected_count.gz"
system2("wget", args = c("-P", output_folder, url_healthy))

file_name <- basename(url_healthy)
file_path <- file.path(output_folder, file_name)
system2("gzip", args = c("-d", file_path))

# Cancer data set - Counts:

url_cancer <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.star_counts.tsv.gz"
system2("wget", args = c("-P", output_folder, url_cancer))

file_name <- basename(url_cancer)
file_path <- file.path(output_folder, file_name)
system2("gzip", args = c("-d", file_path))

# Phenotype data set:

url_pheno <- "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz"
system2("wget", args = c("-P", output_folder, url_pheno))

file_name <- basename(url_pheno)
file_path <- file.path(output_folder, file_name)
system2("gzip", args = c("-d", file_path))

### IMPORTANT !!
# Upload the "cell_membrane_proteins_enst" and "Membrane_lars_ENSG" in the _raw folder manually