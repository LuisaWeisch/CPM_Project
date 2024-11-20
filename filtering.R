tcga_data <- read.delim("_raw/TCGA-PAAD.star_tpm.tsv", header = TRUE, sep = "\t", row.names = 1)
deeploc <- read.delim("_raw/Membrane_lars_ENSG.txt")

rownames(tcga_data) <- gsub("\\.\\d+$", "", rownames(tcga_data))
subset_data <- tcga_data[rownames(tcga_data) %in% deeploc$To,]
healthy_data <- read.delim("data/healthy_data_subset.tsv")
