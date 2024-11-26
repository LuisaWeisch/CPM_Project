#source("eleni_notes.R") #to run
library("readr")
tcga_data <- read.delim("_raw/TCGA-PAAD.star_tpm.tsv", header = TRUE, sep = "\t", row.names = 1)
deeploc <- read.delim("_raw/Membrane_lars_ENSG.txt")

rownames(tcga_data) <- gsub("\\.\\d+$", "", rownames(tcga_data))
subset_data <- tcga_data[rownames(tcga_data) %in% deeploc$To,]

GTEx_pheno <- read.delim("GTEX_phenotype")

healthy <- fromJSON("gtex_RSEM_gene_tpm.json")
