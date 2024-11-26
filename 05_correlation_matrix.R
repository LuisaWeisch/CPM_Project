# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(reshape2)


# Load data

# Check if the variables already exists, otherwise read the file

if (!exists("cancer_tpm")) {
  cancer_tpm <- read.delim("_raw/TCGA-PAAD.star_tpm.tsv", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("membrane_data")) {
  membrane_data <- read.delim("_raw/Membrane_lars_ENSG.txt")
}

if (!exists("filtered_data")) {
  filtered_data<-read_csv('data/filtered_deseq2_results.csv')
}


# Format data
rownames(cancer_tpm) <- gsub("\\.\\d+$", "", rownames(cancer_tpm))

# Match TCGA & Deeploc data
data <- cancer_tpm[rownames(cancer_tpm) %in% membrane_data$To,]

# Match TCGA & DESeq2 data
gene_list <- filtered_data$...1
filtered_counts <- data[row.names(data) %in% gene_list, ]


# Transform to fit correllation:
filtered_counts <- t(filtered_counts)

# Getting correlation matrix:
cor_filtered_counts <- cor(filtered_counts, method = "pearson")


# Getting the data in the right format:
melt_data<-melt(cor_filtered_counts)


# Creating a ggheatmap
ggheatmap <- ggplot(melt_data, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)

#remove correlation with self and NA:
clean_matrix <- melt_data[complete.cases(melt_data), ]

#Getting the top correlated gene_pairs:
top_data <- clean_matrix[clean_matrix[, 3] > 0.90, ]
top_corr_data <- top_data[top_data[, 3] != 1, ]
top_corr_data <- top_corr_data[!duplicated(top_corr_data[, 3]), ]

#Saving as a .csv
write.csv(top_corr_data, "data/gene_matches.csv", row.names = FALSE)

