# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(reshape2)


# Load data

# Check if the variables already exists, otherwise read the file

if (!exists("cancer_count")) {
  cancer_count <- read.delim("_raw/TCGA-PAAD.star_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("membrane_data")) {
  membrane_data <- read.delim("_raw/Membrane_lars_ENSG.txt")
}

if (!exists("filterede_data")) {
  filtered_data<-read_csv('filtered_deseq2_results.csv')
}




# Format data
rownames(cancer_count) <- gsub("\\.\\d+$", "", rownames(cancer_count))

# Match TCGA & Deeploc data
data <- cancer_count[rownames(cancer_count) %in% membrane_data$To,]

# Match TCGA & DESeq2 data

gene_list <- filtered_data$...1
filtered_counts <- data[row.names(data) %in% gene_list, ]


#transform to fit correllation:
filtered_counts <- t(filtered_counts)

# Getting correlation matrix:
cor_filtered_counts <- cor(filtered_counts, method = "spearman")


#getting the data in the right format:
melt_data<-melt(cor_filtered_counts)


# Create a ggheatmap
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






