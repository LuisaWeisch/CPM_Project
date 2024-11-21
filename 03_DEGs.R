# Load libraries
library(tidyverse)
library(readr)
library(DESeq2)

# Identify common genes
common_genes <- intersect(rownames(healthy_data_pancreas_mb), rownames(cancer_mb))

# Subset both datasets to only include common genes
healthy_counts <- healthy_data_pancreas_mb[common_genes, ]
cancer_counts <- cancer_mb[common_genes, ]

# Combine the counts into a single matrix (only the common genes)
counts_matrix <- cbind(healthy_counts, cancer_counts)

# Counts as integers because no decimal points are accepted
counts_matrix <- round(counts_matrix)

# Create a metadata table
metadata <- data.frame(
  Sample = colnames(counts_matrix),                 # Sample names
  Condition = c(rep("Healthy", ncol(healthy_counts)), 
                rep("Cancer", ncol(cancer_counts))) # Group labels
)

# Converting the condition to a factor
metadata$Condition <- factor(metadata$Condition, levels = c("Healthy", "Cancer"))


# Create a DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ Condition # Use "Condition" as the design variable
)

# Pre-filter low-count genes (optional but recommended)
dds <- dds[rowSums(counts(dds)) > 10, ]

#Run the analysis
dds <- DESeq(dds)
res <- results(dds)

#Filter the results by |logFC| > 1 and padj < 0.01
filtered_res <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# Save the results in a csv file
write.csv(as.data.frame(filtered_res), "filtered_deseq2_results.csv", row.names = TRUE)

