# Load libraries
library(tidyverse)
library(readr)
library(DESeq2)

if (!exists("healthy_data_pancreas_mb")) {
  healthy_data_pancreas_mb <- read.delim("data/healthy_data_subset.tsv", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("cancer_mb")) {
  cancer_mb <- read.delim("data/cancer_data_subset.tsv", header = TRUE, sep = "\t", row.names = 1)
}

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

# Run the analysis
dds <- DESeq(dds)
res <- results(dds)

# Check the distribution of padj to define a threshold  
hist(res$padj, 
     breaks = 50, 
     col = "skyblue", 
     main = "padj distribution", 
     xlab = "padj",
     ylab = "Frequency")

res_df <- as_data_frame(res) |> 
  mutate(Significant = padj < 0.001 & abs(res$log2FoldChange) > 2)

# Volcano Plot
ggplot (data = res_df,
        mapping = aes(x = log2FoldChange, 
                      y = -log10(padj),
                      color = Significant)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.001),
             linetype = 'dashed') +
  geom_vline(xintercept = c(log2(0.25), log2(4)),
             linetype = 'dashed') +
  scale_color_manual(values = c("pink", "skyblue")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot") + 
  theme_minimal()

# Filter the results by |logFC| and padj 
filtered_res_signif <- as.data.frame(res[which(res$padj < 0.001 & abs(res$log2FoldChange) > 2), ])

# Save the results in a csv file
write.csv(filtered_res_signif, "data/filtered_deseq2_results.csv", row.names = TRUE)

