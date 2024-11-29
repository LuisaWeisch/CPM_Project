# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(reshape2)


# Load data

# Check if the variables already exists, otherwise read the file

if (!exists("cancer_trans")) {
  cancer_trans <- read.delim("_raw/tcga_rsem_isoform_tpm", header = TRUE, sep = "\t", row.names = 1)
}

if (!exists("filtered_trans")) {
  filtered_trans<-read_csv('data/filtered_deseq2_transcripts_results.csv')
}


rownames(cancer_trans) <- gsub("\\.\\d+$", "", rownames(cancer_trans))

# Match TCGA & DESeq2 data
gene_list <- filtered_trans$...1
filtered_data_trans <- cancer_trans[row.names(cancer_trans) %in% gene_list, ]


# Transform to fit correllation:
filtered_data <- t(filtered_data_trans)

# Getting correlation matrix:
filtered_data <- cor(filtered_data, method = "pearson")

# Getting the data in the right format:
melt_data<-melt(filtered_data)


#remove correlation with self and NA:
clean_matrix <- melt_data[complete.cases(melt_data), ]


#Getting the top 0.9 genes:
top_data <- clean_matrix[clean_matrix[, 3] > 0.80, ]
top_data <- top_data[top_data[, 3] != 1, ]
gene_list2 <- unique(top_data[, 1])


#Filter DATA by second gene list with the highest correlation genes:

filtered_data <- filtered_data_trans[row.names(filtered_data_trans) %in% gene_list2, ]

# Transform to fit correllation:
filtered_data <- t(filtered_data)

# Getting correlation matrix:
filtered_data <- cor(filtered_data, method = "pearson")

# Getting the data in the right format:
melt_data<-melt(filtered_data)


ggheatmap <- ggplot(melt_data, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)

# get list of all genes (unique):

write.csv(gene_list2, "data/gene_list_trans.csv", row.names = FALSE)

#Get only unique pairs:
top_data_90 <- clean_matrix[clean_matrix[, 3] > 0.90, ]
top_data_90 <- top_data_90[top_data_90[, 3] != 1, ]
top_corr_data <- top_data[!duplicated(top_data[, 3]), ]
top_corr_data_90 <- top_data[!duplicated(top_data_90[, 3]), ]
gene_list2 <- unique(top_data_90[, 1])



#Saving as a .csv
write.csv(top_corr_data, "data/gene_matches_trans.csv", row.names = FALSE)
write.csv(top_corr_data_90, "data/gene_matches_trans_90.csv", row.names = FALSE)

write.csv(gene_list2, "data/gene_list_trans_90.csv", row.names = FALSE)


###########################

#          HEALTY         #

###########################
cancer_trans <- NULL

if (!exists("healty_trans")) {
  healty_trans <- read.delim("_raw/gtex_rsem_isoform_tpm", header = TRUE, sep = "\t", row.names = 1)
}

rownames(healty_trans) <- gsub("\\.\\d+$", "", rownames(healty_trans))

write.csv(gene_list2, "data/gene_list_trans.csv", row.names = FALSE)


gene_list2 <- read.csv("data/gene_list_trans.csv")

gene_list2 <- unique(gene_list2[, 1])

filtered_healty_trans <- healty_trans[row.names(healty_trans) %in% gene_list2, ]



# Transform to fit correllation:
filtered_healty <- t(filtered_healty_trans)

# Getting correlation matrix:
filtered_healty <- cor(filtered_healty, method = "pearson")

# Getting the data in the right format:
melt_data_h<-melt(filtered_healty)


#remove correlation with self and NA:
clean_matrix_h <- melt_data_h[complete.cases(melt_data_h), ]


#Getting the top 0.9 genes:
top_data_h <- clean_matrix_h[clean_matrix_h[, 3] > 0.80, ]
top_data_h <- top_data[top_data_h[, 3] != 1, ]
gene_list_h <- unique(top_data_h[, 1])


#Filter DATA by second gene list with the highest correlation genes:

filtered_healty <- filtered_healty_trans[row.names(filtered_healty_trans) %in% gene_list_h, ]

# Transform to fit correllation:
filtered_healty <- t(filtered_healty)

# Getting correlation matrix:
filtered_healty <- cor(filtered_healty, method = "pearson")

# Getting the data in the right format:
melt_data_h<-melt(filtered_healty)


ggheatmap <- ggplot(clean_matrix_h, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)

# get list of all genes (unique):

write.csv(gene_list_h, "data/gene_list_trans_h.csv", row.names = FALSE)

#Get only unique pairs:
top_data_90_h <- clean_matrix_h[clean_matrix_h[, 3] > 0.90, ]
top_data_90_h <- top_data_90_h[top_data_90_h[, 3] != 1, ]
top_corr_data_h <- top_data_h[!duplicated(top_data_h[, 3]), ]
top_corr_data_90_h <- top_data_90_h[!duplicated(top_data_90_h[, 3]), ]
gene_list2_h <- unique(top_data_90_h[, 1])



#Saving as a .csv
write.csv(top_corr_data_h, "data/gene_matches_trans_h.csv", row.names = FALSE)
write.csv(top_corr_data_90_h, "data/gene_matches_trans_90_h.csv", row.names = FALSE)
write.csv(clean_matrix_h, "data/gene_matches_trans_h_all.csv", row.names = FALSE)

write.csv(gene_list2, "data/gene_list_trans_90_h.csv", row.names = FALSE)







