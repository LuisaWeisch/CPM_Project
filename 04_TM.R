library(biomaRt)

### THIS HAS BEEN DONE LOCALLY BECAUSE IT IS NOT WORKING CORRECTLY ON THE SERVER

# Connect to Ensembl server
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# DESeq2 results
degs_df <- read.delim("data/filtered_deseq2_results.csv", header = TRUE, sep = "\t", row.names = 1)
# Ensembl gene IDs 
ensembl_gene_ids <- rownames(degs_df)

# Query Ensembl to get the protein sequences for the given ENSEMBL gene IDs
protein_sequences <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "peptide"),
  filters = "ensembl_gene_id",
  values = ensembl_gene_ids,
  mart = ensembl
)

write.csv(protein_sequences, "data/protein_sequences.csv", row.names = FALSE)

# Protein sequences in FASTA format
fasta_lines <- c()
for (i in 1:nrow(protein_sequences)) {
  fasta_lines <- c(fasta_lines, paste0(">", protein_sequences$peptide[i]))
}

write(fasta_lines, file = "protein_sequences_fasta.txt")

