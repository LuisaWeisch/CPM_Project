library(clusterProfiler)
library(org.Hs.eg.db)

# Load the `res` object
# Assuming `res` is already present in the workspace


# Extract log2FoldChange values and create a ranked list
ranked_genes <- res$log2FoldChange

names(ranked_genes) <- rownames(res)  # Use gene IDs as names

# Remove NA values and sort in decreasing order
ranked_genes <- ranked_genes[!is.na(ranked_genes)]

ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Perform GSEA for Gene Ontology
gsea_results <- gseGO(
  geneList = ranked_genes,              # Ranked list of genes
  OrgDb = org.Hs.eg.db,                 # Annotation database
  ont = "ALL",                          # ALL categories
  keyType = "ENSEMBL",                  # Gene ID type
  minGSSize = 10,                       # Minimum gene set size
  maxGSSize = 500,                      # Maximum gene set size
  pvalueCutoff = 0.05,                  # p-value cutoff
  verbose = TRUE
)
    
# Save results
write.csv(as.data.frame(gsea_results), "data/gsea_results_from_res.csv", row.names = FALSE)

# Print results
gsea_results@result[["Description"]]
gsea_results@result[["setSize"]]


#### GSEA lab 03
library(fgsea)
library(msigdbr)

rankings <- res$stat
names(rankings) <- rownames(res)

rankings <- rankings[!is.na(rankings)]


gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
#gene_sets <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
gene_sets <- split(x = gene_sets$ensembl_gene, f = gene_sets$gs_name)

fgsea_res <- fgsea(pathways = gene_sets, 
                   stats = rankings)

head(fgsea_res)

topPathways <- fgsea_res  |> 
  dplyr::arrange(padj) |> 
  dplyr::slice(1:10)


ggplot(topPathways, 
       mapping = aes(x = reorder(pathway, NES), 
                     y = NES)) +
  geom_col() +
  coord_flip() +
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score (NES)", 
       title = "Top 10 Enriched Pathways") +
  theme_minimal()
