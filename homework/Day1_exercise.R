Extract the new metadata from the filtered Seurat object to go through the same plots as with the unfiltered data

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

Perform all of the same plots as with the unfiltered data and answer the following questions:
# 1. Cell counts: do ctrl and stim group have similar cell counts after filtering?

# 2. UMI counts (transcripts) per cell: do you observe the removal of cells with less than 500 UMI?

# 3. Genes detected per cell: do you observe the removal of cells with less than 250 genes?

# 4. UMIs vs. genes detected: do you observe cells with a high number of UMIs but only a few number of genes?

# 5. Mitochondrial counts ratio: do you observe the removal of cells with more than 0.2 mitochondrial counts ratio?

# 6. Complexity: do you observe the removal of cells with less than 0.8 log10GenesPerUMI?
