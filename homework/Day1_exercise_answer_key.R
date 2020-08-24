#### QC

## Extract the new metadata from the filtered Seurat object to go through the same plots as with the unfiltered data

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

# 1. Cell counts: do ctrl and stim group have similar cell counts after filtering?
# Yes

# 2. UMI counts (transcripts) per cell: do you observe the removal of cells with less than 500 UMI?
# Yes

# 3. Genes detected per cell: do you observe the removal of cells with less than 250 genes?
# Yes

# 4. UMIs vs. genes detected: do you observe cells with a high number of UMIs but only a few number of genes?
# No

# 5. Mitochondrial counts ratio: do you observe the removal of cells with more than 0.2 mitochondrial counts ratio?
# Yes

# 6. Complexity: do you observe the removal of cells with less than 0.8 log10GenesPerUMI?
# Yes

#### scTransform

# 1. Are the same assays available for the "stim" samples within the split_seurat object? What is the code you used to check that?
# Yes they are
split_seurat$stim@assays

# 2. Any observations for the genes or features listed under "First 10 features:" and the "Top 10 variable features:" for "ctrl" versus "stim"?
# "First 10 features:" - It appears that the same genes are present in both "ctrl" and "stim"
# "Top 10 variable features:" - These are different in the the 2 conditions with some overlap between them
