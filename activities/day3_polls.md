1. During quality control of clustered cells, which of the following do we NOT explore:
  i. whether integration of our data resulted in cell type-specific clusters
  i. whether cells cluster by number of genes, cell cycle or mitochondrial gene expression
  i. whether we observe particular PCs driving specific clusters
  i. **whether we need to remove specific cells based on low UMIs and low gene numbers**
  i. whether we see the cell types that we expect to be present using known cell type markers
  i. whether we see clusters representing doublets
  i. whether our resolution is appropriate for identifying cell-type specific clusters


1. Suppose you are analyzing a single-cell RNAseq data with two experimental conditions: wildtype and knockout. You have finished clustering, and you want to find the markers for a specific cluster. What is the **best** function to use in Seurat?
    1. FindAllMarkers 
    1. **FindConservedMarkers**
    1. FindMarkers
    
2. Which of the following statement about the argument of marker identification is **NOT** correct?
    1. **A logfc.threshold value of 1 indicates that the expression is the same between compared clusters.**
    1. A higher value of min.diff.pct will result in more stringent filtering criteria.
    1. We can use a combination of arguments to adjust filtering stringency. 

3. Which of the following statement about the output of marker identification is **NOT** correct?
    1. We can obtain genes with positive and negative expression changes.
    1. The output p-value for each condition should not be used to interpret significance.
    1. **After identifying marker genes, we normally visualize their expressions in PCA plot.**
