### Introduction to scRNA-seq

1. Which of the following does **NOT** contribute to the complexity of scRNA-seq analysis:
    - i. scRNA-seq generates a lot of data
    - i. Low depth of sequencing showing zero counts for many genes
    - i. **Transcripts over-represented in the library due to PCR artifacts**
    - i. Temporal changes causing biological variability
    - i. Batch effects causing technical variability


1. **Choose the most appropriate approach.** You do not expect a lot of cell type heterogeneity in your samples and you would like to determine the effects of a treatment on gene expression.
    - i. **Bulk RNA-seq**
    - i. scRNA-seq
    
1. **Choose the most appropriate approach.** You would like to characterize the cell types present in a tissue.
    - i. Bulk RNA-seq
    - i. **scRNA-seq**
    
1. **Choose the most appropriate approach.** You would like to determine the genes turned on during differentiation of a particular cell type.
    - i. Bulk RNA-seq
    - i. **scRNA-seq**
    
1. **Choose the most appropriate approach.** You would like to identify a strong biomarker for a disease in which we expect a lot of heterogeneity across cell types.
    - i. Bulk RNA-seq
    - i. **scRNA-seq**
    
1. **Choose the most appropriate approach.** You FACS sort your samples for a cell type of interest and would like to determine the differential expression of gene isoforms.
    - i. **Bulk RNA-seq**
    - i. scRNA-seq
    
1. **Choose the most appropriate approach.** You would like to determine differentially expressed genes between 2 conditions in heterogenous tissue; but, your budget can only handle either 1 scRNA-seq sample/condition or multiple bulk RNA-seq samples/condition. 
    - i. **Bulk RNA-seq**
    - i. scRNA-seq
 
### Raw data to count matrix

1. Which of the following statements is **NOT TRUE** about the 3' end sequencing protocols for single cell RNA-seq"?
    - Cheaper per cell cost
    - Larger number of cells sequenced allow for better identity of cell populations
    - **Ideal for detection of isoform-level differences in expression**
    - More accurate quantification due to use UMIs

1. **True**/False. Two reads with the same sample index, cellular barcode and UMI should be counted as a single read.

1. True/**False**. Using single cell RNA-seq we are assaying so many cells per sample, so there is no need to have biological replicates.

1. Which of the following steps are **NOT** involved in creating a cell x gene count matrix?
    - Filtering out unknown cellular barcodes
    - Demultiplexing the samples
    - Mapping/pseudo-mapping to transcriptome
    - **Compute quality metrics for each individual sample**
    - Collapsing UMIs and quantification of reads

