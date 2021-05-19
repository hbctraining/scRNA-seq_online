### Introduction to scRNA-seq

1. There are many challenges associated with scRNA-seq analysis. Which of these challenges do UMIs attempt to resolve?
    - The large amount of data associated with scRNA-seq
    - Low depth of sequencing resulting in zero counts for many genes
    - **Over-represented transcripts in the library due to PCR artifacts**
    - Temporal changes causing biological variability
    - Batch effects causing technical variability
 
### Raw data to count matrix

1. Which of the following statements is **NOT TRUE** about the 3' end sequencing protocols for single cell RNA-seq?
    - Cheaper per cell cost
    - Larger number of cells sequenced allow for better identity of cell populations
    - **Ideal for detection of isoform-level differences in expression**
    - More accurate quantification due to UMI usage

1. **True**/False. Two reads with the same sample index, cellular barcode and UMI should be counted as a single read.

1. True/**False**. Using single cell RNA-seq we are assaying so many cells per sample, so there is no need to have biological replicates.

1. Which of the following steps are **NOT** involved in creating a cell x gene count matrix?
    - Filtering out unknown cellular barcodes
    - Demultiplexing the samples
    - Mapping/pseudo-mapping to transcriptome
    - **Compute quality metrics for each individual sample**
    - Collapsing UMIs and quantification of reads

