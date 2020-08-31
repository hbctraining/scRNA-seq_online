
1. Which of these statements about the quality control metrics is FALSE?
    1. **The cell count number (also the number of unique cellular barcodes) in a sample should be equal to the number of cells you loaded.**
    1. The proportional histogram of genes detected per cell can often be bimodal.
    1. When evaluating the relationship between the nUMIs and the number of genes detected we examine the slope of the line, and any scatter of data points in the bottom right hand quadrant of the plot.
    1. Mitochondrial counts ratio can identify whether there is a large amount of mitochondrial contamination from dead or dying cells.
    1. Cells with a low novelty score are those that have a less complex RNA species than other cells. 

1. True or **False**. If a sample fails to meet the threshold for any one quality metric, it is best to discard it before moving forward.

1. Which of these statements about the processing steps post-QC is FALSE? 
    1. sctransform is the normalization method from Seurat implemented in our worflow which simultaneously performs variance stabilization and regresses out unwanted variation.
    1. **Cell cycle is the most common effect that we regress out from our data, regardless of what we observe when exploring the data.**
    1. Plotting a PCA helps to identify sources of unwanted variation.
    1. A source of unwanted variation in one dataset could be considered a biological effect of interest in another dataset, and so it’s important to have a good idea of expectations for your data.

1. **True** or False. When calculating the principal components, the PC that deals with the largest variation in the dataset is designated PC1.

1. There are as many PCs as there are ______ in a given scRNA-seq analysis.
    1. **Cells**
    2. Genes
