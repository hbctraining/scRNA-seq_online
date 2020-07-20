---
title: "Single-cell RNA-seq: Integration"
author: "Mary Piper, Meeta Mistry, Radhika Khetani"
date: Tuesday, February 25th, 2020
---

Approximate time: 90 minutes

## Learning Objectives:

* Perform integration of cells across conditions using the most variant genes to identify cells most similar to each other


# Single-cell RNA-seq clustering analysis: Integration


<img src="../img/sc_workflow_integration.png" width="800">

***

_**Goals:**_ 

 - _To **align similar cells** across conditions._

_**Challenges:**_
 
 - _**Aligning cells of similar cell types** so that we do not have cell-type specific clustering downstream_

_**Recommendations:**_
 
 - _Go through the analysis without integration first to determine whether integration is necessary_

***

## To integrate or not to integrate?

Generally, we always look at our clustering without integration before deciding whether we need to perform any alignment. If we had performed the normalization on both conditions together in a Seurat object and visualized the similarity between cells, we would have seen condition-specific clustering:

<p align="center">
<img src="../img/unintegrated_umap.png" width="400">
</p>

Condition-specific clustering of the cells indicates that we need to integrate the cells across conditions to ensure that cells of the same cell type cluster together. In this lesson, we will cover the integration of our samples across conditions, which is adapted from the [Seurat v3 Guided Integration Tutorial](https://satijalab.org/seurat/v3.0/immune_alignment.html).

> _**NOTE:** Seurat has a [vignette](https://satijalab.org/seurat/v3.1/sctransform_vignette.html) for how to run through the workflow without integration. The workflow is fairly similar to this workflow, but the samples would not necessarily be split in the beginning and integration would not be performed._
>
> _It can help to first run conditions individually if unsure what clusters to expect or expecting some different cell types between conditions (e.g. tumor and control samples), then run them together to see whether there are condition-specific clusters for cell types present in both conditions. Oftentimes, when clustering cells from multiple conditions there are condition-specific clusters and integration can help ensure the same cell types cluster together._


## **Integrate** or align samples across conditions using shared highly variable genes

_**If cells cluster by sample, condition, batch, dataset, or modality, this integration step can greatly improve the clustering and the downstream analyses**._

To integrate, we will use the shared highly variable genes (identified using SCTransform) from each group, then, we will "integrate" or "harmonize" the groups to overlay cells that are similar or have a "common set of biological features" between groups. For example, we could integrate across:

- Different **conditions** (e.g. control and stimulated)
	<img src="../img/seurat_condition_integ.png" width="800">

- Different **datasets** (e.g. scRNA-seq from datasets generated using different library preparation methods on the same samples)
	<img src="../img/seurat_dataset_integ.png" width="800">

- Different **modalities** (e.g. scRNA-seq and scATAC-seq)
	<img src="../img/seurat_modality_integ.png" width="800">

Integration is a powerful method that uses these shared sources of greatest variation to identify shared subpopulations across conditions or datasets [[Stuart and Bulter et al. (2018)](https://www.biorxiv.org/content/early/2018/11/02/460147)]. The goal of integration is to ensure that the cell types of one condition/dataset align with the same celltypes of the other conditions/datasets (e.g. control macrophages align with stimulated macrophages).

Specifically, this integration method expects "correspondences" or **shared biological states** among at least a subset of single cells across the groups. The steps in the integration analysis are outlined in the figure below:

<p align="center">
<img src="../img/integration.png" width="600">
</p>

_**Image credit:** Stuart T and Butler A, et al. Comprehensive integration of single cell data, bioRxiv 2018 (https://doi.org/10.1101/460147)_

The different steps applied are as follows:

1. Perform **canonical correlation analysis (CCA):**
	
	CCA identifies shared sources of variation between the conditions/groups. It is a form of PCA, in that it **identifies the greatest sources of variation** in the data, but only if it is **shared or conserved** across the conditions/groups (using the 3000 most variant genes from each sample).
	
	This step roughly aligns the cells using the greatest shared sources of variation.

	> _**NOTE:** The shared highly variable genes are used because they are the most likely to represent those genes distinguishing the different cell types present._

2. **Identify anchors** or mutual nearest neighbors (MNNs) across datasets (sometimes incorrect anchors are identified):
	
	MNNs can be thought of as 'best buddies'. For each cell in one condition:
	- The cell's closest neighbor in the other condition is identified based on gene expression values - it's 'best buddy'.
	- The reciprical analysis is performed, and if the two cells are 'best buddies' in both directions, then those cells will be marked as **anchors** to 'anchor' the two datasets together.
	
	> "The difference in expression values between cells in an MNN pair provides an estimate of the batch effect, which is made more precise by averaging across many such pairs. A correction vector is obtained and applied to the expression values to perform batch correction." [[Stuart and Bulter et al. (2018)](https://www.biorxiv.org/content/early/2018/11/02/460147)]. 

3. **Filter anchors** to remove incorrect anchors:
	
	Assess the similarity between anchor pairs by the overlap in their local neighborhoods (incorrect anchors will have low scores) - do the adjacent cells have 'best buddies' that are adjacent to each other?

4. **Integrate** the conditions/datasets:

	Use anchors and corresponding scores to transform the cell expression values, allowing for the integration of the conditions/datasets (different samples, conditions, datasets, modalities)

	> _**NOTE:** Transformation of each cell uses a weighted average of the two cells of each anchor across anchors of the datasets. Weights determined by cell similarity score (distance between cell and k nearest anchors) and anchor scores, so cells in the same neighborhood should have similar correction values._

	**If cell types are present in one dataset, but not the other, then the cells will still appear as a separate sample-specific cluster.**


Now, using our SCTransform object as input, let's perform the integration across conditions.

First, we need to specify that we want to use all of the 3000 most variable genes identified by SCTransform for the integration. By default, this function only selects the top 2000 genes.

```r
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
```

Then, we need to **prepare the SCTransform object** for integration.

```r                                          
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
```

Now, we are going to **perform CCA, find the best buddies or anchors and filter incorrect anchors**. For our dataset, this will take up to 15 minutes to run. *Also, note that the progress bar in your console will stay at 0%, but know that it is actually running.*

```r
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
```

Finally, we can **integrate across conditions**.

```r
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
```

This would often be a good place to **save the R object**.

```r
# Save integrated seurat object
saveRDS(seurat_integrated, "results/integrated_seurat.rds")
```

### UMAP visualization

After integration, to visualize the integrated data we can use dimensionality reduction techniques, such as PCA and Uniform Manifold Approximation and Projection (UMAP). While PCA will determine all PCs, we can only plot two at a time. In contrast, UMAP will take the information from any number of top PCs to arrange the cells in this multidimensional space. It will take those distances in multidimensional space, and try to plot them in two dimensions. In this way, the distances between cells represent similarity in expression.

To generate these visualizations we need to first run PCA and UMAP methods. Let's start with PCA.

```r
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")  
```

<p align="center">
<img src="../img/integrated_pca.png" width="600">
</p>

We can see with the PCA mapping that we have a good overlay of both conditions by PCA. 

Now, we can also visualize with UMAP. Let's run the method and plot.

```r
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
			     reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)                             
```

<p align="center">
<img src="../img/integrated_umap.png" width="600">
</p>

Again, we see good alignment of the two conditions using both methods. Sometimes it's easier to see whether all of the cells align well if we **split the plotting between conditions**, which we can do by adding the `split.by` argument to the `DimPlot()` function:

```r
DimPlot(seurat_integrated,
        split.by = "sample")  
```

<p align="center">
<img src="../img/SC_umap_split_int.png" width="600">
</p>

> When we compare to the unintegrated dataset, it is clear that this dataset benefitted from the integration!
> 
> <p align="center">
> <img src="../img/unintegrated_umap.png" width="400">
> </p>


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*
