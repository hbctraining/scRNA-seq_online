---
title: "Single-cell RNA-seq: Post-QC workflow"
author: "Mary Piper, Meeta Mistry, Radhika Khetani"
date: Tuesday, February 25th, 2020
---

Approximate time: 20 minutes

## Learning Objectives:

* Describe the workflow for single-cell RNA-seq analysis after the quality control step.


# Single-cell RNA-seq Clustering Workflow

Now that we have our high quality cells, we can move forward with the workflow. Ultimately, we want to cluster cells and identify different potential celltypes however there are a few steps to walk-through before we get there. **The green boxes in our workflow schematic below correspond to the steps taken post-QC and together consistute the clustering workflow.**

<p align="center">
<img src="../img/sc_workflow_2022.jpg" width="630">
</p>

## Clustering workflow

For something to be informative, it needs to exhibit variation, but not all variation is informative. The goal of our clustering analysis is to keep the major sources of variation in our dataset that should define our cell types, while restricting the variation due to uninteresting sources of variation (sequencing depth, cell cycle differences, mitochondrial expression, batch effects, etc.). Then, to determine the cell types present, we will perform a clustering analysis using the most variable genes to define the major sources of variation in the dataset. 

The workflow for this analysis is adapted from the following sources:

- Satija Lab: [Seurat v3 Guided Integration Tutorial](https://satijalab.org/seurat/v3.0/immune_alignment.html)
- Paul Hoffman: [Cell-Cycle Scoring and Regression](http://satijalab.org/seurat/cell_cycle_vignette.html)

To identify clusters, the following steps will be performed:


### 1. Explore sources of unwanted variation

The first step in the workflow is to see if our data contains any unwanted variability. The most common biological effect that is evaluated in single-cell RNA-seq data is the effect of cell cycle on the transcriptome. Another known biological effect is mitochondrial gene expression, which is interpreted as an indication of cell stress. This step of the workflow involves exploring our data to identify which covariates we would like to regress out. 

### 2. Normalization and regressing out sources of unwanted variation

Seurat recently introduced a new method called `sctransform` which performs multiple processing steps on scRNA-seq data. Normalization is required to scale the raw count data to obtain correct relative gene expression abundances between cells. The `sctransform` function implements an advanced normalization and variance stabilization of the data. The `sctransform` function also regresses out sources of unwanted variation in our data. In the previous step, we had identified these sources of variability, and here we specify what those covariates are. 

### 3. Integration

Often with single cell RNA-seq we are working with multiple samples which correspond to different sample groups, multiple experiments or different modalities. If we want to ultimately compare celltype expression between groups it is recommended to integrate the data. Integration is a powerful method that uses these shared sources of greatest variation to identify shared sub-populations across conditions or datasets [Stuart and Butler et al. (2018)]. There are several steps involved in performing intergration in Seurat. Once complete, we use visualization methods to ensure a good integration before we proceed to cluster cells.

> **NOTE:** Integration is optional. We recommend going through the workflow without integration to decide whether or not it is necessary for your data. 

### 4. Clustering cells

Clusters of cells are obtained by grouping cells based on the similarity of their gene expression profiles. Expression profile similarity is determined via distance metrics, which often take dimensionality‐reduced representations as input. Seurat assigns cells to clusters based on their PCA scores derived from the expression of the integrated most variable genes. 

### 5. Cluster quality evaluation

The clusters identified in our data represent groups of cells that presumably belong to a similar cell type. Before we can confirm the celltype of a group of member cells, the following steps are taken:

   * **a.** Check to see that clusters are not influenced by sources of uninteresting variation.
   * **b.** Check to see whether the major principal components are driving the different clusters.
   * **c.** Explore the cell type identities by looking at the expression for known markers across the clusters. 


***


*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*
