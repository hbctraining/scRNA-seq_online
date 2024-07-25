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


<details> 
<summary><b>Click here for details on Implementing Harmony within the Seurat workflow</b></summary>
In practice, we can easily use Harmony within our Seurat workflow. To perform integration, <code>Harmony</code> takes as input a <i>merged</i> Seurat object, containing data that has been appropriately normalized (i.e. here, normalized using <code>SCTransform</code>) and for which highly variable features and PCs are defined.

There are <b>2 ways to create the input</b>:

<ol><li><b>Merge the <i>raw</i> Seurat objects</b> for all samples to integrate; then perform normalization, variable feature selection and PC calculation on this merged object (workflow recommended by <code>Harmony</code> developers)<br></li>

<li>Perform (SCT) normalization independently on each sample and find integration features across samples using Seurat; then <b>merge these normalized Seurat objects</b>, set variable features manually to integration features, and finally calculate PCs on this merged object (workflow best reflecting recommendations for application of <code>SCTransform</code>)</li></ol><br>

In the first scenario, assuming <code>raw_seurat_list</code> is a list of N samples containing raw data that have only undergone QC filtering, we would thus run the following code:

<pre>
# Merge raw samples
merged_seurat <- merge(x = raw_seurat_list[[1]],
		       y = raw_seurat_list[2:length(raw_seurat_list)],
		       merge.data = TRUE)<br>

# Perform log-normalization and feature selection, as well as SCT normalization on global object
merged_seurat <- merged_seurat %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData() %>%
    SCTransform(vars.to.regress = c("mitoRatio"))<br>

# Calculate PCs using variable features determined by SCTransform (3000 by default)
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)
</pre>

In the second scenario, assuming <code>norm_seurat_list</code> is a list of N samples similar to our <code>split_seurat</code> object, i.e. containing data that have been normalized as demonstrated in the previous lecture on SCT normalization, we would thus run the following code:

<pre>
# Find most variable features across samples to integrate
integ_features <- SelectIntegrationFeatures(object.list = norm_seurat_list, nfeatures = 3000)<br>

# Merge normalized samples
merged_seurat <- merge(x = norm_seurat_list[[1]],
		       y = norm_seurat_list[2:length(raw_seurat_list)],
		       merge.data = TRUE)
DefaultAssay(merged_seurat) <- "SCT"<br>

# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- integ_features<br>

# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)
</pre>

<blockquote><i><b>NOTE:</b> As mentioned above, there is active discussion within the community regarding which of those 2 approaches to use (see for example <a href="https://github.com/immunogenomics/harmony/issues/41">here</a> and <a href="https://github.com/satijalab/sctransform/issues/55#issuecomment-633843730">here</a>). We recommend that you check GitHub forums to make your own opinion and for updates.</i></blockquote>

Regardless of the approach, we now have a merged Seurat object containing normalized data for all the samples we need to integrate, as well as defined variable features and PCs. 

One last thing we need to do before running <code>Harmony</code> is to <b>make sure that the metadata of our Seurat object contains one (or several) variable(s) describing the factor(s) we want to integrate on</b> (e.g. one variable for <code>sample_id</code>, one variable for <code>experiment_date</code>). 

We're then ready to run <code>Harmony</code>!

<pre>
harmonized_seurat <- RunHarmony(merged_seurat, 
				group.by.vars = c("sample_id", "experiment_date"), 
				reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
</pre>

<blockquote><i><b>NOTE:</b> You can specify however many variables to integrate on using the <code>group.by.vars</code> parameter, although we would recommend keeping these to the minimum necessary for your study.</i></blockquote>

The line of code above adds a new reduction of 50 "harmony components" (~ corrected PCs) to our Seurat object, stored in <code>harmonized_seurat@reductions$harmony</code>.

To make sure our </code>Harmony</code> integration is reflected in the data visualization, we still need to generate a UMAP derived from these harmony embeddings instead of PCs:

<pre>
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:40)
</pre>

Finally, when running the clustering analysis later on (see next lecture for details), we will also need to set the reduction to use as "harmony" (instead of "pca" by default).

<pre>
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
</pre>

The rest of the <code>Seurat</code> workflow and downstream analyses after integration using <code>Harmony</code> can then proceed without further amendments.
</details>

***


*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*
