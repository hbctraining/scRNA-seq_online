---
title: "Single-cell RNA-seq: Marker identification"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Monday, July 15, 2019
---

Approximate time: 45 minutes

## Learning Objectives:

* Understand how to determine markers of individual clusters
* Understand the iterative processes of clustering and marker identification

# Single-cell RNA-seq marker identification

Now that we have identified our desired clusters, we can move on to marker identification, which will allow us to verify the identity of certain clusters and help surmise the identity of any unknown clusters. 

<img src="../img/sc_workflow.png" width="800">

***

_**Goals:**_ 
 
 - _To **determine the gene markers** for each of the clusters_
 - _To **identify cell types** of each cluster using markers_
 - _To determine whether there's a need to **re-cluster based on cell type markers**, perhaps clusters need to be merged or split_

_**Challenges:**_
 
 - _Over-interpretation of the results_
 - _Combining different types of marker identification_

_**Recommendations:**_
 
 - _Think of the results as hypotheses that need verification. Inflated p-values can lead to over-interpretation of results (essentially each cell is used as a replicate). Top markers are most trustworthy._
 _Identify all markers conserved between conditions for each cluster_
 - _Identify markers that are differentially expressed between specific clusters_

***

Our clustering analysis resulted in the following clusters:

<p align="center">
<img src="../img/SC_umap_loadimg.png" width="800">
</p>

Remember that we had the following questions from the clustering analysis:

1. *What is the cell type identity of cluster 7?*
2. *Is cluster 6 a CD8+ T cell or an NK cell?* *Is cluster 13 a T cell or an NK cell?*
3. *Do the clusters corresponding to the same cell types have biologically meaningful differences? Are there subpopulations of these cell types?*
4. *Can we acquire higher confidence in these cell type identities by identifying other marker genes for these clusters?*

There are a few different types of marker identification that we can explore using Seurat to get to the answer of these questions. Each with their own benefits and drawbacks:

1. **Identification of all markers for each cluster:** this analysis compares each cluster against all others and outputs the genes that are differentially expressed/present. 
	- *Useful for identifying unkown clusters and improving confidence in hypothesized cell types.*

2. **Identification of conserved markers for each cluster regardless of condition:** This analysis looks for those genes that are conserved in the cluster across all conditions. This analysis will output genes that are consistently differentially expressed/present for all of the sample groups. These genes can help to figure out the identity for the cluster. Often, this analysis is performed only for those clusters whose identity is uncertain or novel.
	- *Useful with more than one condition to identify cell type markers that are conserved across conditions.*  	
3. **Marker identification between specific clusters:** this analysis explores differentially expressed genes between specific clusters. 
	- *Useful for determining differences in gene expression between clusters with markers that are similar in the above analyses.*

Since we only are exploring a single sample, we do not need to explore conserved markers, but we will investigate all markers for each cluster and tease out markers differentiating specific clusters.

## Identification of all markers for each cluster

For this analysis we are comparing each cluster against all other clusters to identify cluster markers using the ` FindAllMarkers()` function. This function has two important arguments which provide thresholds for determining whether a gene is a marker:

- `logfc.threshold`: minimum log2 foldchange for average expression of gene in cluster relative to the average expression in all other clusters combined
	- **Cons:** 
		- could miss those cell markers that are expressed in the cluster being compared, but not in the other clusters, if the average log2FC doesn't meet the threshold
		- could return a lot of metabolic/ribosomal genes due to slight differences in metabolic output by different cell types, which are not as useful to distinguish cell type identities
- `min.diff.pct`: minimum percent difference between the percent of cells expressing the gene in the cluster and the percent of cells expressing gene in all other clusters combined
	- **Cons:** could miss those cell markers that are expressed in all cells, but are highly up-regulated in this specific cell type
	
You could use one or the other of these arguments or both. We will be a bit lenient and use only the log2 fold change threshold greater than 0.25. We will also specify to return only the positive markers for each cluster

> **NOTE:** This command can quite take long to run (~15 minutes). In the interest of time and for us to get through materials in a timely manner, we ask that you **download the data** using the [link provided](https://raw.githubusercontent.com/hbctraining/scRNA-seq/master/data/control_all_markers.csv). *Right-click and "Save as.." into your `results` folder.* Now you can load this data into your R environment:
> 
> `markers <- read.csv("results/control_all_markers.csv")`
>


```r
## DO NOT RUN THIS CODE ##
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_control, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)
                          
View(markers)                          
```

The order of the columns doesn't seem the most intuitive, so we will reorder the columns with the `cluster` first followed by the `gene`.

```r
# Combine markers with gene descriptions 
ann_markers <- inner_join(x = markers, 
                          y = annotations[, c("gene_name", "description")],
                          by = c("gene" = "gene_name")) %>%
        unique()

# Rearrange the columns to be more intuitive
ann_markers <- ann_markers[ , c(6, 7, 2:4, 1, 5,8)]

# Order the rows by p-adjusted values
ann_markers <- ann_markers %>%
        dplyr::arrange(cluster, p_val_adj)

View(ann_markers)
```

<p align="center">
<img src="../img/marker_table_loadObj.png" width="800">
</p>

**Usually the top markers are relatively trustworthy, but because of inflated p-values, many of the less significant genes are not so trustworthy as markers.**

When looking at the output, we suggest looking for markers with large differences in expression between `pct.1` and `pct.2` and larger fold changes. For instance if `pct.1` = 0.90 and `pct.2` = 0.80, I might not be as excited about that marker. However, if `pct.2` = 0.1 instead, then I would be much more excited about it. Also, I look for the majority of cells expressing marker in my cluster of interest. If `pct.1` is low, such as 0.3, I again might not be as interested in it.

- **cluster:** number corresponding to cluster
- **gene:** gene id
- **avg_logFC:** average log2 fold change. Positive values indicate that the gene is more highly expressed in the cluster.
- **pct.1**: The percentage of cells where the gene is detected in the cluster
- **pct.2**: The percentage of cells where the gene is detected on average in the other clusters
- **p_val:** p-value not adjusted for multiple test correction
- **p_val\_adj:** Adjusted p-value, based on bonferroni correction using all genes in the dataset, used to determine significance


If the format looks good, we can save our marker analysis results to file.

```r
# Save markers to file
write.csv(ann_markers, 
          file = "results/control_all_markers_annotated.csv", 
          quote = FALSE, 
          row.names = FALSE)
```

We can also output the top 5 markers by log2 fold change for each cluster for a quick perusal.

```r
# Extract top 5 markers per cluster
top5 <- ann_markers %>% 
        group_by(cluster) %>% 
        top_n(n = 5, 
              wt = avg_logFC)

# Visualize top 5 markers per cluster
View(top5)

```

<p align="center">
<img src="../img/top5_markers_loadObj.png" width="800">
</p>

Based on these marker results, we can determine whether the markers make sense for our hypothesized identities for each cluster:

| Cell Type | Clusters |
|:---:|:---:|
| CD14+ monocytes | 0, 15 | 
| FCGR3A+ monocytes | 8 |
| Conventional dendritic cells | 12 |
| Plasmacytoid dendritic cells | 12 |
| B cells | 4, 11 |
| T cells | 1, 2, 3, 6, 9, 10, 13, 14 |
| CD4+ T cells | 1, 2, 3, 9, 10, 13, 14 |
| CD8+ T cells| 6 |
| NK cells | 5,6, 13 |
| Megakaryocytes | 10 |
| Erythrocytes | - |
| Unknown | 7 |

If there were any questions about the identity of any clusters, exploring the cluster's markers would be the first step. Let's look at the `ann_markers`, filtering for cluster 7 and see if we can **obtain any hints about our unknown cluster**.

<p align="center">
<img src="../img/cluster7_markers_loadObj.png" width="800">
</p>

We see a lot of heat shock and DNA damage genes appear. Based on these markers, it is likely that these are **stressed or dying cells**. However, we also see T cell-associated genes and markers of activation. It is possible that these could be activated (cytotoxic) T cells. We could explore the quality metrics for these cells in more detail before removing just to support that argument.

We also had questions regarding the identity of cluster 6. Is cluster 6 a CD8+ T cell, an NK cell, or an NK T cell?

We can look at the markers of **cluster 6 to try to resolve the identity**:

<p align="center">
<img src="../img/cluster6_markers_loadObj.png" width="800">
</p>

There are definitely T cell receptors that are enriched among cluster 6. Since NK cells cannot have expression of the T cell receptor genes we can therefore conclude that these cannot be NK cells. On the other hand CD8+ T cells can have expression of killer cell receptors. So, could these be NK T cells? Probably not, since NK T cells are usually a rare population and in our case we have many cells here. Thus, we **hypothesize that cluster 6 represents activated CD8+ T cells (cytotoxic T cells)**.

To get a better idea of cell type identity we can explore the expression of different identified markers by cluster using the `FeaturePlot()` function. For example, we can look at the cluster 6 markers:

```r
# Plot top 5 markers for cluster 6
FeaturePlot(object = seurat_control, 
            features = top5[top5$cluster == 6, "gene"] %>%
                    pull(gene))
```

<p align="center">
<img src="../img/fig_cluster6_loadObj.png" width="800">
</p>

We can also explore the range in expression of specific markers by using violin plots:

> **Violin plots** are similar to box plots, except that they also show the probability density of the data at different values, usually smoothed by a kernel density estimator. A violin plot is more informative than a plain box plot. While a box plot only shows summary statistics such as mean/median and interquartile ranges, the violin plot shows the full distribution of the data. The difference is particularly useful when the data distribution is multimodal (more than one peak). In this case a violin plot shows the presence of different peaks, their position and relative amplitude.

```r
# Vln plot - cluster 6
VlnPlot(object = seurat_control, 
        features = top5[top5$cluster == 6, "gene"] %>%
                    pull(gene))
```        

<p align="center">
<img src="../img/fig_cluster6_loadObj_violin.png" width="800">
</p>

These results and plots can help us determine the identity of these clusters or verify what we hypothesize the identity to be after exploring the canonical markers of expected cell types previously.

***

**[Exercises](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_marker_analysis.html)**

1. **Download [the file](https://raw.githubusercontent.com/hbctraining/scRNA-seq/master/data/stim_all_markers.csv) contaning the results from `FindAllMarkers()` on the stimulated sample**. You will want to right-click and "Save as.." into your `results` folder. Read that file into your R session to a variable called `markers_stim`. 
2. **Annotate** the markers with gene descriptions.
3. **Reorder the columns** to be in the order shown below
4. **Arrange rows** by cluster, then by p-adjusted values
5. **Save** our rearranged marker analysis results to a file called `stim_all_markers.csv` in the `results` folder.
6. **Extract** the top 5 markers by log2 fold change for each cluster.
7. **Visualize** top 5 markers for cluster 0 using the `FeaturePlot()` and `VlnPlot()` functions.

***

## Identifying gene markers for each cluster

The last set of questions we had regarding the analysis involved whether the clusters corresponding to the same cell types have biologically meaningful differences. Sometimes the list of markers returned don't sufficiently separate some of the clusters. For instance, we had previously identified clusters 0, and 15 as CD14+ monocytes, but are there biologically relevant differences between these clusters of cells? We can use the `FindMarkers()` function to determine the genes that are differentially expressed between two specific clusters. 

```r
# Determine differentiating markers for CD14+ monocytes - clusters 0 versus 15
cd14_monos <- FindMarkers(seurat_control,
                          ident.1 = 0,
                          ident.2 = 15)                     

# Add gene symbols to the DE table
cd14_monos_markers <- cd14_monos %>%
        rownames_to_column("gene") %>%
        inner_join(y = annotations[, c("gene_name", "description")],
                   by = c("gene" = "gene_name")) %>%
        unique()

# Reorder columns and sort by log2 fold change        
cd14_monos_markers <- cd14_monos_markers[, c(1, 3:5,2,6:7)]

cd14_monos_markers <- cd14_monos_markers %>%
        dplyr::arrange(avg_logFC)
        
# View data
View(cd14_monos_markers)
```

<p align="center">
<img src="../img/sc_mono14_markers.png" width="800">
</p>

When looking through the results, we see quite a few T cell-specific markers, such as the T cell marker, CD3D, and T cell receptor genes. We also see lower expression of the CD14 and LYZ monocyte cell markers. It's possible that this cluster could represent doublets of CD14+ monocytes and T cells. 

We are not going to explore these genes in more depth, although, you would probably want to explore the expression of these genes in more depth visually using feature plots and violin plots before deciding on a label.

We would also like to determine how the CD4+ T cell clusters are different from each other. We could again explore these differences:

```r
# Determine differentiating markers for CD4+ T cell
cd4_tcells <- FindMarkers(seurat_control,
                          ident.1 = 1,
                          ident.2 = c(2, 3, 9, 10, 13, 14),
                          only.pos = TRUE)                     

# Add gene symbols to the DE table
cd4_tcells <- cd4_tcells %>%
        rownames_to_column("gene") %>%
        inner_join(y = annotations[, c("gene_name", "description")],
                   by = c("gene" = "gene_name")) %>%
        unique()

# Reorder columns and sort by log2 fold change        
cd4_tcells <- cd4_tcells[, c(1, 3:5,2,6:7)]

cd4_tcells <- cd4_tcells %>%
        dplyr::arrange(dplyr::desc(abs(avg_logFC))) 

# View data
View(cd4_tcells)

```

<p align="center">
<img src="../img/sc_cd4t_markers.png" width="800">
</p>

Of these top genes the **CREM gene** stands out as a marker of activation. We know that another marker of activation is CD69, and markers of naive or memory cells include the SELL and CCR7 genes. Let's explore activation status a bit visually using these new cell state markers:

| Cell State | Marker |
|:---:|:---:|
| Naive T cells | CCR7, SELL | 
| Activated T cells | CREM, CD69 |

```r
# Plot gene markers of activated and naive/memory T cells
FeaturePlot(seurat_control, 
            reduction = "umap", 
            features = c("CREM", "CD69", "CCR7", "SELL"))
```

<p align="center">
<img src="../img/sc_cd4t_act-mem_markers.png" width="800">
</p>

The activated CD4+ T cells correspond to clusters 1 and 9, while the naive or memory CD4+ T cells represent clusters 2, 3, and 14.

Now taking all of this information, we can surmise the cell types of the different clusters and plot the cells with cell type labels.


| Cluster ID	| Cell Type |
|:-----:|:-----:|
|0	| CD14+ Monocytes|
|1	| Activated CD4+ T cells |
|2	| Naive or memory CD4+ T cells|
|3	| Naive or memory CD4+ T cells|
|4	| B cells |
|5	| NK cells |
|6	| CD8+ T cells |
|7	| Stressed / dying cells |
|8	| FCGR3A+ monocytes |
|9	| Activated CD4+ T cells |
|10	| Megakaryocytes |
|11	| B cells |
|12	| Dendritic cells |
|13	| NK cells |
|14	| Naive or memory CD4+ T cells |
|15| CD14+ monocytes / T cell doublets |


We can then reassign the identity of the clusters to these cell types:

```r
# Rename all identities
seurat_control <- RenameIdents(object = seurat_control, 
                               "0" = "CD14+ monocytes",
                               "1" = "Activated CD4+ T cells",
                               "2" = "Naive or memory CD4+ T cells",
                               "3" = "Naive or memory CD4+ T cells",
                               "4" = "B cells",
                               "5" = "NK cells",
                               "6" = "CD8+ T cells",
                               "7" = "Stressed / dying cells",
                               "8" = "FCGR3A+ monocytes",
                               "9" = "Activated CD4+ T cells",
                               "10" = "Megakaryocytes",
                               "11" = "B cells",
                               "12" = "Dendritic cells",
                               "13" = "NK cells",
                               "14" = "Naive or memory CD4+ T cells",
                               "15" = "CD14+ monocytes / T cell doublets")


# Plot the UMAP
DimPlot(object = seurat_control, 
        reduction = "umap", 
        label = TRUE,
        label.size = 6,
        repel = TRUE)
```

<p align="center">
<img src="../img/umap_labelled_subset_loadObj.png" width="800">
</p>

If we wanted to remove the stressed cells, we could use the `subset()` function:

```r
# Remove the stressed or dying cells
control_labelled <- subset(seurat_control,
                               idents = "Stressed / dying cells", invert = TRUE)

# Re-visualize the clusters
DimPlot(object = control_labelled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 6)
```

<p align="center">
<img src="../img/umap_control_labelled_subset_loadObj.png" width="800">
</p>

Now we would want to save our final labelled Seurat object:

```r        
# Save final R object
write_rds(control_labelled,
          file = "results/seurat_control_labelled.rds")       
```

We have completed the first round of the analysis for the `control` sample. Based on our results, we may need to proceed back to prior steps to optimize our parameters:

- Clusters were not separated well enough (i.e. CD4+ and CD8+ T cells were in same cluster):  adjust the clustering resolution or number of PCs used.
- Too many low quality cells observed: go back and perform more stringent QC by filtering out the low quality cells. 

If we validated the stressed/dying cells were indeed low quality, we may want to remove those cells in the QC and re-cluster. However, we are going to just leave those cells in and move on to include the `stimulated` sample.

[Click here for next lesson](08_SC_clustering_analysis_integration.md)

***


*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*
