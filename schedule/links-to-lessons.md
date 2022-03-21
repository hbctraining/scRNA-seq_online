# Single-cell RNA-seq data analysis workshop 

### Learning Objectives

- Understand the considerations when designing a single-cell RNA-seq experiment
- Discuss the steps involved in taking raw single-cell RNA-sequencing data and generating a count (gene expression) matrix
- Compute and assess QC metrics at every step in the workflow
- Cluster cells based on expression data and derive the identity of the different cell types present
- Perform integration of different sample conditions

## Installations

1. [Follow the instructions linked here](../README.md#installation-requirements) to download R and RStudio + Install Packages from CRAN and Bioconductor

1. [Download this project](https://www.dropbox.com/s/we1gmyb9c8jej2u/single_cell_rnaseq.zip?dl=1)

## Lessons

### Part 1
1. [Introduction to scRNA-seq](../lessons/01_intro_to_scRNA-seq.md)
1. [Raw data to count matrix](../lessons/02_SC_generation_of_count_matrix.md)

***

### Part II
1. [Quality control set-up](../lessons/03_SC_quality_control-setup.md)
1. [Quality control](../lessons/04_SC_quality_control.md)
1. [Overview of Clustering Workflow](../lessons/postQC_workflow.md)
1. [Theory of Normalization and PCA](../lessons/05_normalization_and_PCA.md)
1. [Normalization and regressing out unwanted variation](../lessons/06_SC_SCT_normalization.md)

* [Solution to exercises in above lessons](../homework/day1_hw_answer-key.R)
          
***

### Part III
1. [Integration](../lessons/06_integration.md)
1. [Clustering](../lessons/07_SC_clustering_cells_SCT.md)
1. [Clustering quality control](../lessons/08_SC_clustering_quality_control.md)
1. [Marker identification](../lessons/09_merged_SC_marker_identification.md)

* [Solution to exercises in above lessons](../homework/Day2_exercise_answer_key.R)

***

## Building on this workshop

* Downstream analysis
  - [Differential expression between conditions](../lessons/pseudobulk_DESeq2_scrnaseq.md)

* Other online scRNA-seq courses:
  - [http://bioconductor.org/books/release/OSCA/](http://bioconductor.org/books/release/OSCA/)
  - [https://liulab-dfci.github.io/bioinfo-combio/](https://liulab-dfci.github.io/bioinfo-combio/)
  - [https://hemberg-lab.github.io/scRNA.seq.course/](https://hemberg-lab.github.io/scRNA.seq.course/)
  - [https://github.com/SingleCellTranscriptomics](https://github.com/SingleCellTranscriptomics)
  - [https://broadinstitute.github.io/2020_scWorkshop/](https://broadinstitute.github.io/2020_scWorkshop/)

* Resources for scRNA-seq Sample Prep:
  - [https://www.protocols.io/](https://www.protocols.io/)
  - [https://support.10xgenomics.com/single-cell-gene-expression/sample-prep](https://support.10xgenomics.com/single-cell-gene-expression/sample-prep)
  - [https://community.10xgenomics.com/](https://community.10xgenomics.com/)

***

## Resources
We have covered the analysis steps in quite a bit of detail for scRNA-seq exploration of cellular heterogeneity using the Seurat package. For more information on topics covered, we encourage you to take a look at the following resources:

* [Seurat vignettes](https://satijalab.org/seurat/vignettes.html)
* [Seurat cheatsheet](https://satijalab.org/seurat/essential_commands.html)
* [Satija Lab: Single Cell Genomics Day](https://satijalab.org/scgd21/)
* ["Principal Component Analysis (PCA) clearly explained"](https://www.youtube.com/watch?v=_UVHneBUBW0), a video from [Josh Starmer](https://twitter.com/joshuastarmer)
* [Additional information about cell cycle scoring](../lessons/cell_cycle_scoring.md)
* [Using R on the O2 cluster](https://hbctraining.github.io/Intro-to-Unix-QMB/lessons/R_on_o2.html)
* Highlighted papers for sample processing steps (pre-sequencing):
  - ["Sampling time-dependent artifacts in single-cell genomics studies."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02032-0) *Massoni-Badosa et al.* 2019
  - ["Dissociation of solid tumor tissues with cold active protease for single-cell RNA-seq minimizes conserved collagenase-associated stress responses."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1830-0) *O'Flanagan et al.* 2020
  - ["Systematic assessment of tissue dissociation and storage biases in single-cell and single-nucleus RNA-seq workflows."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02048-6) *Denisenko et al.* 2020

****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
