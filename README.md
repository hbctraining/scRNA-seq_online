# Single-cell RNA-seq analysis workshop 

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R/) | 2-day workshop (~10 hours of trainer-led time)|

### Description

This repository has teaching materials for a **2-day**, hands-on **Introduction to single-cell RNA-seq analysis** workshop. This 2-day hands-on workshop will instruct participants on how to design a single-cell RNA-seq experiment, and how to efficiently manage and analyze the data starting from count matrices. This will be a hands-on workshop in which we will focus on using the Seurat package using R/RStudio. Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/). 


### Learning Objectives

- Undertand the considerations when designing a single-cell RNA-seq experiment
- Discuss the steps involved in taking raw single-cell RNA-sequencing data and generating a count (gene expression) matrix
- Compute and assess QC metrics at every step in the workflow
- Cluster cells based on expression data and derive the identity of the different cell types present
- Perform integration of different sample conditions

> These materials are developed for a trainer-led workshop, but also amenable to self-guided learning.

### Lessons

Below are links to the lessons and suggested schedules:

* [Click here for schedule](schedule)


### Installation Requirements

#### Applications
Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) **(version 3.6.0 or above)**
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Packages for R

> **Note 1: Install the packages in the order listed below.**

> **Note 2:  When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable.**
 
> **Note 3: All the package names listed below are case sensitive!**

**(1)** Install the 7 packages listed below from **CRAN** using the `install.packages()` function. 

1. `tidyverse`
1. `Matrix`
1. `RCurl`
1. `scales`
1. `cowplot`
1. `BiocManager`
1. `Seurat`**

**Please install them one-by-one as follows:**

```r
install.packages("tidyverse")
install.packages("Matrix")
install.packages("RCurl")
& so on ...
```

> ** If you have trouble installing `Seurat`, please install `multtest` using the following lines of code, then try installing `Seurat` again:
>
> `install.packages("BiocManager")`
>
> `BiocManager::install("multtest")`

**(2)** Install the 3 packages listed below from **Bioconductor** using the the `BiocManager::install()` function.

1. `SingleCellExperiment`
1. `AnnotationHub`
1. `ensembldb`

**Please install them one-by-one as follows:**

```r
BiocManager::install("SingleCellExperiment")
BiocManager::install("AnnotationHub")
& so on ...
```

**(3)** Finally, please check that all the packages were installed successfully by loading them **one at a time** using the `library()` function.  

```r
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(SingleCellExperiment)
library(AnnotationHub)
library(ensembldb)
```

**(4)** Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```

****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
