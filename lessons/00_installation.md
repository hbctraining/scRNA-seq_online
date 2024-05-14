# Single-cell RNA-seq data analysis workshop 

### Installation Requirements

#### Applications
Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) **(version 4.0.0 or above)**
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Packages for R

> **Note 1: Install the packages in the order listed below.**

> **Note 2:  All the package names listed below are case sensitive!**
 
> **Note 3**: At any point (especially if you’ve used R/Bioconductor in the past), in the console **R may ask you if you want to update any old packages by asking Update all/some/none? [a/s/n]:**. If you see this, **type "a" at the prompt and hit Enter** to update any old packages. _Updating packages can sometimes take quite a bit of time to run, so please account for that before you start with these installations._  

> **Note 4:** If you see a message in your console along the lines of “binary version available but the source version is later”, followed by a question, **“Do you want to install from sources the package which needs compilation? y/n”, type n for no, and hit enter**.


**(1)** Install the 4 packages listed below from **Bioconductor** using the the `BiocManager::install()` function.

1. `AnnotationHub`
1. `ensembldb`
1. `multtest`
1. `glmGamPoi`

**Please install them one-by-one as follows:**

```r
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
& so on ...
```

**(2)** Install the 8 packages listed below from **CRAN** using the `install.packages()` function. 

1. `tidyverse`
1. `Matrix`
1. `RCurl`
1. `scales`
1. `cowplot`
1. `BiocManager`
1. `Seurat`
1. `metap`

**Please install them one-by-one as follows:**

```r
install.packages("tidyverse")
install.packages("Matrix")
install.packages("RCurl")
& so on ...
```

**(3)** Finally, please check that all the packages were installed successfully by **loading them one at a time** using the `library()` function.  

```r
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(metap)
library(AnnotationHub)
library(ensembldb)
library(multtest)
library(glmGamPoi)
```

**(4)** Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```

---
