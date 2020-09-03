# Schedule for the single-cell RNA-seq data analysis workshop

## Pre-reading

* [Introduction to scRNA-seq](../lessons/01_intro_to_scRNA-seq.md)
* [Raw data to count matrix](../lessons/02_SC_generation_of_count_matrix.md)
* [Download this project](https://www.dropbox.com/sh/pgyjaf507mxixuy/AABMhBOcqyoY147aOh3gx4HTa?dl=1)

## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 09:45 | [Workshop introduction](../slides/Intro_to_workshop.pdf) | Radhika |
| 09:45 - 11:00| [Single-cell RNA-seq design and methods](https://www.dropbox.com/s/c3z0l5mvnfis3xd/SCC_HBC_082020_final.pdf?dl=1) | [Dr. Luciano Martelotto](https://singlecellcore.hms.harvard.edu/people/luciano-martelotto) |
| 11:00 - 11:20 | scRNA-seq pre-reading discussion | Meeta |
| 11:20 - 12:00 | [Quality control set-up](../lessons/03_SC_quality_control-setup.md) | Jihe |


### Self Learning #1

* [Quality control](../lessons/04_SC_quality_control.md)
* [Overview of Clustering Workflow](../lessons/postQC_workflow.md)
* [Theory of Normalization and PCA](../lessons/05_normalization_and_PCA.md)
* [Normalization and regressing out unwanted variation](../lessons/06_SC_SCT_normalization.md)

### Assignment #1
* All exercises from above lessons have been put together in [R script format](../homework/Day1_exercise.R) (download for local access).
* Add your solutions to the exercises in the downloaded `.R` file and **upload the saved file** to [Dropbox](https://www.dropbox.com/request/nQ6oZ63K0GOl8j4T5Mqe) **day before the next class**.
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* **[Answer Key](../homework/Day1_exercise_answer_key.R)**

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:45 | Self-learning lessons discussion | Meeta/Jihe |
| 10:45 - 12:00| [Integration](../lessons/06_integration.md) | Radhika |

### Self Learning #2
* [Clustering](../lessons/07_SC_clustering_cells_SCT.md)
* [Clustering quality control](../lessons/08_SC_clustering_quality_control.md)
* [Marker identification](../lessons/09_merged_SC_marker_identification.md)

### Assignment #2
* All exercises from above lessons have been put together in [R script format](../homework/Day2_exercise.R) (download for local access).
* Add your solutions to the exercises in the downloaded `.R` file and **upload the saved file** to [Dropbox](https://www.dropbox.com/request/bVrUuHc3AaRZOe22CXvH) **day before the next class**.
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* **[Answer Key](../homework/Day2_exercise_answer_key.R)**

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:30 | Self-learning lessons discussion | All |
| 10:30 - 10:45 | Workflow summary | Meeta |
| 10:45 - 11:30 | Discussion, Final Q & A | All |
| 11:30 - 12:00 | [Wrap up](../slides/Workshop_wrapup.pdf) | Radhika |

***

## Downstream analysis tutorials

* [scRNA-seq Differential Expression Analysis between Conditions with DESeq2](../lessons/pseudobulk_DESeq2_scrnaseq.md)

## Resources
We have covered the analysis steps in quite a bit of detail for scRNA-seq exploration of cellular heterogeneity using the Seurat package. For more information on topics covered, we encourage you to take a look at the following resources:

* [Seurat vignettes](https://satijalab.org/seurat/vignettes.html)
* [Seurat cheatsheet](https://satijalab.org/seurat/essential_commands.html)
* ["Principal Component Analysis (PCA) clearly explained"](https://www.youtube.com/watch?v=_UVHneBUBW0), a video from [Josh Starmer](https://twitter.com/joshuastarmer)
* [Additional information about cell cycle scoring](../lessons/cell_cycle_scoring.md)
* [Using R on the O2 cluster](https://hbctraining.github.io/Intro-to-Unix-QMB/lessons/R_on_o2.html)
* Highlighted papers for sample processing steps (pre-sequencing):
  - ["Sampling time-dependent artifacts in single-cell genomics studies."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02032-0) *Massoni-Badosa et al.* 2019
  - ["Dissociation of solid tumor tissues with cold active protease for single-cell RNA-seq minimizes conserved collagenase-associated stress responses."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1830-0) *O'Flanagan et al.* 2020
  - ["Systematic assessment of tissue dissociation and storage biases in single-cell and single-nucleus RNA-seq workflows."](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02048-6) *Denisenko et al.* 2020


## Building on this workshop

* Other online scRNA-seq courses:
  - [https://hemberg-lab.github.io/scRNA.seq.course/](https://hemberg-lab.github.io/scRNA.seq.course/)
  - [https://github.com/SingleCellTranscriptomics](https://github.com/SingleCellTranscriptomics)

* Resources for scRNA-seq Sample Prep:
  - [https://www.protocols.io/](https://www.protocols.io/)
  - [https://support.10xgenomics.com/single-cell-gene-expression/sample-prep](https://support.10xgenomics.com/single-cell-gene-expression/sample-prep)
  - [https://community.10xgenomics.com/](https://community.10xgenomics.com/)

## Other helpful links
* [Online learning resources](https://hbctraining.github.io/bioinformatics_online/lists/online_trainings.html)
* [All hbctraining training materials](https://hbctraining.github.io/main)


