# Schedule for the single-cell RNA-seq data analysis workshop

## Pre-reading

* [Introduction to scRNA-seq](../lessons/01_intro_to_scRNA-seq.md)
* [Raw data to count matrix](../lessons/02_SC_generation_of_count_matrix.md)
* [Download this project](https://www.dropbox.com/s/we1gmyb9c8jej2u/single_cell_rnaseq.zip?dl=1)

## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 09:45 | [Workshop introduction](../slides/Intro_to_workshop_all.pdf) | Meeta |
| 09:45 - 11:00| [Single-cell RNA-seq design and methods]() | [Dr. Mandovi Chatterjee](https://singlecellcore.hms.harvard.edu/people/mandovi-chatterjee-phd) |
| 11:00- 11:05 | Break |
| 11:05 - 11:15 | scRNA-seq pre-reading discussion | Meeta |
| 11:15 - 11:55 | [Quality control set-up](../lessons/03_SC_quality_control-setup.md) | Radhika |
| 11:55 - 12:00 | Overview of self-learning materials and homework submission | Meeta |


### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:
   1. [Quality control](../lessons/04_SC_quality_control.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Before you start any analysis, it’s important to know whether or not you have good quality cells. At these early stages you can flag or remove samples that could produce erroneous results downstream.
 <br><br>In this lesson you will:<br>
             - Compute essential QC metrics for each sample<br>
             - Plot the metrics to visualize the metrics per sample<br>
             - Critically evaluate each plot and learn what each QC metric means<br><br>
        </details>

1. Please **study the contents** and **work through all the code** within the following lessons:
   * [](../lessons/04_SC_quality_control.md)
   * [Overview of Clustering Workflow](../lessons/postQC_workflow.md)
   * [Theory of Normalization and PCA](../lessons/05_normalization_and_PCA.md)
   * [Normalization and regressing out unwanted variation](../lessons/06_SC_SCT_normalization.md)

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your R code from the exercises to this [(downloadable) R script](../homework/Day1_exercise.R)
   * **Upload the saved R script file** to [Dropbox](https://www.dropbox.com/request/kIqpFLAIDCix9eGZ6vDh) **day before the next class**.

3. **Run the code in this [script](https://github.com/hbctraining/scRNA-seq_online/raw/master/scripts/integration_code.R)** to perform the steps of integration. We will discuss the code and theory in class.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:40 | Self-learning lessons discussion | Meeta |
| 10:40 - 10:45 | Break |
| 10:45 - 12:00| [Integration](../lessons/06_integration.md) | Radhika |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:
   * [Clustering](../lessons/07_SC_clustering_cells_SCT.md)
   * [Clustering quality control](../lessons/08_SC_clustering_quality_control.md)
   * [Marker identification](../lessons/09_merged_SC_marker_identification.md)

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your R code from the exercises to this [(downloadable) R script](../homework/Day2_exercise.R)
   * **Upload the saved R script file** to [Dropbox](https://www.dropbox.com/request/gEmJDyFK0pY7hpMzPu0X) **day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:30 | Self-learning lessons discussion | All |
| 10:30 - 10:40 | Workflow summary | Radhika |
| 10:40 - 10:45 | Break |
| 10:45 - 11:30 | Discussion, Final Q & A | All |
| 11:30 - 12:00 | [Wrap up](../slides/Workshop_wrapup.pdf) | Meeta |

***

## Answer Keys

* **[Answer key - assignment #1](../homework/Day1_exercise_answer_key.md)**
* **[Answer key - assignment #2](../homework/Day2_exercise_answer_key.R)**

## Downstream analyses

[Differential expression between conditions](../lessons/pseudobulk_DESeq2_scrnaseq.md)

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
* [Azimuth reference-based analysis](https://azimuth.hubmapconsortium.org/)
* [CellMarker resource](http://biocc.hrbmu.edu.cn/CellMarker/)
* Highlighted papers for single-nuclei RNA-seq:
  - [Single-nucleus and single-cell transcriptomes compared in matched cortical cell types](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6306246/)
  - [A single-cell and single-nucleus RNA-Seq toolbox for fresh and frozen human tumors](https://www.nature.com/articles/s41591-020-0844-1)
* [Ligand-receptor analysis with CellphoneDB](https://www.nature.com/articles/s41576-020-00292-x)


## Building on this workshop

* Other online scRNA-seq courses:
  - [http://bioconductor.org/books/release/OSCA/](http://bioconductor.org/books/release/OSCA/)
  - [https://liulab-dfci.github.io/bioinfo-combio/](https://liulab-dfci.github.io/bioinfo-combio/)
  - [https://hemberg-lab.github.io/scRNA.seq.course/](https://www.singlecellcourse.org/)
  - [https://github.com/SingleCellTranscriptomics](https://github.com/SingleCellTranscriptomics)
  - [https://broadinstitute.github.io/2020_scWorkshop/](https://broadinstitute.github.io/2020_scWorkshop/)
  - [https://nbisweden.github.io/workshop-scRNAseq/schedule.html](https://nbisweden.github.io/workshop-scRNAseq/schedule.html)

* Resources for scRNA-seq Sample Prep:
  - [https://www.protocols.io/](https://www.protocols.io/)
  - [https://support.10xgenomics.com/single-cell-gene-expression/sample-prep](https://support.10xgenomics.com/single-cell-gene-expression/sample-prep)
  - [https://community.10xgenomics.com/](https://community.10xgenomics.com/)

## Other helpful links
* [Online learning resources](https://hbctraining.github.io/bioinformatics_online/lists/online_trainings.html)
* [All hbctraining training materials](https://hbctraining.github.io/main)


