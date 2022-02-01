#### QC

## 1. Extract the new metadata from the filtered Seurat object using the code provided below:

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

## 2. Perform all of the same QC plots using the filtered data.


## 3. Report the number of cells left for each sample, and comment on whether the number of cells removed is high or low. Can you give reasons why this number is still not ~12K (which is how many cells were loaded for the experiment)?

The number of cells drops from just over 15K down to just udner 15K. We did not lose a lot of cells during the filtering process. We still don't have 12K cells. The number of cells we have can be higher due to:

- barcoded beads in the droplet with no actual cell present
- more than one barcode in the droplet; results in a single cell's mRNA represented as two cells
- dead or dying cells encapsulated into the droplet 


## 4. After filtering for nGene per cell, you should still observe a small shoulder to the right of the main peak. What might this shoulder represent?

The shoulder to the right represents a population of cells that have a higher number of genes identified. This set of cells could represent a biologically different type of cells (i.e. quiescent cell populations, less complex cells of interest), but it could also just mean that they are a set of bigger celltypes (i.e. cells with high counts may be cells that are larger in size). This threshold should be assessed with other metrics to see if it is necessary to distinguis this set of cells.


## 5. When plotting the nGene against nUMI do you observe any data points in the bottom right quadrant of the plot? What can you say about these cells that have been removed?

The data points in the bottom right quadrant have been removed due to the nGene threhold applied uring filtering. Those data points represent cells with a high number of UMIs but only a few number of genes. These are cells in which a small set of of genes were sequenced over and over again. These low complexity cells could represent a specific cell type (i.e. red blood cells which lack a typical transcriptome). But this could also represent dying cells, or could be due to some other strange artifact or contamination.


#### Exploring sources of unwanted variation

## 1. First, turn the mitochondrial ratio variable into a new categorical variable based on quartiles (using the code below):

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                   breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                   labels=c("Low","Medium","Medium high", "High"))

## 2. Next, plot the PCA similar to how we did with cell cycle regression. Hint: use the new mitoFr variable to split cells and color them accordingly.


## 3. 

## a. Evaluate the PCA plot generated in #2.


## b. Describe what you see.


## c. Would you regress out mitochndrial fraction as a source of unwanted variation?


#### SCTransform

## 1. Are the same assays available for the "stim" samples within the split_seurat object? What is the code you used to check that?

## 2. Any observations for the genes or features listed under "First 10 features:" and the "Top 10 variable features:" for "ctrl" versus "stim"?
