#### QC

## 1. Extract the new metadata from the filtered Seurat object using the code provided below:

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

## 2. Perform all of the same QC plots using the filtered data.


## 3. Report the number of cells left for each sample, and comment on whether the number of cells removed is high or low. Can you give reasons why this number is still not ~12K (which is how many cells were loaded for the experiment)?




## 4. After filtering for nGene per cell, you should still observe a small shoulder to the right of the main peak. What might this shoulder represent?




## 5. When plotting the nGene against nUMI do you observe any data points in the bottom right quadrant of the plot? What can you say about these cells that have been removed?



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
