# Answer key -  SCT, regression and normalization 

**1. First, turn the mitochondrial ratio variable into a new categorical variable based on quartiles (using the code below)::**

```
# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                   breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                   labels=c("Low","Medium","Medium high", "High"))
```

**2. Next, plot the PCA similar to how we did with cell cycle regression.** _Hint: use the new mitoFr variable to split cells and color them accordingly._

```
# Plot the PCA colored by mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")
```


**3. Evaluate the PCA plot generated in #2.**

<p align="center">
<img src="../img/pre_mitoFr_pca.png" width="600">
</p>

  
  **b. Describe what you see.** 
  
Based on this plot, we can see that there is a different pattern of scatter for the plot containing cells with "High" mitochondrial expression. We observe that the lobe of cells on the left-hand side of the plot is where most of the cells with high mitochondrial expression are. For all other levels of mitochondrial expression we see a more even distribution of cells across the PCA plot. 

  
  **c. Would you regress out mitochndrial fraction as a source of unwanted variation?**
  
  Since we see this clear difference, we will regress out the 'mitoRatio' when we identify the most variant genes.
  
  
**4. Are the same assays available for the "stim" samples within the `split_seurat` object? What is the code you used to check that?**

Yes they are available. The code use is:

```
split_seurat$stim@assays
```

**5. Any observations for the genes or features listed under *"First 10 features:"* and the *"Top 10 variable features:"* for "ctrl" versus "stim"?**

For the first 10 features, it appears that the same genes are present in both "ctrl" and "stim"

For the top 10 variable features, these are different in the the 2 conditions with some overlap between them.
