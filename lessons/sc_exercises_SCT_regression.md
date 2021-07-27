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
  
  **c. Would you regress out mitochndrial fraction as a source of unwanted variation?**
