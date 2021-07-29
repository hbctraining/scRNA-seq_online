#### Clustering
# After loading seurat_integrated.RData, assign the identity of the clusters with different resolution (0.4, 0.6, 0.8, 1.0, 1.4), and plot the corresponding UMAP. How many clusters are present for each resolution? Which resolution do you think makes sense?

#### Clustering quality control
# Hypothesize the clusters corresponding to each of the different clusters in the table (fill in the question marks):
# Cell Type	Clusters
# CD14+ monocytes	1, 3, 14
# FCGR3A+ monocytes	9
# Conventional dendritic cells	15
# Plasmacytoid dendritic cells	19
# Marcrophages	-
# B cells	?
# T cells	?
# CD4+ T cells	?
# CD8+ T cells	?
# NK cells	?
# Megakaryocytes	?
# Erythrocytes ?
# Unknown	?

#### Marker identification
# In the previous lesson, we identified cluster 9 as FCGR3A+ monocytes by inspecting the expression of known cell markers FCGR3A and MS4A7. Use FindConservedMarkers() function to find conserved markers for cluster 9. What do you observe? Do you see FCGR3A and MS4A7 as highly expressed genes in cluster 9?
