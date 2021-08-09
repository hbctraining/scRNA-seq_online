#### Clustering

# After loading seurat_integrated.RData, assign the identity of the clusters with different resolution (0.4, 0.6, 0.8, 1.0, 1.4), and plot the corresponding UMAP. How many clusters are present for each resolution? Which resolution do you think makes sense?

# Resolution 0.4: 14 clusters
# Resolution 0.6: 17 clusters
# Resolution 0.8: 21 clusters
# Resolution 1.0: 22 clusters
# Resolution 1.4: 27 clusters

#### Clustering quality control
# Answer:
# Hypothesize the clusters corresponding to each of the different clusters in the table:
# Cell Type	Clusters
# CD14+ monocytes	1, 3, 14
# FCGR3A+ monocytes	9
# Conventional dendritic cells	15
# Plasmacytoid dendritic cells	19
# Marcrophages	-
# B cells	6, 11, 17
# T cells	0, 2, 4, 5, 10, 13, 18
# CD4+ T cells	0, 2, 4, 10, 18
# CD8+ T cells	5, 13
# NK cells	8, 12
# Megakaryocytes	16
# Erythrocytes -
# Unknown	7, 20

#### Marker identification
# In the previous lesson, we identified cluster 9 as FCGR3A+ monocytes by inspecting the expression of known cell markers FCGR3A and MS4A7. Use FindConservedMarkers() function to find conserved markers for cluster 9. What do you observe? Do you see FCGR3A and MS4A7 as highly expressed genes in cluster 9?
# Answer: Yes, FCGR3A and MS4A7 are among the top highly expressed genes (2nd and 3rd in our analysis) when performing FindConservedMarkers() for cluster 9. This observation confirms cluster 9 as FCGR3A+ monocytes.
