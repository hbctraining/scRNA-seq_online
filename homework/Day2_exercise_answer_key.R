#### Clustering
# After loading seurat_integrated.RData, set the resolution to 0.4, and plot the UMAP. How many clusters are present in our data?
# Answer: 13

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
