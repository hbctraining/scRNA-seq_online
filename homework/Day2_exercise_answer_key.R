#### Clustering

# After loading `seurat_integrated.RData`, check the object clusters with different resolution (0.4, 0.6, 0.8, 1.0, 1.4). For each resolution plot the corresponding UMAP and report how many clusters you observe. Which resolution do you think makes sense?


# Resolution 0.4: 13 clusters
# Resolution 0.6: 15 clusters
# Resolution 0.8: 17 clusters
# Resolution 1.0: 22 clusters
# Resolution 1.4: 27 clusters

#### Clustering quality control
# Answer:
# Hypothesize the clusters corresponding to each of the different clusters in the table:
# Cell Type	Clusters
# CD14+ monocytes	1, 3
# FCGR3A+ monocytes	10
# Conventional dendritic cells	14
# Plasmacytoid dendritic cells	16
# Marcrophages	-
# B cells	11, 7, 13
# T cells	0, 2, 6, 4, 5, 9
# CD4+ T cells	4, 0, 6, 2
# CD8+ T cells	5, 9
# NK cells	8, 12
# Megakaryocytes	15
# Erythrocytes -

#### Marker identification
# In the previous lesson, we identified cluster 9 as FCGR3A+ monocytes by inspecting the expression of known cell markers FCGR3A and MS4A7. Use FindConservedMarkers() function to find conserved markers for cluster 10. What do you observe? Do you see FCGR3A and MS4A7 as highly expressed genes in cluster 10?
# Answer: Yes, FCGR3A and MS4A7 are among the top highly expressed genes when performing FindConservedMarkers() for cluster 9. This observation confirms cluster 10 as FCGR3A+ monocytes.
