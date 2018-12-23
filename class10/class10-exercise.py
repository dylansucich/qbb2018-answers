#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import datasets
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from sklearn.cluster import KMeans
from scipy.stats import ttest_ind

df = pd.read_table(sys.argv[1], header=0, index_col=0)
# scipy.cluster.hierarchy.linkage(, method='single', metric='euclidean', optimal_ordering=False)

cell_types = df.columns.values.tolist()
gene_names = df.index.tolist()

matrix = df.values


# gene_names = matrix[:, 0]
# print(gene_names)
# print(cell_types)
# Transpose axis of matrix

Z = linkage(matrix, "ward")
ZT = linkage(df.T, "ward")

fig, ax = plt.subplots()
plt.title("Dendrogram")
plt.xlabel("sample")
plt.ylabel("distance")
dendrogram(
    ZT,
    show_leaf_counts=False,
    leaf_rotation=90.,
    leaf_font_size=12.,
    show_contracted=True)

fig.savefig("dendrogram_" + sys.argv[2] + ".png")
plt.close(fig)

idx_rows = leaves_list(Z)
data = matrix[idx_rows, :]
idx_columns = leaves_list(ZT)
data = data[:, idx_columns]

X = (data-np.average(data,axis=0))/np.std(data,axis=0)

m = np.max(np.abs(X))
fig, ax = plt.subplots(figsize=(8, 6)) 

ax.set_title("Heatmap of Genes by Cell Type") # Add a title to the top
im = ax.pcolor(                              # Treat the values like pixel intensities in a picture
	X,                                       # ... Using X as the values
	cmap="RdBu",                             # ... Use the Red-white-blue colormap to assign colors to your pixel values
	vmin=-1*m,                               # ... Set the lowest value to show on the scale
	vmax=m,                                  # ... Set the highest value to show on the scale. Since we are using a 'diverging' colormap, these should match.
	)
    
ax.grid(False) 
ax.set_xticks(                      # Edit the xticks being shown
	np.arange(0.5, X.shape[1]+0.5), # ... use the values centered on each column of pixels
	)
ax.set_xticklabels(                 # Label the ticks
	cell_types,                         # ... at position which correspond to the indices of our labels
	rotation=50,                    # ... and rotate the labels 50 degrees counter-clockwise
	)
# ax.set_yticks(idx_rows, rotation=50)                   # Edit the ticks on the y-axis to show....NOTHING

ax.set_yticklabels(                 # Label the ticks
	idx_rows,                         # ... at position which correspond to the indices of our labels
	rotation=50,                    # ... and rotate the labels 50 degrees counter-clockwise
	)

cbar = fig.colorbar(im, ax=ax)      # Add a bar to the right side of the plot which shows the scale correlating the colors to the pixel values

fig.subplots_adjust( # Adjust the spacing of the subplots, to help make everything fit
    left = 0.05,     # ... the left edge of the left-most plot will be this percent of the way across the width of the plot
    bottom = 0.15,   # ... the bottom edge of the bottom-most plot will be this percent of the way up the canvas
    right = 1.0,     # ... the right edge of the right-most plot will be this percent of the way across the width
    top = 0.95,      # ... the top edge of the top-most plot will be this percent of the way from the bottom
)

fig.savefig("2xClustered_" + sys.argv[2] + "_heatmap.png") # Save the image
plt.close(fig) # Close the canvas


cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)
ax = sns.clustermap(X , metric="euclidean", standard_scale=1, method="ward", cmap="Reds")
# plt.show()
ax.savefig("cluster_seaborn.png")
plt.close()



kmeans = KMeans(n_clusters = 4)
kmeans.fit(X)
y_means = kmeans.predict(X)

fig, ax = plt.subplots()
plt.scatter(X[:, 0], X[:, 1], c = y_means, s = 5, cmap = "viridis")

plt.ylabel("poly")
plt.xlabel("CFU")

fig.savefig("kmeans.png")
plt.close(fig)

early = ["CFU", "mys"]
late = ["poly", "unk"]

t_stat, p_val = ttest_ind(df[early], df[late], axis = 1)
df["p_value"] = p_val

df = df.mask(df["p_value"] > 0.05).dropna(how = "any").sort_values ("p_value")

#list of differentially expressed genes
diffExp = df.ix[:,4]
df = pd.DataFrame(diffExp)
df.to_csv("diffExp.csv", sep="\t")




# fpkms = matrix[:, 1:]

# for row in matrix:
     # pass
     # gene_names = genes.append.matrix[row,0]