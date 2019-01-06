import numpy as np
from sklearn.cluster import KMeans, DBSCAN

# get vectors somehow
inputVectors = ....

# perform K-means clustering:
kmeans = KMeans( n_clusters=K )
model1 = kmeans.fit( inputVectors )

# extract output:
labels = model1.labels_
centers = model1.cluster_centers_

# perform DBSCAN clustering:
dbscan = DBSCAN( metric="correlation" )
model2 = dbscan.fit( inputVectors )

# extract output:
labels = model2.labels_
