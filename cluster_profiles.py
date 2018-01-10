from sklearn.cluster import KMeans
from sklearn import metrics
import numpy as np


def analyzeProfileClusters(profilesArray, n_init=10000, max_rms_distance_to_merge_clusters=0.04 ):
    if( profilesArray.shape[0] < 2 ):
        raise Exception("Can't cluster %d profiles" % profilesArray.shape[0])

    #range_n_clusters = (1, 2, 3)
    range_n_clusters = (1, 2)

    results = {}
    
    for n_clusters in range_n_clusters:

        if n_clusters >= profilesArray.shape[0]:  # Number of clusters must be >= number of points
            continue
        
        kmeans = KMeans(n_clusters=n_clusters, n_init = n_init )
        model = kmeans.fit(profilesArray)
        labels = model.labels_
        centers = model.cluster_centers_

        dist = 0
        if len(centers) == 2:
            dist = metrics.mean_squared_error(centers[0], centers[1])
        assert(len(centers) <= 2 )

        # Sort the centeroids so they appear in a consistent order; update the labels accordingly
        counts = [sum([1 for x in labels if x==i]) for i in range(len(centers))]  # count how many profiles belong to each cluster
        centers_order = np.argsort( -np.array(counts) ) # sort in descending order
        assert(centers_order.shape == (n_clusters,))
        centers = centers[centers_order,:]  # sort the centeroids
        labels = [centers_order[x] for x in labels]  # update the labels to match the sorted centroids
        
        
        results[n_clusters] = (centers, labels, dist)
        #print("%d\t%g\t%g" % (n_clusters,
        #                      metrics.silhouette_score(profilesArray, labels ),
        #                      metrics.calinski_harabaz_score( profilesArray, labels) ) )

    if( 2 in results and results[2][2] >= max_rms_distance_to_merge_clusters ):
        return results[2]
    else:
        return results[1]
        
        
        
