# Cluster profiles to find a small number of "representative profiles" (this is only meant for plotting, not downstream analysis)
# We'd like to make sure each returned centroid is similar to all members of its cluster (and is not a simply the average of a number of sub-groups).
# Outliers are not "noise" and should be represented by a cluster of their own.
# Note that the clustering approach is used for convenience - we don't claim the clusters represent actual distinct groups (e.g., all profiles might be part of a continuous range; that range should be broken up to "clusters" to illustrate the variation in the data).
# One challenge when applying clustering algorithms in this setting is that N might be as small as 1 (so e.g., density-based methods are not applicable).
from sklearn.cluster import KMeans
from sklearn import metrics
import numpy as np

# Configuration
#metric = metrics.mean_squared_error
defaultMetric = metrics.pairwise.paired_euclidean_distances

"""
For a given K-means clustering solution, find the maximum distance between any point to the centroid of its assigned cluster.
"""
def findMaxDistanceToCentroid( points, clusters, centers, metric=defaultMetric ):
    assert( len(clusters) == points.shape[0] ) 
    assert( len(centers) == len(set(clusters)) )

    maxDistanceToCentroid = 0.0
    
    for i in range(len(clusters)):  # iterate over each point
        clusterOfThisPoint = clusters[i]
        centroidOfThisPoint = centers[clusterOfThisPoint]
        thisPoint = points[i]
        distToCentroid = metric( np.array(thisPoint).reshape(1,-1), np.array(centroidOfThisPoint).reshape(1,-1) )[0]
        #print(">> {}".format(distToCentroid))
        if distToCentroid > maxDistanceToCentroid:
            maxDistanceToCentroid = distToCentroid

        for otherCentroid in centers:  # sanity test only
            x2 = metric( np.array(thisPoint).reshape(1,-1), np.array(otherCentroid).reshape(1,-1) )[0]
            #print("({})".format(x2))
            if( x2 < distToCentroid ):
                #print("Warning: Distance to other centroid ({}) < dist to our centroid ({}) for point {} ({}), otherCentroid {}. Points: {} clusters{} centers {} ".format(x2, distToCentroid, i, thisPoint, otherCentroid, points, clusters, centers ))
                print("Warning: Distance to other centroid ({}) < dist to our centroid ({}) for point {} ({}), otherCentroid {}. clusters{}".format(x2, distToCentroid, i, thisPoint, otherCentroid, clusters ))
                maxDistanceToCentroid = 1e9  # Make sure this solution isn't selected
                
            #assert( metric( np.array(thisPoint).reshape(1,-1), np.array(otherCentroid).reshape(1,-1) )[0] >= distToCentroid ) # TODO RESTORE THIS !!!
            
            
    return maxDistanceToCentroid
    

def analyzeProfileClusters(profilesArray, n_init=10000, max_permissible_distance_centroid=0.2, max_clusters = 15, metric=defaultMetric ):
    if( profilesArray.shape[0] < 2 ):
        raise Exception("Can't cluster {} profiles".format( profilesArray.shape[0]) )

    range_n_clusters = range(1,max_clusters+1)

    results = None

    maxDistanceToCentroid = 0.0
    
    for K in range_n_clusters:

        if K > profilesArray.shape[0]:  # Number of clusters must be >= number of points
            break
        
        attempt = 1
        maxDistanceToCentroid = 1e9
        
        while( maxDistanceToCentroid > max_permissible_distance_centroid and attempt <= 10 ):
        
            print("K={} N={} (attempt={})".format(K, profilesArray.shape[0], attempt ))

            kmeans = KMeans(n_clusters=K, n_init = n_init, max_iter=1000, tol=1e-5, verbose=0 )
            model = kmeans.fit(profilesArray)
            labels = model.labels_
            centers = model.cluster_centers_

            # Sort the centeroids so they appear in a consistent order; update the labels accordingly
            counts = [sum([1 for x in labels if x==i]) for i in range(len(centers))]  # count how many profiles belong to each cluster
            centers_order = np.argsort( -np.array(counts) ) # sort in descending order
            assert(centers_order.shape == (K,))
            centers = centers[centers_order,:]  # sort the centeroids
            labels = [centers_order[x] for x in labels]  # update the labels to match the sorted centroids


            #print("%d\t%g\t%g" % (K,
            #                      metrics.silhouette_score(profilesArray, labels ),
            #                      metrics.calinski_harabaz_score( profilesArray, labels) ) )

            maxDistanceToCentroid = findMaxDistanceToCentroid( profilesArray, labels, centers, metric=metric )
            
            if K==1:
                break   # for K==1, no use retrying (we should try larger K immediately)
            else:
                print("...... {}".format( maxDistanceToCentroid ) )
                attempt += 1
        
        
        if maxDistanceToCentroid <= max_permissible_distance_centroid:
            return (centers, labels, maxDistanceToCentroid)
        else:
            print("... Criterion not met: {} > {}".format( maxDistanceToCentroid, max_permissible_distance_centroid ) )

    raise Exception("No clustering solution met criterion of maxDistance < %.3g (actual maxDistance: %.3g)" % (max_permissible_distance_centroid, maxDistanceToCentroid) )

