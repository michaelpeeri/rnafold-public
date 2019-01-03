# Cluster profiles to find a small number of "representative profiles" (this is only meant for plotting, not downstream analysis)
# We'd like to make sure each returned centroid is similar to all members of its cluster (and is not a simply the average of a number of sub-groups).
# Outliers are not "noise" and should be represented by a cluster of their own.
# Note that the clustering approach is used for convenience - we don't claim the clusters represent actual distinct groups (e.g., all profiles might be part of a continuous range; that range should be broken up to "clusters" to illustrate the variation in the data).
# One challenge when applying clustering algorithms in this setting is that N might be as small as 1 (so e.g., density-based methods are not applicable).
import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn import metrics
import scipy.spatial.distance
from sklearn import decomposition
#from scipy.stats import spearmanr
from collections import Counter
#from mfe_plots import saveHistogram

# Configuration
#metric = metrics.mean_squared_error
defaultMetric = metrics.pairwise.paired_euclidean_distances
correlationMetric = scipy.spatial.distance.correlation

"""
For a given K-means clustering solution, find the maximum distance between any point to the centroid of its assigned cluster.
[N -- number of points]
[K -- number of clusters]
points   - raw data points (before clustering) [N-vector]
clusters - assigned cluster for each point [N-vector of indices in the ranges 0..K-1]
centers  - centroid for each cluster  [K-vector of vectors]
"""
def findMaxDistanceToCentroid( points, clusters, centers, metric=defaultMetric ):
    assert( len(clusters) == points.shape[0] )    # 
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

        messageAlreadyPrinted = False

        for otherCentroid in centers:  # sanity test only
            x2 = metric( np.array(thisPoint).reshape(1,-1), np.array(otherCentroid).reshape(1,-1) )[0]
            #print("({})".format(x2))
            if( x2 < distToCentroid ):
                #print("Warning: Distance to other centroid ({}) < dist to our centroid ({}) for point {} ({}), otherCentroid {}. Points: {} clusters{} centers {} ".format(x2, distToCentroid, i, thisPoint, otherCentroid, points, clusters, centers ))
                #print("Warning: Distance to other centroid ({}) < dist to our centroid ({}) for point {} ({}), otherCentroid {}. clusters{}".format(x2, distToCentroid, i, thisPoint, otherCentroid, clusters ))
                if not messageAlreadyPrinted:
                    print("Warning: Distance to other centroid ({}) < dist to our centroid ({}) for point {}".format(x2, distToCentroid, i ))
                    messageAlreadyPrinted = True
                    
                maxDistanceToCentroid = 1e9  # Make sure this solution isn't selected
                
            #assert( metric( np.array(thisPoint).reshape(1,-1), np.array(otherCentroid).reshape(1,-1) )[0] >= distToCentroid ) # TODO RESTORE THIS !!!
            
            
    return maxDistanceToCentroid


def sortVectorsUsingPCA(vectors, labels=None):

    #print(type(vectors))
    
    #dictForLabelsTesting = {}
    #if not labels is None:
    #    dictForLabelsTesting = dict(zip(map(tuple, vectors.tolist()), labels))

    X = vectors # .copy() #np.vstack(vectors) # convert dict of vectors to matrix

    pca = decomposition.PCA()
    pca.fit(X)

    pca.n_components = 1  # force 1 components
    X_reduced = pca.fit_transform(X)
    #print(X_reduced.shape)   # must be (K,1)

    K = X_reduced.shape[0]
    decorated = zip(X_reduced.tolist(), X.tolist(), range(K))
    sorted_ = sorted(decorated, key=lambda x: x[0])
    print(sorted_)
    centroids = [x[1] for x in sorted_]
    if not labels is None:
        labelsMapping = dict( zip( [x[2] for x in sorted_], range(K) ) )
        labels = [labelsMapping[x] for x in labels]

    return np.asarray(centroids), labels


def calcDiversityMetrics( points, metric=defaultMetric, saveHistogramAs=None ):
    pairwiseDistances = []
    absDeviations = []
    
    points = points[~np.any(np.isnan(points), axis=1)]  # remove points containig NaN
    if points.shape[0] < 2:
        return 0.0

    center = np.mean( points, axis=0 ).reshape(1,-1)
    
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            assert(i<j)
            pairwiseDistances.append( metric( np.asarray(points[i]).reshape(1,-1),
                                              np.asarray(points[j]).reshape(1,-1) ) )

        absDeviations.append( metric( center,
                                      np.asarray(points[i]).reshape(1,-1) ) )

    #if (not saveHistogramAs is None) and len(pairwiseDistances)>1:
    #    saveHistogram(np.asarray(pairwiseDistances), saveHistogramAs)
        
    #if (not saveHistogramAs is None) and len(absDeviations)>1:
    #    saveHistogram(np.asarray(absDeviations), "absDev_" + saveHistogramAs)

    return np.median( np.asarray(absDeviations) )


"""
Note - 'metric' is used only for the criterion for increasing the number of clusters; For the clustering itself, Euclidean distance is always used (since this is K-means...)
"""
def analyzeProfileClustersUsingKMeans(profilesArray, n_init=100, max_permissible_distance_centroid=0.2, max_clusters=15, metric=defaultMetric ):
    if( profilesArray.shape[0] < 2 ):
        raise Exception("Can't cluster {} profiles".format( profilesArray.shape[0]) )

    range_n_clusters = range(1,max_clusters+1)

    results = None

    maxDistanceToCentroid = 0.0

    profilesArray = profilesArray[~np.any(np.isnan(profilesArray), axis=1)]  # remove points containig NaN

    for K in range_n_clusters:

        if K > profilesArray.shape[0]:  # Number of clusters must be >= number of points
            break
        
        attempt = 1
        maxDistanceToCentroid = 1e9
        
        while( maxDistanceToCentroid > max_permissible_distance_centroid and attempt <= 15 ):
        
            print("K={} N={} (attempt={})".format(K, profilesArray.shape[0], attempt ))

            kmeans = KMeans( n_clusters=K,
                             n_init = n_init if K>1 else 100,
                             max_iter=10000*attempt,
                             tol=1e-2,
                             verbose=False )
            
            model = kmeans.fit(profilesArray)
            labels = model.labels_
            centers = model.cluster_centers_

            # ----------------------------------------------------------------------------------------------------
            # Sort the centeroids so they appear in a consistent order; update the labels accordingly
            # ----------------------------------------------------------------------------------------------------
            # TODO - fix bug causing this to disrupt the correct order in some cases!
            # ----------------------------------------------------------------------------------------------------
            #counts = [sum([1 for x in labels if x==i]) for i in range(len(centers))]  # count how many profiles belong to each cluster
            #centers_order = np.argsort( -np.array(counts) ) # sort in descending order
            #assert(centers_order.shape == (K,))
            #centers = centers[centers_order,:]  # sort the centeroids
            #labels = [centers_order[x] for x in labels]  # update the labels to match the sorted centroids
            # ----------------------------------------------------------------------------------------------------


            #print("%d\t%g\t%g" % (K,
            #                      metrics.silhouette_score(profilesArray, labels ),
            #                      metrics.calinski_harabaz_score( profilesArray, labels) ) )

            maxDistanceToCentroid = findMaxDistanceToCentroid( profilesArray, labels, centers, metric=metric )
            
            if K==1 or maxDistanceToCentroid < 1e8: # this will be the last clustering attempt
                break   # for K==1, no use retrying (we should try larger K immediately)
            
            else:
                print("...... {}".format( maxDistanceToCentroid ) )
                attempt += 1 # We may repeat the clustering in the hope of finding a better solution

        #print("centers={}".format(centers))
        #print("labels={}".format(labels))
        # Clustering using the current K is done; did we find an acceptable solution?
        if maxDistanceToCentroid <= max_permissible_distance_centroid:
            before = centers.shape
            beforec = sorted(Counter(labels).values())
            #print("Before: {}".format(labels))
            #print(beforec)
            centers, labels = sortVectorsUsingPCA(centers, labels)
            #print("After:  {}".format(labels))
            #print( sorted(Counter(labels).values()) )
            assert( centers.shape == before )  # labels should appear in the new order, but no label should be reassigned
            assert( sorted(Counter(labels).values()) == beforec )
            
            return (centers, labels, maxDistanceToCentroid, max_permissible_distance_centroid)
        else:
            print("... Criterion not met: {} > {}".format( maxDistanceToCentroid, max_permissible_distance_centroid ) )

    raise Exception("No clustering solution met criterion of maxDistance < %.3g (actual maxDistance: %.3g)" % (max_permissible_distance_centroid, maxDistanceToCentroid) )

# def calcCorrelationDistances(vectors):
#     N = vectors.shape[0]
#     out = np.zeros((N,N))

#     for i in range(N):
#         for j in range(i+1,N):
#             assert(j>i)
#             corr = spearmanr(vectors[i,:],vectors[j,:])
#             out[i,j] = 1 - corr[0]
#             out[j,i] = 1 - corr[0]
            
#     return out
            

def analyzeProfileClustersUsingAggClus(profilesArray, distThreshold=0.5):
    profilesArray = profilesArray[~np.any(np.isnan(profilesArray), axis=1)]  # remove points containing NaN

    print("Calculating distances...")
    #distances = calcCorrelationDistances(profilesArray)
    distances = metrics.pairwise_distances(profilesArray, metric='correlation')

    N = profilesArray.shape[0]  # Number of clusters must be >= number of points
    
    for nClusters in range(2, max(15, N)):
    
        #aggclus = AgglomerativeClustering( n_clusters=max( 2, min( 8, N ) ),
        aggclus = AgglomerativeClustering( n_clusters=nClusters,
                                           linkage='complete',
                                           affinity='precomputed' )

        model = aggclus.fit(distances)
        labels = model.labels_

        count = 0
        centers = []

        print("ALL) {} >= corr >= {} (mu={})".format(np.min(distances), np.max(distances), np.mean(distances) ))

        maxInnerDistance = 0.0

        for currCluster in range(max(labels)+1):
            ixs = [i for i,x in enumerate(labels) if x==currCluster]
            clusterMembers = profilesArray.take(ixs, axis=0)
            count += clusterMembers.shape[0]
            centers.append(np.mean(clusterMembers, axis=0) )

            innerDistances = metrics.pairwise_distances(clusterMembers, metric='correlation')
            print("{}) {} >= corr >= {} (mu={})".format(currCluster, np.min(innerDistances), np.max(innerDistances), np.mean(innerDistances) ))
            maxd = np.max(innerDistances)
            maxInnerDistance = max(maxd, maxInnerDistance)

        assert(count==N)
        
        if maxInnerDistance <= distThreshold:
            break

    centers, labels = sortVectorsUsingPCA(np.asarray(centers), labels)
        
    return (centers, labels, maxd, distThreshold)
    
def plotDistancesDistribution(profilesArray, savePlotAs=None):
    
    profilesArray = profilesArray[~np.any(np.isnan(profilesArray), axis=1)]  # remove points containig NaN

    print("Plotting all distance for {} profiles...".format(profilesArray.shape[0]))
    
    distances = metrics.pairwise_distances(profilesArray, metric='correlation')

    data = []
    for i in range(profilesArray.shape[0]):
        for j in range(i, profilesArray.shape[0]):
            data.append(distances[i,j])

    #if not savePlotAs is None:
    #    saveHistogram(np.asarray(data), savePlotAs)



def analyzeProfileClusters(profilesArray, method="KMeans", *args, **kw):
    if method=="KMeans":
        return analyzeProfileClustersUsingKMeans(profilesArray, *args, **kw)
    elif method=="AggClus":
        return analyzeProfileClustersUsingAggClus(profilesArray, *args, **kw)
    else:
        raise Exception("Unknown method: '{}'".format(method))
