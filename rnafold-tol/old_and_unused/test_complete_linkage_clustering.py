import numpy as np
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from complete_linkage_clustering import CompleteLinkageClustering

# modified from: https://stackoverflow.com/a/14314054
def moving_average(a, n=3, axis=None, dtype=None) :
    ret = np.cumsum(a, axis=axis, dtype=dtype)
    #print(ret)
    ret[:,n:] = ret[:,n:] - ret[:,:-n]
    return ret[:,n - 1:] / n


#print(moving_average(np.array([[1, 1, 1, 1, 7, 1, 4, 1, 1], [2, 2, 2, 2, 0, 2, 2, 0, 0]]), 3, axis=1, dtype=float))

def makeRandomProfiles(numProfiles = 100, profileLen = 31, mean = 0.0, sd=1.0, numPivots=4):
    breadth = profileLen/(numPivots-1)
                        
    pivotPositions = range(0, profileLen+1, breadth )
    assert( len(pivotPositions) == numPivots )
    
    pivots = np.random.rand( numProfiles * numPivots ).reshape( numProfiles, numPivots) * sd + mean
    assert( pivots.shape == (numProfiles, numPivots) )

    data = np.repeat(pivots, breadth, axis=1)
    assert( data.shape[0] == numProfiles )

    #print([1/float(breadth)]*breadth )
    data = moving_average(data, breadth, axis=1 )
    
    noise = np.random.rand( data.shape[0] * data.shape[1] ).reshape( data.shape[0], data.shape[1] ) * 0.1
    data = data + noise
    assert(data.shape[0] == numProfiles )
    # TODO - will data.shape[1] always equal profileLen???

    return data
    

data = makeRandomProfiles(numProfiles=600)

#clust = AgglomerativeClustreing( eps=0.48, min_samples=1 ).fit(data)
clust = CompleteLinkageClustering( threshold=1.15 )
result = clust.fit(data)
labels = result.labels_
print(result.labels_)

pca = PCA()
pca.fit(data)
print(pca.explained_variance_)
pca.n_components = 2
data_reduced = pca.fit_transform(data)
print(data_reduced.shape)
fig, ax = plt.subplots()
plt.scatter( data_reduced[:,0], data_reduced[:,1], c=labels )
plt.savefig("test_dbscan_clustering.out.pdf")
plt.close()
