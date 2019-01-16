import numpy as np
from sklearn import metrics


class CompleteLinkageClustering(object):
    def __init__(self, threshold, metric=metrics.pairwise.paired_euclidean_distances):
        self._threshold = threshold
        self._metric = metric


    def findPeers(self, cluster):
        return [idx for idx,clusterId in enumerate(self._assignments) if clusterId==cluster]
        
    def fit(self, data):
        N = data.shape[0]

        self._assignments = list(range(N))
        assignments = self._assignments

        allGroups = set(assignments)
        iter = 1

        while( True ):
            print(iter)

            groupsForMerging = []
            allGroups = set(assignments)

            for n in allGroups:
                for m in allGroups - frozenset((n,)):
                    assert(n != m)
                    if m < n:
                        continue

                    peers_n = self.findPeers(n)
                    peers_m = self.findPeers(m)

                    pairs = np.array([[x,y] for x in peers_n for y in peers_m])
                    distances = self._metric( data[pairs[:,0],:], data[pairs[:,1],:] )
                    if max(distances) < self._threshold:
                        groupsForMerging.append( (n,m) )

            if not groupsForMerging:
                break
            
            alreadyMerged = set()
            for x,y in groupsForMerging:
                if x in alreadyMerged:
                    continue
                if y in alreadyMerged:
                    continue

                for i in range(N):
                    if assignments[i] == y:
                        assignments[i] = x

            iter += 1
            
        # Debug only
        self.labels_ = assignments
        return self
