# Try to learn the parameters of an unknown distribution, based on a large number of small samples.
# This is different than the more usual task of learning the parameters using a single large sample (which can be done using ML estimators).
# Note - Why not just combine the samples into a single large sample?
#        The idea is that the samples all come from the same family of distributions, but not necessarily from the same distribution.
#        e.g., some of the parameters may be different.

from bz2 import BZ2File
from random import randint
from scipy import stats, optimize
import numpy as np

hiddenDistribution = stats.gamma(a=13.5, scale=2.0)

def getSample():
    return hiddenDistribution.rvs(size=50)

# Data from which to learn
#dataSet = list(map(lambda x: getSample(), range(1000)))
filename = 'rawdata_taxid_511145.csv.bz2'

def linesSource(filename):
    with BZ2File(filename, 'r') as f:
        for line in f:
            yield line

def convertLineFields(fields):
    protId = fields[0]
    position = int(fields[1])
    MFEscores = np.array( map(lambda x: None if x=='None' else float(x), fields[2:]), dtype=np.dtype('float32'), ndmin=2) * -1
    if( np.any(np.isnan(MFEscores)) ):
        return (protId, position, None)
    
    MFEscores[MFEscores < 1e-6] = 1e-6
    return (protId, position, list(MFEscores))

dataSet = []
testSet = []
useNpartsOf = 600
#testData = pd.DataFrame(columns=range(50))
for line in linesSource(filename):
    # Line format: ProtId,CDSPos,Shuffle0,....,Shuffle49
    data = line.split(',')
    assert(len(data)==52)
    
    # select a sample of 1/useNPartsOf lines to be included
    r = randint( 0, useNpartsOf )
    if( r == 0 ):
        (_, _, MFEscores) = convertLineFields(data)
        if( MFEscores is not None ):
            dataSet.append(MFEscores)
        

    if( r == 1 ):
        (_, _, MFEscores) = convertLineFields(data)
        if( MFEscores is not None ):
            testSet.append(MFEscores)

def calcDistanceToDistribution(samples, cdf):
    return stats.kstest(samples, cdf).pvalue

def calcCombinedDistance(cdf):
    a = map(lambda sample: calcDistanceToDistribution(sample, cdf), dataSet)
    #return stats.gmean(a)
    return 1.0 - np.mean(a)

d0 = stats.gamma(a=3.0, loc=0.0, scale=2.0)
print(calcCombinedDistance(d0.cdf))

def objective(x):
    return calcCombinedDistance( stats.gamma( a=x[0], scale=x[1] ).cdf )

#print(test( [12.0, 3.0] ))

#d1b = stats.gamma(a=12.2, scale=2.8)
#print(calcCombinedDistance(d1b.cdf))

#d2 = stats.gamma(a=13.49, scale=1.95)
#print(calcCombinedDistance(d2.cdf))

#d2b = stats.gamma(a=13.5, scale=2.0)
#print(calcCombinedDistance(d2b.cdf))

#d3 = stats.norm(loc=10.0, scale=2.0)
#print(calcCombinedDistance(d3.cdf))

x0 = [3.0, 2.0]
res = optimize.minimize(objective, x0, method='nelder-mead', options={'xtol': 1e-4, 'disp': True})
print(res)

print("Test against training set: ")
print stats.kstest( np.array(dataSet).ravel(), stats.gamma( a=res[0], scale=res[1] ).cdf )

print("Test against test set: ")
print stats.kstest( np.array(testSet).ravel(), stats.gamma( a=res[0], scale=res[1] ).cdf )
