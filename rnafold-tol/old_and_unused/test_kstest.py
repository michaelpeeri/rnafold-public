#from math import log10
#from bz2 import BZ2File
#from random import randint
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
#plt.style.use('ggplot') # Use the ggplot style

targetDistribution = stats.gamma(a=4.0, scale=2.0)

N = 10000
K = 200
c = 0

def calcPvalue( statistic, distribution ):
    dist = stats.gamma( a=distribution[0], loc=0.0, scale=distribution[2] )
    statistics = np.array( [ stats.kstest( dist.rvs( size=50 ), dist.cdf ).statistic for i in range(K) ] )
    global c
    print(c); c += 1
    return float(len(statistics[statistics > statistic])) / K

# Construct samples of size 50
trainingData = [ targetDistribution.rvs(size=50) for i in range(N) ]

ksTestStatistics = [ stats.kstest(trainingData[i], targetDistribution.cdf).statistic for i in range(N) ]

# Fit a distribution for each sample
fittedDistributions = [ stats.gamma.fit( trainingData[i], floc=0.0 ) for i in range(N) ]

# 
EmpiricalPvalues = [ calcPvalue( ksTestStatistics[i], fittedDistributions[i] ) for i in range(1000) ]
#print(Pvalues)

# Plot the distribution of test statistics for 50-samples against the correct distribution
fig, ax = plt.subplots()
quants = np.arange(0.0, 1.0, 0.002)

data = pd.DataFrame(ksTestStatistics)
qdata = data.quantile(q = quants)

empscores = pd.DataFrame( EmpiricalPvalues )
qempscores = empscores.quantile(q = quants)

plt.plot(qdata, label='K-S test statistic')

plt.plot( qempscores, label='Empirical scores')

ax.legend(fontsize='small', loc='upper left')

plt.grid()
plt.title('title')
plt.xlabel('Cumulative Probability')
plt.ylabel('K-S test statistic')

plt.savefig('ksStatisticSimulation.pdf')




