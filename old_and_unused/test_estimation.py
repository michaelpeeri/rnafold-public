import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
#plt.style.use('ggplot') # Use the ggplot style

targetDistribution = stats.gamma(a=4.0, scale=2.0)

# Number of simulations
N = 10000

# Sample size
S = 50000


simulatedSamples = [ targetDistribution.rvs(size=S) for i in range(N) ]  # -> len(simulatedSamples[n<N]) == S

# Fit a distribution for each sample
estimatedDistributions = [ stats.gamma.fit( simulatedSamples[i], floc=0.0 ) for i in range(N) ]

percentiles = np.arange(0.01, 1.0, 0.01)
percentileValues = targetDistribution.ppf( percentiles ) # ppf is the inverse of the cdf
# dist.cdf( dist.ppf( x )) == x (for 0>x>1)

resultingPercentiles = np.array( [ stats.gamma(a=dist[0], loc=0.0, scale=dist[2]).cdf( percentileValues )  for dist in estimatedDistributions ] )

#print(len(resultingPercentiles))
#print(len(resultingPercentiles[0]))
#print(resultingPercentiles[0])

#xsorted = np.sort(resultingPercentiles, axis=1 )
print('----- 0.05 ------')
x05 = np.sort(resultingPercentiles[:,4])
#print(x05)
print(x05[int(0.05*N)])
print(np.mean(x05))
print(x05[int(0.95*N)])

print('----- 0.50 ------')
x50 = np.sort(resultingPercentiles[:,49])
print(x50[int(0.05*N)])
print(np.mean(x50))
print(x50[int(0.95*N)])

print('----- 0.95 ------')
x95 = np.sort(resultingPercentiles[:,94])
print(x95[int(0.05*N)])
print(np.mean(x95))
print(x95[int(0.95*N)])
print('.')
print(x95)

#print(xsorted)
#

fig, ax = plt.subplots()
quants = np.arange(0.0, 1.0, 0.001)

data = pd.DataFrame(resultingPercentiles, columns=map(lambda x:'%g'%x, percentiles))
#data = pd.DataFrame(resultingPercentiles, columns=percentiles)
#print(xsorted.shape)
#print(len(xsorted[:,4]))
#print(len(quants))
#data['0.50'] = np.quantile
#data['0.05'] = xsorted[:, 4]
#data['0.50'] = xsorted[:,49]
#data['0.95'] = xsorted[:,94]

#print(data)
#print(data.columns)
#print(data['0.95'])

q05 = data['0.05'].quantile(q=quants)
q25 = data['0.25'].quantile(q=quants)
q50 = data['0.5'].quantile(q=quants)
q75 = data['0.75'].quantile(q=quants)
q95 = data['0.95'].quantile(q=quants)
#print(d50)

data2 = pd.DataFrame(index=quants, columns=(0.05, 0.25, 0.50, 0.75, 0.95))
data2[0.05] = q05
data2[0.25] = q25
data2[0.50] = q50
data2[0.75] = q75
data2[0.95] = q95
#d50b = data.quantile(0.50, axis=1)
#print(d50b)



plt.plot(data2[0.95], label='0.95')
plt.plot(data2[0.75], label='0.75')
plt.plot(data2[0.50], label='0.50')
plt.plot(data2[0.25], label='0.25')
plt.plot(data2[0.05], label='0.05')
plt.grid()
plt.title('S = %d' % S)
plt.xlabel('Cumulated probability')
plt.ylabel('ylabel')
plt.legend()

plt.savefig('estimated_minmax.pdf')






