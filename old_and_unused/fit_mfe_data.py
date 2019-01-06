from math import log10
from bz2 import BZ2File
from random import randint
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
#plt.style.use('ggplot') # Use the ggplot style


filename = 'rawdata_taxid_511145.csv.bz2'
useNpartsOf = 150   # Fraction of lines to be included. 602 = Number of lines in rawdata_taxid_511145.csv.bz2 (601799) / 1000 (so that we get 50*1000 values for fitting the distribution)
c = 0
u = 0
v = 0

def linesSource(filename):
    with BZ2File(filename, 'r') as f:
        for line in f:
            yield line

def convertLineFields(fields):
    protId = fields[0]
    position = int(fields[1])
    MFEscores = np.array( map(lambda x: None if x=='None' else float(x), fields[2:]), dtype=np.dtype('float32'), ndmin=2) * -1
    MFEscores[MFEscores < 1e-9] = 1e-9
    return (protId, position, MFEscores)

        
trainData = pd.DataFrame(columns=range(50))
testData = pd.DataFrame(columns=range(50))
gammaScores = []
gammaScores2 = []
normScores = []

for line in linesSource(filename):
    # Line format: ProtId,CDSPos,Shuffle0,....,Shuffle49
    data = line.split(',')
    assert(len(data)==52)
    c += 1
    
    # select a sample of 1/useNPartsOf lines to be included
    r = randint( 0, useNpartsOf )
    if( r == 0 ):
        u += 1

        (_, _, MFEscores) = convertLineFields(data)

        newdata = pd.DataFrame(MFEscores, columns=range(50))
        trainData = trainData.append(newdata)

        gammaParams = stats.gamma.fit( MFEscores, floc=0.0 )
        if( not( gammaParams[0] > 0.02 and gammaParams[0] < 100 )):
            print("Warning: gamma fit failed (%s)" % str(gammaParams))
            continue
        if( not( gammaParams[2] > 0.1 and gammaParams[2] < 100 )):
            print("Warning: gamma fit failed (%s)" % str(gammaParams))
            continue

        normParams = stats.norm.fit( MFEscores )
        if( not( normParams[1] > 0.01 and normParams[1] < 100 )):
            print("Warning: normal fit failed (%s)" % str(normParams))
            continue

        gammaQual = stats.kstest( newdata, stats.gamma( a=gammaParams[0], scale=gammaParams[2] ).cdf ).statistic
        gammaScores.append( gammaQual)

        # TEST ONLY #####  TEST ONLY #####  TEST ONLY #####  TEST ONLY #####  TEST ONLY ##### 
        gammaQual2 = stats.kstest( newdata, stats.gamma( a=gammaParams[0]*1.2, scale=gammaParams[2] ).cdf ).statistic
        gammaScores2.append( gammaQual2)
        # TEST ONLY #####  TEST ONLY #####  TEST ONLY #####  TEST ONLY #####  TEST ONLY ##### 
        
        #print(gammaQual, log10(gammaQual) )
        
        #print(normParams)
        normQual = stats.kstest( newdata, stats.norm( loc=normParams[0], scale=normParams[1] ).cdf ).statistic
        normScores.append( normQual)


    elif( r == 1 ):
        v += 1

        (_, _, MFEscores) = convertLineFields(data)

        newdata = pd.DataFrame(MFEscores)
        testData = testData.append(newdata)

# Test for equal variances
assert(len(trainData.get_values()[0]) == 50)
print(stats.fligner( *trainData.get_values() ))

fig, ax = plt.subplots()
quants = np.arange(0.0, 1.0, 0.002)

dfgs = pd.DataFrame(gammaScores)
qdfgs = dfgs.quantile(q = quants)

dfgs2 = pd.DataFrame(gammaScores2)
qdfgs2 = dfgs2.quantile(q = quants)


dfns = pd.DataFrame(normScores) 
qdfns = dfns.quantile(q = quants)

#dflo = dfgs / dfns
#qdflo = dflo.quantile(q = quants)

plt.plot(qdfgs, label='gamma')
plt.plot(qdfgs2, label='gamma2')
plt.plot(qdfns, label='norm')
#plt.plot(qdflo, label='lo')

ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('title')
plt.xlabel('Cumulative Probability')
plt.ylabel('log P')

plt.savefig('qdfgs.pdf')



print(trainData.shape)
print(len(trainData))
print(testData.shape)
print(len(testData))

print("Used %d out of %d lines as training data" % (u, c))
print("Used %d out of %d lines as test data"     % (v, c))

td = trainData.iloc[:,1:].stack()
print(td.shape)
print(td.describe())
fitParams = stats.gamma.fit(td, floc=0.0)
fitted = stats.gamma( a= fitParams[0], scale=fitParams[2] )
print(fitParams)

print("Scoring fitted distribution against training-data:")
print(stats.kstest(td, fitted.cdf ) )

print("Scoring fitted distribution against test-data:")
xd = testData.iloc[:1,:].stack()
print(stats.kstest(xd, fitted.cdf ) )

print("Scoring fitted distribution against data sampled from the fitted distribution:")
print(stats.kstest(fitted.rvs(size=xd.shape[0]), fitted.cdf ) )
