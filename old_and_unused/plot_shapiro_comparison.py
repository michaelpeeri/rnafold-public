import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
#plt.style.use('ggplot') # Use the ggplot style
import pandas as pd
import numpy as np
from scipy import stats

taxons = (223926, 511145, 284811, 1125630, 93061, 44056)

fig, ax = plt.subplots()
quants = np.arange(0.0, 1.0, 0.001)

for taxid in taxons:
    fShapiroPval = 'shapiropval_taxid_%d.csv' % taxid
    dfShapiroPval = pd.read_csv(fShapiroPval, header=None, names=('ProtId',) + tuple(range(150)), na_values='None')


    x1 = dfShapiroPval.iloc[:,1:].stack()
    q1 = x1.quantile(q = quants)
    plt.plot(q1, label='taxid=%d' % taxid)

    
ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('Shapiro-Wilk normality test (N=50 samples)')
plt.xlabel('Cumulative Probability')
plt.ylabel('Shapiro-Wilk p-value')

plt.savefig('shapiro_comparison.pdf')

#----------------------------------------------------------
fig2, ax = plt.subplots()

for taxid in taxons:
    fSkew        = 'skew_taxid_%d.csv' % taxid
    dfSkew        = pd.read_csv(fSkew,        header=None, names=('ProtId',) + tuple(range(150)), na_values='None')

    x4 = dfSkew.iloc[:,1:].stack()
    q4 = x4.quantile(q = quants)
    plt.plot(q4, label='taxid=%d' % taxid)



ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('Skew (N=50 samples)')
plt.xlabel('Cumulative Probability')
plt.ylabel('Skew')

plt.savefig('skew_comparison.pdf')

#----------------------------------------------------------
fig3, ax = plt.subplots()

for taxid in taxons:
    fKurtosis    = 'kurtosis_taxid_%d.csv' % taxid
    dfKurtosis    = pd.read_csv(fKurtosis,    header=None, names=('ProtId',) + tuple(range(150)), na_values='None')

    x10 = dfKurtosis.iloc[:,1:].stack()
    q10 = x10.quantile(q = quants)
    plt.plot(q10, label='taxid=%d' % taxid)

ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('Kurtosis (N=50 samples)')
plt.xlabel('Cumulative Probability')
plt.ylabel('Kurtosis (Normal=0)')

plt.savefig('kurtosis_comparison.pdf')


plt.close()
