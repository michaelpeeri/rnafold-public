import sys
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
#plt.style.use('ggplot') # Use the ggplot style
import pandas as pd
import numpy as np
from scipy import stats

#taxid = 223926
taxid = int(sys.argv[1])
fShapiroPval = 'shapiropval_taxid_%d.csv' % taxid
fSkew        = 'skew_taxid_%d.csv' % taxid
fKurtosis    = 'kurtosis_taxid_%d.csv' % taxid

dfShapiroPval = pd.read_csv(fShapiroPval, header=None, names=('ProtId',) + tuple(range(150)), na_values='None')
dfSkew        = pd.read_csv(fSkew,        header=None, names=('ProtId',) + tuple(range(150)), na_values='None')
dfKurtosis    = pd.read_csv(fKurtosis,    header=None, names=('ProtId',) + tuple(range(150)), na_values='None')

ncount = 100
df_ref_norm       = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))
df_ref_logistic   = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))
df_ref_gamma      = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))
df_ref_chi        = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))

df_ref_skew_norm  = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))
df_ref_skew_gamma = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))
df_ref_skew_chi   = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))

df_ref_kurt_norm  = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))
df_ref_kurt_gamma = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))
df_ref_kurt_chi   = pd.DataFrame(columns=('ProtId',)+tuple(range(ncount)), index=range(101), dtype=np.dtype('float32'))

mygamma = stats.gamma(a=13.5)

print("Creating simulated values...")
for i in range(df_ref_norm.shape[0]):
    for j in range(df_ref_norm.shape[1]):
        # Note - the parameters shouldn't matter for normal and logistic distributions
        df_ref_norm.iloc[i,j]     = stats.shapiro( stats.norm.rvs(size=50) )[1]
        df_ref_logistic.iloc[i,j] = stats.shapiro( stats.logistic.rvs(size=50) )[1]
        df_ref_gamma.iloc[i,j]    = stats.shapiro( -1.0 * mygamma.rvs(size=50) )[1]
        df_ref_chi.iloc[i,j]      = stats.shapiro( -1.0 * stats.chi.rvs(3, size=50) )[1]
        
        df_ref_skew_norm.iloc[i,j]     = stats.skewtest(        stats.norm.rvs(size=50) ).statistic
        df_ref_skew_gamma.iloc[i,j]    = stats.skewtest( -1.0 * mygamma.rvs(size=50) ).statistic
        df_ref_skew_chi.iloc[i,j]      = stats.skewtest( -1.0 * stats.chi.rvs(3, size=50) ).statistic

        df_ref_kurt_norm.iloc[i,j]     = stats.kurtosis(        stats.norm.rvs(size=50) )
        df_ref_kurt_gamma.iloc[i,j]    = stats.kurtosis( -1.0 * mygamma.rvs(size=50) )
        df_ref_kurt_chi.iloc[i,j]      = stats.kurtosis( -1.0 * stats.chi.rvs(3, size=50) )
print("Done!")

fig, ax = plt.subplots()

#cax = plt.imshow(df, interpolation='none', cmap=matplotlib.cm.coolwarm)
#cbar = fig.colorbar(cax, ticks=[0.0, 0.05, 1.0], orientation='horizontal')
#plt.savefig('shapiro.png', dpi=900)
#plt.close()

quants = np.arange(0.0, 1.0, 0.01)

x1 = dfShapiroPval.iloc[:,1:].stack()
q1 = x1.quantile(q = quants)
plt.plot(q1, label='MFE distribution (taxid=%d)' % taxid)

x2 = df_ref_norm.iloc[:,1:].stack()
q2 = x2.quantile(q = quants)
plt.plot(q2, label='Normal')

x3 = df_ref_logistic.iloc[:,1:].stack()
q3 = x3.quantile(q = quants)
plt.plot(q3, label='Logistic')

x3b = df_ref_gamma.iloc[:,1:].stack()
q3b = x3b.quantile(q = quants)
plt.plot(q3b, label='-1*Gamma(a=13.5)')

x3c = df_ref_chi.iloc[:,1:].stack()
q3c = x3c.quantile(q = quants)
plt.plot(q3c, label='-1*Chi(df=3)')

ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('Shapiro-Wilk normality test (N=50 samples)')
plt.xlabel('Cumulative Probability')
plt.ylabel('Shapiro-Wilk p-value')

plt.savefig('shapiro_hist_taxid_%d.pdf' % taxid)

fig2, ax = plt.subplots()
quants = np.arange(0.0, 1.0, 0.001)

x4 = dfSkew.iloc[:,1:].stack()
q4 = x4.quantile(q = quants)
plt.plot(q4, label='MFE distribution (taxid=%d)' % taxid)

x5 = df_ref_skew_norm.iloc[:,1:].stack()
q5 = x5.quantile(q = quants)
plt.plot(q5, label='Normal')

x6 = df_ref_skew_gamma.iloc[:,1:].stack()
q6 = x6.quantile(q = quants)
plt.plot(q6, label='-1*Gamma(a=13.5)')

x7 = df_ref_skew_chi.iloc[:,1:].stack()
q7 = x7.quantile(q = quants)
plt.plot(q7, label='-1*Chi(df=3)')


ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('Distribtion skewness (N=50 samples)')
plt.xlabel('Cumulative Probability')
plt.ylabel('Skewness')

plt.savefig('skew_hist_taxid_%d.pdf' % taxid)
#-----------------------------------------------------------------

fig3, ax = plt.subplots()
quants = np.arange(0.0, 1.0, 0.001)

x10 = dfKurtosis.iloc[:,1:].stack()
q10 = x10.quantile(q = quants)
plt.plot(q10, label='MFE distribution (taxid=%d)' % taxid)

x11 = df_ref_kurt_norm.iloc[:,1:].stack()
q11 = x11.quantile(q = quants)
plt.plot(q11, label='Normal')

x11 = df_ref_kurt_gamma.iloc[:,1:].stack()
q11 = x11.quantile(q = quants)
plt.plot(q11, label='-1*Gamma(a=13.5)')

x11 = df_ref_kurt_chi.iloc[:,1:].stack()
q11 = x11.quantile(q = quants)
plt.plot(q11, label='-1*Chi(df=3)')


ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('Distribtion kurtosis (N=50 samples)')
plt.xlabel('Cumulative Probability')
plt.ylabel('Kurtosis')

plt.savefig('kurt_hist_taxid_%d.pdf' % taxid)


#-----------------------------------------------------------------

plt.close()
