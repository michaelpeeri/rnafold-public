import csv
import sys
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style



taxid = 559292
#taxid = 511145
fn1 = "wilcoxon_taxid_%d.csv" % taxid

allpoints = []

# "%s,%.4g,%.3g\n"
dfWilcoxon = pd.read_csv(fn1, header=None, names=('ProtId','Statistic1', 'Pval1', 'Statistic2', 'Pval2'), na_values='None')

allpointsx = np.log10(dfWilcoxon['Pval1'].values)

allpointsy = np.log10(dfWilcoxon['Pval2'].values)




fig, ax = plt.subplots()
# Plot the diagonal (equality) line
plt.plot(range(-15,1), range(-15,1))
plt.plot(range(-15,1), [np.log10(0.05)] * 16)
plt.plot([np.log10(0.05)] * 16, range(-15,1))

# Plot the actual data
plt.scatter(allpointsx, allpointsy, s=10)

plt.xlabel("method 1")
plt.ylabel("method 2")

#plt.xlim([-2,0])
#plt.ylim([-2,0])

plt.grid(True)
#plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
plt.savefig("wilcoxon_method_comparison_taxid_%d.pdf" % taxid)
plt.savefig("wilcoxon_method_comparison_taxid_%d.svg" % taxid)
plt.close(fig)

