import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
#plt.style.use('ggplot') # Use the ggplot style

taxid = 511145
fn1 = "wilcoxon_taxid_%d.csv" % taxid

allpoints = []

# "%s,%.4g,%.3g\n"
dfWilcoxon = pd.read_csv(fn1, header=None, names=('ProtId','Statistic1', 'Pval1', 'Statistic2', 'Pval2'), na_values='None')


fig, ax = plt.subplots()

quants = np.arange(0.0, 1.0, 0.01)


q1 = dfWilcoxon['Pval1'].quantile(q = quants)
plt.plot(q1, label='Wilcoxon method 1, P-val distribution (taxid=%d)' % taxid)

q2 = dfWilcoxon['Pval2'].quantile(q = quants)
plt.plot(q2, label='Wilcoxon method 2, P-val distribution (taxid=%d)' % taxid)

plt.plot(q1, np.ones(quants.size)*0.05, label='0.05')


ax.legend(fontsize='small', loc='upper left')
plt.grid()
plt.title('Wilcoxon signed-rank test')
plt.xlabel('Cumulative Probability')
plt.ylabel('Wilcoxon p-value')

plt.savefig('wilcoxon_pval_taxid_%d.pdf' % taxid)
