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



#taxid = 559292
taxid = int(sys.argv[1])
fn1 = "wilcoxon_taxid_%d.csv" % taxid
fn2 = sys.argv[2]

allpoints = []

# NP_416900.1,5393,0.522,2540,0.0511,-8.05,0.475,-0.4
dfWilcoxon = pd.read_csv(fn1, header=None, names=('ProtId','Statistic1', 'SelectPval1', 'Statistic2', 'SelectPval2', 'NativeMFE', 'GC', 'Select'), na_values='None')

dfPA = pd.read_csv(fn2, header=None, names=('ProtId', 'PApercentile'))

#s = np.log10(dfWilcoxon['Pval1'].values)

selection_logPval = np.log10(dfWilcoxon['SelectPval1'].values)
selection_logPval_signed = np.abs( np.log10(dfWilcoxon['SelectPval1'].values) ) * np.sign( dfWilcoxon['Select'] )


nativeMFE = dfWilcoxon['NativeMFE']
gc = dfWilcoxon['GC']

x = pd.merge(dfWilcoxon, dfPA, how='left')
#print(x)
print(len(gc))
#print(x.iloc[20:40])
#print(x.iloc[40:60])
#print(x.iloc[60:80])
#print(x.iloc[80:100])
nopaidx = x['PApercentile'].isnull()
paidx = x['PApercentile'].notnull()
print(len(x[nopaidx]))
print(len(x[paidx]))

print(x[paidx]['PApercentile'].describe())

#papercent = x[paidx]['PApercentile']

#print(dfWilcoxon['GC'].describe())


#--------------------------------------------------------------------------------


fig, ax = plt.subplots()


xvals = gc
yvals = dfWilcoxon['Select']

# Linear correlation and factors
pearson = pearsonr(xvals, yvals)
spearman = spearmanr(xvals, yvals)
kendall = kendalltau(xvals, yvals)
l = linregress(xvals, yvals)

abline_x = np.arange(min(xvals), max(xvals), 0.1)
abline_y = abline_x * l.slope + l.intercept
plt.plot(abline_x, abline_y)

# plot the linear approximation
plt.annotate(s="Pearson r: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                 xy=(0.2, -2.50),  fontsize=6 )
plt.annotate(s="Pearson r^2: %1.3f"  % (pearson[0]**2,),                              xy=(0.2, -2.75),  fontsize=6 )
plt.annotate(s="Spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue),  xy=(0.2, -3.00),  fontsize=6 )
plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(0.2, -3.25),  fontsize=6 )
plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(0.2, -3.50),  fontsize=6 )

# class MyNorm(matplotlib.colors.Normalize):
#     def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#         matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

#     def __call__(self, value, clip=None):
#         print(value)
#         for i in range(len(value)):
#             if value[i] is None:
#                 value[i] = 0.0
#             else:
#                 value[i] = value[i]*0.5 + 0.5
#         return value

print(x.describe())
#print(len(x['PApercentile'<0.50]))

print(x['PApercentile'].values)
plt.plot(np.arange(0.25, 0.75, 0.1), [0.0] * 5)
#xcmp = plt.cm.plasma
#norm1 = MyNorm(vmin=0.0, vmax=1.0)
#a = x['PApercentile'].values.tolist()
#plt.scatter(gc, dfWilcoxon['Select'], s=7, c=x['PApercentile'].values, cmap=xcmp, norm=norm1, vmin=0.0, vmax=1.0)
#plt.scatter(gc, dfWilcoxon['Select'], s=7, c=a, cmap=xcmp, norm=norm1, vmin=0.0, vmax=1.0)


plt.scatter(gc[nopaidx], dfWilcoxon[nopaidx]['Select'], s=7, color="black")
plt.scatter(gc[paidx], dfWilcoxon[paidx]['Select'], s=8, c=x[paidx]['PApercentile'], cmap=plt.cm.copper, vmin=0.0, vmax=1.0)

plt.xlabel("GC% (mean)")
plt.ylabel("Selection (delta MFE)")
#plt.xlim([-2,0])
#plt.ylim([-2,0])
plt.grid(True)
#plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
plt.savefig("genelevel_gc_selection_taxid_%d.pdf" % taxid)
plt.savefig("genelevel_gc_selection_taxid_%d.svg" % taxid)
plt.close(fig)

#--------------------------------------------------------------------------------
fig, ax = plt.subplots()

xvals = gc
yvals = nativeMFE

# Linear correlation and factors
pearson = pearsonr(xvals, yvals)
spearman = spearmanr(xvals, yvals)
kendall = kendalltau(xvals, yvals)
l = linregress(xvals, yvals)

abline_x = np.arange(min(xvals), max(xvals), 0.1)
abline_y = abline_x * l.slope + l.intercept
plt.plot(abline_x, abline_y)

# plot the linear approximation
plt.annotate(s="Pearson r: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                 xy=(0.65, -0.5),  fontsize=6 )
plt.annotate(s="Pearson r^2: %1.3f"  % (pearson[0]**2,),                              xy=(0.65, -1.0),  fontsize=6 )
plt.annotate(s="Spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue),  xy=(0.65, -1.5),  fontsize=6 )
plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(0.65, -2.0),  fontsize=6 )
plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(0.65, -2.5),  fontsize=6 )


plt.scatter(gc, nativeMFE, s=10)
plt.xlabel("GC% (mean)")
plt.ylabel("Native MFE")
#plt.xlim([-2,0])
#plt.ylim([-2,0])
plt.grid(True)
#plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
plt.savefig("genelevel_gc_nativeMFE_taxid_%d.pdf" % taxid)
plt.savefig("genelevel_gc_nativeMFE_taxid_%d.svg" % taxid)
plt.close(fig)

#--------------------------------------------------------------------------------


fig, ax = plt.subplots()


xvals = nativeMFE
yvals = dfWilcoxon['Select']

# Linear correlation and factors
pearson = pearsonr(xvals, yvals)
spearman = spearmanr(xvals, yvals)
kendall = kendalltau(xvals, yvals)
l = linregress(xvals, yvals)

abline_x = np.arange(min(xvals), max(xvals), 0.1)
abline_y = abline_x * l.slope + l.intercept
plt.plot(abline_x, abline_y)

# plot the linear approximation
plt.annotate(s="Pearson r: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                 xy=(-2.0, -0.5),  fontsize=6 )
plt.annotate(s="Pearson r^2: %1.3f"  % (pearson[0]**2,),                              xy=(-2.0, -1.0),  fontsize=6 )
plt.annotate(s="Spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue),  xy=(-2.0, -1.5),  fontsize=6 )
plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(-2.0, -2.0),  fontsize=6 )
plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(-2.0, -2.5),  fontsize=6 )


plt.plot(abline_x, abline_x*0.0)
plt.scatter(nativeMFE, dfWilcoxon['Select'], s=10)
plt.xlabel("Native MFE")
plt.ylabel("Selection (delta MFE)")
#plt.xlim([-2,0])
#plt.ylim([-2,0])
plt.grid(True)
#plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
plt.savefig("genelevel_nativeMFE_selection_taxid_%d.pdf" % taxid)
plt.savefig("genelevel_nativeMFE_selection_taxid_%d.svg" % taxid)
plt.close(fig)

#--------------------------------------------------------------------------------


fig, ax = plt.subplots()

#fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', gridspec_kw={'height_ratios': [4, 1], 'width_ratios': [4,1]})

#fig.add_subplot(1,1,1)

plt.plot([min(gc),max(gc)], [0.0, 0.0], color="black")

plt.setp(ax.get_xticklabels(), fontsize=16)
plt.setp(ax.get_yticklabels(), fontsize=16)
#plt.scatter(gc, selection_logPval_signed, s=10)

plt.scatter(gc[nopaidx], selection_logPval_signed[nopaidx], s=11, color="gray")
#plt.scatter(gc[paidx], selection_logPval_signed[paidx], s=13, c=x[paidx]['PApercentile'], cmap=plt.cm.summer, vmin=0.0, vmax=1.0)
plt.scatter(gc, selection_logPval_signed, s=15, color=(0.1,0.2,0.5), alpha=0.4)


#plt.annotate(s="mean=%g"%(np.mean(selection_logPval_signed[selection_logPval_signed<0.0])), xy=(0.5, -5))
#plt.annotate(s="mean=%g"%(np.mean(selection_logPval_signed[selection_logPval_signed>0.0])), xy=(0.5, -10))
plt.annotate(s="n= %d" % (len(selection_logPval_signed)), xy=(0.40, 20), fontsize=20)

plt.xlabel("GC% (mean)", fontsize=16)
plt.ylabel("Signed selection (log(P-val) * sign)", fontsize=16)
#plt.colorbar()
plt.xlim([min(gc), max(gc)])
plt.ylim([min(selection_logPval_signed), max(selection_logPval_signed)])
plt.grid(True)
#plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')


#fig.add_subplot(1,1,2)
#hist, histbins = np.histogram( selection_logPval_signed, bins = 50, normed=True )
#ax3.bar(histbins[:-1], hist)


#fig.add_subplot(1,1,3)
#hist, histbins = np.histogram( gc, bins = 50, normed=True )
#ax2.bar(histbins[:-1], hist)

plt.savefig("genelevel_gc_signed_selection_taxid_%d.pdf" % taxid)
plt.savefig("genelevel_gc_signed_selection_taxid_%d.svg" % taxid)
plt.close(fig)
