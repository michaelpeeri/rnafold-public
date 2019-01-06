import sys
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress, wilcoxon
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style


def readChipsFile(filename):
    return pd.read_csv(filename, sep='\s+', index_col=0, usecols=(0,3), header=0, names=('gene', 'ENc'), na_values='None')

df1 = readChipsFile(sys.argv[1])
df1.rename(columns={'ENc':'ENc_1'}, inplace=True)
df2 = readChipsFile(sys.argv[2])
df2.rename(columns={'ENc':'ENc_2'}, inplace=True)

df3 = pd.concat([df1, df2], axis=1)
df3.dropna(axis=0, how='any', inplace=True)
print(df3)

fig, ax = plt.subplots()




# Linear correlation and factors
pearson = pearsonr(df3['ENc_1'], df3['ENc_2'])
spearman = spearmanr(df3['ENc_1'], df3['ENc_2'])
kendall = kendalltau(df3['ENc_1'], df3['ENc_2'])
l = linregress(df3['ENc_1'], df3['ENc_2'])

abline_x = np.arange(20, 65, 1)
abline_y = abline_x * l.slope + l.intercept
plt.plot(abline_x, abline_y, '--', zorder=10, lw=0.3, alpha=0.5)

plt.plot([20,65], [20,65], c='black', zorder=9, lw=0.3, alpha=0.5)


df3.plot(x='ENc_1', y='ENc_2', kind='scatter', alpha=0.32, s=6, zorder=11, ax=ax)

sorted = (df3['ENc_1']-df3['ENc_2']).sort_values()
print(sorted.head(10))
print(sorted.tail(10))
print(np.abs(df3['ENc_1']-df3['ENc_2']).sort_values().head(10))


topr = 30
scaler = 1.5
# plot the linear approximation
plt.annotate(s="y = %.3g*x + %.3g"  % (l.slope, l.intercept),                         xy=(54, topr-scaler*0),  fontsize=6 )
plt.annotate(s="Pearson r: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                 xy=(54, topr-scaler*1),  fontsize=6 )
plt.annotate(s="Pearson r^2: %1.3f"  % (pearson[0]**2,),                              xy=(54, topr-scaler*2),  fontsize=6 )
plt.annotate(s="Spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue),  xy=(54, topr-scaler*3),  fontsize=6 )
plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(54, topr-scaler*4),  fontsize=6 )
plt.annotate(s="n= %d"  % len(df3['ENc_2']),                                          xy=(54, topr-scaler*5),  fontsize=6 )

plt.xlim(20,65)
plt.ylim(20,65)


plt.savefig("ENc_comparison.pdf")
plt.savefig("ENc_comparison.svg")
plt.close(fig)


