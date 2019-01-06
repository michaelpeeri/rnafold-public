import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style

x = np.arange(0.3, 0.7, 0.005)
y = x*0.08 - 0.04 + norm.rvs(0.0, 0.015, size=len(x))
print(x)
print(y)
#print(norm.rvs(size=10))

plt.scatter(x, y)


pearson = pearsonr(x,y)
spearman = spearmanr(x,y)
kendall = kendalltau(x,y)
#print(spearmanr(x,y))
l = linregress(x,y)

abline_x = np.arange(0.3, 0.8, 0.1)
abline_y = abline_x * l.slope + l.intercept
plt.plot(abline_x, abline_y)

plt.annotate(s="pearson r: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                xy=(0.6, 0.36) )
plt.annotate(s="pearson r^2: %1.3f"  % (pearson[0]**2,),                xy=(0.6, 0.042) )
plt.annotate(s="spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue), xy=(0.6, 0.048) )
plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(0.6, 0.054) )



plt.grid(True)
plt.savefig("test_pearson.pdf")
plt.savefig("test_pearson.svg")
plt.close()

