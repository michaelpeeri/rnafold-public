import csv
import sys
from scipy import stats
import matplotlib
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style

fn1 = sys.argv[1]
fn2 = sys.argv[2]

allpointsx = []
allpointsy = []

with open(fn1, 'rb') as f1, open(fn2, 'rb') as f2:
    a = csv.reader(f1)
    b = csv.reader(f2)

    for row1, row2 in zip(a,b):
        rowid = row1[0]
        data1 = list(map(float, row1[1:]))
        data2 = list(map(float, row2[1:]))

        N = min(len(data1), len(data2))
        assert(N>=150)

        #print(len(data1[:N]), len(data2))
        print(stats.spearmanr(data1[:N], data2[:N]))
        #l1 = a.next()
        #l2 = b.next()

        allpointsx.extend(data1[:N])    
        allpointsy.extend(data2[:N])    

fig, ax = plt.subplots()
# Plot the diagonal (equality) line
plt.plot(range(-25,1), range(-25,1))
# Plot the actual data
plt.scatter(allpointsx, allpointsy, s=10)

plt.xlabel(fn1)
plt.ylabel(fn2)
plt.grid(True)
#plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
plt.savefig("mfe_impl_comparison.pdf")
plt.savefig("mfe_impl_comparison.svg")
plt.close(fig)

