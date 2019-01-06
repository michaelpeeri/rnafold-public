# Test recovering "typical profile" from bimodal distribution, e.g., given two classes of profiles, find the mode of the distribution
# (Since profiles are continuous, use kde)
#
# Ref: https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import KFold, LeaveOneOut

K = 30  # Profile length
N = 500 # Number of profiles

actualValues = np.linspace( -10, 10, K )

def make_data(N=500, K=150, f=0.3):
    #x = np.random.randn(N, K) + actualValues
    x = np.random.randn(N, K) + actualValues
    y = np.random.randn(N, K)*2.0 + 8 - actualValues*0.7
    choice = np.expand_dims( np.random.binomial( 1, 1-f, N ), 1 )
    #print(choice)
    x = x * choice
    y = y * (1-choice)
    
    x = x + y
    return x

x = make_data(N,K,0.3)

print(x.shape)
print(x[:,None].shape)
#hist = plt.hist(x.flatten(), bins=30, normed=True)

estimatedValues = np.zeros(K)

#for k in range(0,K,1):
#    bandwidths = 10 ** np.linspace(-1, 1, 100)
#    x1 = np.expand_dims(x[:,k], 1)
#    
#    cv = KFold(len(x1), n_folds=10)
#    #cv = LeaveOneOut(len(x1))
#    
#    grid = GridSearchCV(KernelDensity(kernel='gaussian'),
#                        {'bandwidth': bandwidths},
#                        cv=cv )
#    grid.fit(x1)
#    print(grid.best_params_)
#
#import sys
#sys.exit(0)

x_d = np.linspace( -20, 20, 5000 )
#print(x_d.shape)

for k in range(0,K,1):
    kde = KernelDensity(bandwidth=0.37, kernel='gaussian')
    kde.fit(np.expand_dims(x[:,k], 1))

    logprob = kde.score_samples( np.expand_dims(x_d, 1) )

    peakEstimate = x_d[np.argmax(logprob)]
    estimatedValues[k] = peakEstimate
    
    if( True ): #k%1==0 ):
        plt.fill_between(x_d, np.exp(logprob), alpha=0.2)
        #plt.plot(x, np.full_like(x, -0.01), '|k', markeredgewidth=1)
        #plt.ylim(-0.02, 0.25)

print( (estimatedValues - actualValues) ** 2 )


plt.grid(True)
plt.savefig("test_kde.out.pdf")
plt.close()


