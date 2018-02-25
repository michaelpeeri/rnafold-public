# Test recovering "typical profile" from bimodal distribution, e.g., given two classes of profiles, find the mode of the distribution
# (Since profiles are continuous, use kde)
#
# Ref: https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns; sns.set()
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import KFold, LeaveOneOut
from mfe_plots import loadProfileData, getProfileHeatmapTile, getLegendHeatmapTile


# Configuration
#K = 30  # Profile length
#N = 500 # Number of profiles
fileglob = "gcdata_v2_taxid_*_profile_1000_10_begin_0.h5"
profileLengthSamples = 31
profileStep = 10
typicalProfileOutputCSV = "find_typical_profile.out.profile.csv"
fixedYrange = (-2.9, 2.9)
#bwEstimationPoints = 1000
bwEstimationPoints = 150

def getFileNames():
    from glob import glob
    return glob(fileglob)

(_, _, _, _, _, _, _, biasProfiles, _, _) = loadProfileData(getFileNames())

print( "Loaded %d profiles" % len(biasProfiles) )

#import sys
#sys.exit(0)

x = np.vstack( [x[:profileLengthSamples] for x in biasProfiles.values()] )
print(x.shape)

N, K = x.shape

sns.distplot( x[:,0] )
plt.title("Distribution (pos=0nt)")
plt.savefig("find_typical_profile.out.1.pdf")
plt.close()

sns.distplot( x[:,2] )
plt.title("Distribution (pos=20nt)")
plt.savefig("find_typical_profile.out.2.pdf")
plt.close()

sns.distplot( x[:,4] )
plt.title("Distribution (pos=40nt)")
plt.savefig("find_typical_profile.out.3.pdf")
plt.close()

sns.distplot( x[:,6] )
plt.title("Distribution (pos=60nt)")
plt.savefig("find_typical_profile.out.4.pdf")
plt.close()



# referenceProfile = biasProfiles[511145][:profileLengthSamples]

# Determine the optimal bandwidth, using cross-validation
optimalBW = []
for k in range(0,K,1):
    bandwidths = 10 ** np.linspace(-2.0, 1, bwEstimationPoints)
    x1 = np.expand_dims(x[:,k], 1)
    
    cv = KFold(len(x1), n_folds=10)
    #cv = LeaveOneOut(len(x1))
    
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                        {'bandwidth': bandwidths},
                        cv=cv )
    grid.fit(x1)
    print("%d: %s" % (k, grid.best_params_))
    optimalBW.append( grid.best_params_['bandwidth'] )

#import sys
#sys.exit(0)

x_d = np.linspace( -2, 2, 5000 )
#print(x_d.shape)
estimatedValues = np.zeros(K)


for k in reversed(range(0,K,1)):  # go over profiles (for plotting ease, do it in reversed)
    #bandwidth = 0.05
    bandwidth = optimalBW[k]
    print("%d: bw=%g" % (k, bandwidth))
    kde = KernelDensity(bandwidth=bandwidth, kernel='gaussian')
    kde.fit(np.expand_dims(x[:,k], 1))

    logprob = kde.score_samples( np.expand_dims(x_d, 1) )

    peakEstimate = x_d[np.argmax(logprob)]
    estimatedValues[k] = peakEstimate
    
    if( True ): #k%1==0 ):
        plt.fill_between(x_d, np.exp(logprob), alpha=0.2)
        #plt.plot(x, np.full_like(x, -0.01), '|k', markeredgewidth=1)
        #plt.ylim(-0.02, 0.25)

#print( (estimatedValues - actualValues) ** 2 )


plt.grid(True)
plt.savefig("find_typical_profile.out.5.pdf")
plt.close()


u, v = x.shape
plotdata = np.moveaxis( np.stack( (np.broadcast_to([range(v)], (u,v)), x ) ), 0, 2)
assert(plotdata.shape == (u,v,2))

fig, ax1 = plt.subplots()
# Plot density using KDE
ax = sns.kdeplot( plotdata[:,:,0].flatten(), plotdata[:,:,1].flatten(), n_levels=50, bw=0.65, cmap="Blues", shade=True, shade_lowest=True, legend=True )
# Plot raw data (jittered for clarity)
l = plt.plot( plotdata[:,:,0].flatten() + np.random.randn(u*v)*0.1, plotdata[:,:,1].flatten(), ".", alpha=1.0 )
plt.setp( l, "markersize", 3)

# Plot neutral line
plt.plot([0,v], [0,0], '-k')

# Plot window width
plt.plot([0.5, 0.5 + 40/profileStep], [-1.9,-1.9], '-r', linewidth=4)
plt.annotate( s="window", xy=(1.0, -1.83) )

# Annotate number of items
plt.annotate( s='n = %d' % u, xy=(26.0, 1.8) )
plt.xlim([-0.5, v-0.5])
plt.ylim([-2, 2])

# Plot 'typical profile'
plt.plot( np.linspace(0, v, len(estimatedValues)), estimatedValues, '-y', linewidth=1.0 )

# Convert tick labels from array positions to CDS positions (in steps determined by profileStep, e.g., 10)
def format_fn(tick_val, tick_pos):
    return str(int(tick_val*profileStep))
ax.xaxis.set_major_formatter(FuncFormatter(format_fn))


plt.savefig("find_typical_profile.out.6.pdf")
plt.close()





typicalProfilePd = pd.DataFrame( { "CDS_position": np.linspace(0, v, len(estimatedValues)),  "Profile": estimatedValues, "Optimal_KDE_BW": optimalBW } )
typicalProfilePd.to_csv( typicalProfileOutputCSV )


print(getProfileHeatmapTile( 99990010, pd.DataFrame( { 99990010 : estimatedValues } ), fixedYrange, True, profileStep ))

print(getLegendHeatmapTile( fixedYrange ))





