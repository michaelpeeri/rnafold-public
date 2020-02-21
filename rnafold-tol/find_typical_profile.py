# RNAFold - Analyze mRNA folding bias (Local Folding Energy) based on randomizations.
# Copyright (C) 2016-2020 Michael Peeri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# Test recovering "typical profile" from bimodal distribution, e.g., given two classes of profiles, find the mode of the distribution
# (Since profiles are continuous, use kde)
#
# Ref: https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html
import matplotlib # ; matplotlib.use('Agg')#; matplotlib.rc('text', usetex=True)
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns; sns.set()
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import KFold, LeaveOneOut
from mfe_plots import loadProfileData, getProfileHeatmapTile, getLegendHeatmapTile, getHeatmaplotProfilesValuesRange



#----------------------------------------------------------------------------------------------------------
# Configuration
#K = 30  # Profile length
#N = 500 # Number of profiles
#fileglob = "gcdata_v2_taxid_*_profile_1000_10_begin_0.h5"
fileglob = "gcdata_v2_taxid_*_profile_310_10_begin_0_t11.h5"
#fileglob = "gcdata_v2_taxid_*_profile_310_10_end_0_t11.h5"
profileLengthSamples = 33
profileStep = 10
typicalProfileOutputCSV = "find_typical_profile.out.profile.csv"

# Externally-supplied Y range (for plotting the modal profile as a tile). Must match the range used for other profiles (if all of them are used together)
fixedYrange = (-2.8436242765763069, 2.8436242765763069)
#fixedYrange = None

#bwEstimationPoints = 1000
bwEstimationPoints = 150     # Number of bandwidth values to test when estimating the optimal KDE bandwidth
#----------------------------------------------------------------------------------------------------------

def getFileNames():
    from glob import glob
    return glob(fileglob)

(_1, _2, _3, _4, _5, _6, _7, biasProfiles, _9, _10, _11, _12) = loadProfileData(getFileNames())

print( "Loaded %d profiles" % len(biasProfiles) )

if fixedYrange is None:
    fixedYrange = getHeatmaplotProfilesValuesRange( biasProfiles )
print("Using profile scale: {}".format(fixedYrange))

#x = np.vstack( [x[:profileLengthSamples] for x in biasProfiles.values()] )
x = np.vstack( [x for x in biasProfiles.values()] )
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
    #bandwidths = 10 ** np.linspace(-3.0, 2, bwEstimationPoints)
    #x1 = np.expand_dims(x[:,k], 1)
    x1 = x[:,k]
    x1 = x1[~np.isnan(x1)]
    x1 = np.expand_dims(x1, 1)
    
    cv = KFold(len(x1), n_folds=10)
    #cv = LeaveOneOut(len(x1))
    
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                        {'bandwidth': bandwidths},
                        cv=cv )
    grid.fit(x1)
    print("%d: %s" % (k, grid.best_params_))
    optimalBW.append( grid.best_params_['bandwidth'] )


x_d = np.linspace( -2, 3, 5000 )  # This must cover the range of values.. TODO - fix this...
                                  # Note: precision = 5/5000 = 0.001
#print(x_d.shape)
estimatedValues = np.zeros(K)

# 0.0, 0.25, 0.5, 0.75, 1.0
quantileValues = np.zeros((5, K), dtype=np.float32)


def findQuantiles(x_d, logDensity):
    density = np.exp(logDensity)
    total = sum(density)

    quantiles = [0.25, 0.50, 0.75]
    nextPoint = 0
    sigma = 0.0

    inner = []

    for xi, di in zip(x_d, density):
        if sigma/total > quantiles[nextPoint]:
            inner.append( xi )
            nextPoint += 1
            if nextPoint>=len(quantiles):
                break
        sigma += di
    IQR = inner[2]-inner[0]
    IQRpoint = [ inner[0] - 1.5*IQR, inner[2] + 1.5*IQR ]
    print("IQR: {}".format(IQRpoint))

    

    

    print(IQR)
        
    return (inner, IQRpoint)
    


# Calculate the KDE for each point along the profile
for k in reversed(range(0,K,1)):  # go over profiles (for plotting ease, do it in reversed)

    bandwidth = optimalBW[k]  # get the optimal bandwidth (calculated above)
    print("%d: bw=%g" % (k, bandwidth))
    kde = KernelDensity(bandwidth=bandwidth, kernel='gaussian')
    xi = x[:,k]
    xi = xi[~np.isnan(xi)]
    print(xi.shape)
    xi = np.expand_dims(xi, 1)
    print(xi.shape)

    kde.fit(xi)  # calculate the KDE

    logprob = kde.score_samples( np.expand_dims(x_d, 1) )  # score the KDE over a range of values

    print( kde.score_samples( np.expand_dims([-2,-1,0,1, 2], 1) ) )  # test only....

    peakEstimate = x_d[np.argmax(logprob)]
    estimatedValues[k] = peakEstimate
    print( "peak: {}".format(peakEstimate))

    (qs, IQRs) = findQuantiles(x_d, logprob)
    IQRs[0] = max( IQRs[0], min(x[:,k]) )
    IQRs[1] = min( IQRs[1], max(x[:,k]) )
    assert(IQRs[0]     <= qs[0])
    assert(min(x[:,k]) <= qs[0])
    assert(qs[0]       <= qs[1])
    assert(qs[1]       <= qs[2])
    assert(qs[2]       <= max(x[:,k]) )
    assert(qs[2]       <= IQRs[1] )
    quantileValues[:,k] = [IQRs[0], qs[0], qs[1], qs[2], IQRs[1]]
    
    if( True ): #k%1==0 ):
        plt.fill_between(x_d, np.exp(logprob), alpha=0.2)
        #plt.plot(x, np.full_like(x, -0.01), '|k', markeredgewidth=1)
        #plt.ylim(-0.02, 0.25)

#print( (estimatedValues - actualValues) ** 2 )


plt.grid(True)
plt.savefig("find_typical_profile.out.5.pdf")
plt.close()

u, v = x.shape   # u = number of species;  v = profile length
plotdata = np.moveaxis( np.stack( (np.broadcast_to([range(v)], (u,v)), x ) ), 0, 2)
assert(plotdata.shape == (u,v,2))

fig, ax1 = plt.subplots(figsize=(8, 2.2))
# Plot density using KDE
ax = sns.kdeplot( plotdata[:,:,0].flatten(), plotdata[:,:,1].flatten(), n_levels=10, bw=0.65, cmap="Blues", shade=True, shade_lowest=True, legend=True )
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


fig, ax1 = plt.subplots(figsize=(8, 2.2))
# Plot density using KDE

#sns.kdeplot( plotdata[:,:,0].flatten(), plotdata[:,:,1].flatten(), n_levels=10, bw=0.65, cmap="Blues", shade=True, shade_lowest=True, legend=True, ax=ax1 )

#sns.lvplot(data=x, ax=ax1, color="#50d080")
#sns.boxplot(data=x, color="#50d080", notch=True, ax=ax1, linewidth=0.3, boxprops={'zorder': 10})
plt.boxplot(x)


# Plot neutral line
plt.plot([-1,v], [0,0], '-k')

# Plot window width
plt.plot([0.5, 0.5 + 40/profileStep], [-1.9,-1.9], '-r', linewidth=4)
plt.annotate( s="window", xy=(1.0, -1.83) )

# Annotate number of items
plt.annotate( s='n = %d' % u, xy=(26.0, 1.8) )
plt.xlim([-0.5, v-0.5])
plt.ylim([-2, 2.5])


plt.title("lvplot")
plt.savefig("find_typical_profile.out.7.pdf")
plt.close()



def findValueOutliers(biasProfiles, region, N=10):
    #geneLevelScatter = pd.DataFrame({'gc':pd.Series(dtype='float'), 'logpval':pd.Series(dtype='float'), 'abslogpval':pd.Series(dtype='float'), 'protid':pd.Series(dtype='string')})

    #profiles = 
    
    #profiles = geneLevelScatter = geneLevelScatter.append(pd.DataFrame({'gc':pd.Series([meanGC]), 'logpval': pd.Series([directedLogPval]), 'abslogpval': pd.Series([pvalue]), 'protid':pd.Series([protId]), 'pa':pd.Series([paval]), 'cds_length_nt':pd.Series([cds_length_nt])}))


    keys = biasProfiles.keys()
    mtx = np.vstack( [biasProfiles[taxId][(region[0]):(region[1])] for taxId in keys] )
    df = pd.DataFrame( mtx )
    df.index = pd.Index(keys)
    meanValues = df.mean(axis=1)
    sortedValues = meanValues.sort_values()

    topN    = sortedValues[-N: ]
    bottomN = sortedValues[  :N]

    if len(biasProfiles)>=N:
        assert(len(topN)==N)
        assert(len(bottomN)==N)
    else:
        assert(len(topN)==len(biasProfiles))
        assert(len(bottomN)==len(biasProfiles))

    return (topN.index.values, bottomN.index.values)
    

def createProfilesBoxplot( profileData, plotId, namedProfiles=None, highlightProfiles=None, showBoxplot=False, title=None ):

    if title is None:
        title = plotId
        
    fig, ax1 = plt.subplots(figsize=(8, 2.2))

    plt.grid(False)

    # Plot density using KDE
    #sns.kdeplot( plotdata[:,:,0].flatten(), plotdata[:,:,1].flatten(), n_levels=10, bw=0.65, cmap="Blues", shade=True, shade_lowest=True, legend=True, ax=ax1 )

    #sns.lvplot(data=x, ax=ax1, color="#50d080")
    if showBoxplot:
        sns.boxplot(data=profileData, color="#7070d0", linewidth=1.1, ax=ax1 ) #, boxprops={'zorder': 10})

    
    plt.fill_between( range(v), quantileValues[0,:], quantileValues[1,:], facecolor='#80b0f0', edgecolor='none', alpha=0.5 )
    plt.fill_between( range(v), quantileValues[1,:], quantileValues[3,:], facecolor='#1540f0', edgecolor='none', alpha=0.5, label="Q1-Q3" )
    plt.fill_between( range(v), quantileValues[3,:], quantileValues[4,:], facecolor='#80b0f0', edgecolor='none', alpha=0.5 )
    plt.plot(  range(v), quantileValues[2,:], linewidth=2, color='#1030f0',  alpha=0.8, label="Median" )

    p0 = mpatches.Patch(color='#80b0f0', label="1.5*IQR",  alpha=0.5)
    p1 = mpatches.Patch(color='#1540f0', label="Q1-Q3",    alpha=0.5)
    p2 = mpatches.Patch(color='#1030f0', label="Median",   alpha=0.8)
    plt.legend(handles=(p2,p1,p0))
    
    # plot outliers (values outside the 1.5*IQR interval)
    for pos in range(v):
        samples = profileData[:,pos]
        outliers = samples[samples > quantileValues[4,pos]]
        if len(outliers)>0:
            plt.scatter([pos]*len(outliers), outliers, marker="D", s=16.0, c="#202020", alpha=0.4)
        outliers = samples[samples < quantileValues[0,pos]]
        if len(outliers)>0:
            plt.scatter([pos]*len(outliers), outliers, marker="D", s=16.0, c="#202020", alpha=0.4)

    # Plot neutral line
    plt.plot([-1,v], [0,0], '-k', linewidth=1, zorder=1)


    # Plot window width
    #plt.plot([0.5, 0.5 + 40/profileStep], [-1.8, -1.8], '-r', linewidth=4)
    #plt.annotate( s="window", xy=(1.0, -1.70) )
    plt.plot([18.5, 18.5 + 40/profileStep], [2.3, 2.3], '-r', linewidth=4)
    plt.annotate( s="window", xy=(19.0, 2.40) )


    if not highlightProfiles is None:
        for taxId in highlightProfiles:
            profile = namedProfiles[taxId]
            plt.plot( range(len(profile)), profile, linewidth=2.8, alpha=0.5, color="#ff5030" )
        print("{} - {}".format(title, highlightProfiles))
            

    # Annotate number of items
    plt.annotate( s='n = %d' % u, xy=(10.0, 2.3) )
    plt.xlim([-0.5, v-0.5])
    plt.ylim([-1.5, 3.0])

    xticks = list(range(0, v+1, 2))
    ax1.set_xticks(xticks)
    ax1.set_xticklabels([str(x_*10) for x_ in xticks])
    ax1.set_yticks( (-1, 0, 1, 2, 3) )
    ax1.yaxis.grid()

    ax1.set_xlabel( "Window start position (nt)",   fontname="DejaVu Sans Mono" )  # The default font ("Bitstream sans") is missing many unicode ranges
    ax1.set_ylabel( u"\u0394LFE", fontname="DejaVu Sans Mono" )  

    plt.title(  title )
    #plt.ylabel( u"\u0394LFE" )
    #plt.xlabel( "Window start position (nt)" )


    #plt.title("lvplot")
    plt.savefig("find_typical_profile.out.{}.pdf".format(plotId))
    plt.close()


createProfilesBoxplot(profileData=x, plotId="8a", showBoxplot=True)
createProfilesBoxplot(profileData=x, plotId="8b", showBoxplot=False, title="")


(top, bottom) = findValueOutliers(biasProfiles, (0,2))
createProfilesBoxplot(profileData=x, plotId="8.0.2.top",      namedProfiles=biasProfiles, highlightProfiles=tuple(top),    title="Top 10 (0-10nt)")
createProfilesBoxplot(profileData=x, plotId="8.0.2.bottom",   namedProfiles=biasProfiles, highlightProfiles=tuple(bottom), title="Bottom 10 (0-10nt)")

(top, bottom) = findValueOutliers(biasProfiles, (11,32))
createProfilesBoxplot(profileData=x, plotId="8.11.32.top",    namedProfiles=biasProfiles, highlightProfiles=tuple(top),    title="Top 10 (100-300nt)")
createProfilesBoxplot(profileData=x, plotId="8.11.32.bottom", namedProfiles=biasProfiles, highlightProfiles=tuple(bottom), title="Bottom 10 (100-300nt)")





fig, ax1 = plt.subplots(figsize=(8, 2.2))
# Plot density using KDE
#sns.kdeplot( plotdata[:,:,0].flatten(), plotdata[:,:,1].flatten(), n_levels=10, bw=0.65, cmap="Blues", shade=True, shade_lowest=True, legend=True, ax=ax1 )

sns.lvplot(data=x, ax=ax1, color="#50d080")
#sns.boxplot(data=x, color="#50d080", ax=ax1, boxprops={'zorder': 10})


# Plot neutral line
plt.plot([-1,v], [0,0], '-k')

# Plot window width
plt.plot([0.5, 0.5 + 40/profileStep], [-1.9,-1.9], '-r', linewidth=4)
plt.annotate( s="window", xy=(1.0, -1.83) )

# Annotate number of items
plt.annotate( s='n = %d' % u, xy=(26.0, 1.8) )
plt.xlim([-0.5, v-0.5])
plt.ylim([-2, 2.5])


plt.title("lvplot")
plt.savefig("find_typical_profile.out.9.pdf")
plt.close()






# Save the typical (modal) profile values as a csv file
typicalProfilePd = pd.DataFrame( { "CDS_position": np.linspace(0, v, len(estimatedValues)),  "Profile": estimatedValues, "Optimal_KDE_BW": optimalBW } )
typicalProfilePd.to_csv( typicalProfileOutputCSV )

# Create a profile tile with the modal profile
print(getProfileHeatmapTile( 99990010, pd.DataFrame( { 99990010 : estimatedValues } ), fixedYrange, True, profileStep ))
# Create a scale for the modal profile (using the specified Y range)
print(getLegendHeatmapTile( fixedYrange ))





