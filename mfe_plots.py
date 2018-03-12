import sys
import numpy as np
import pandas as pd
from math import log10, sqrt
from bisect import bisect_left
from scipy.stats import pearsonr, spearmanr, kendalltau, linregress, wilcoxon
import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from data_helpers import getSpeciesName, getSpeciesFileName, getGenomicGCContent, getSpeciesProperty, getSpeciesShortestUniqueNamesMapping
from sklearn import decomposition
import seaborn as sns
from pyqtree import Index
from ncbi_entrez import getTaxonomicGroupForSpecies
from rate_limit import RateLimit


def plotMFEProfileWithGC(taxId, profileId, data):
    fig, (ax1,ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    data[['native', 'shuffled']].plot(ax=ax1)
    data[['shuffled75', 'shuffled25']].plot(ax=ax1, style='--')

    speciesName = getSpeciesName(taxId)

    plt.title(speciesName)

    plt.xlabel('Position (nt, window start, from cds start)')

    ax1.set_title("Mean LFE for %s" % speciesName)
    ax1.set_ylabel('Mean LFE')
    ax1.legend(fontsize=8)
    ax1.grid(True)


    data['gc'].plot(ax=ax2)
    ax2.set_title("GC%")
    ax2.set_ylabel('GC% (in window)')
    ax2.grid(True)


    #profileId = "tbd" # str(args.profile.ProfileId).replace(':', '-')
    plt.savefig("mfe_v2_40nt_cds_%s_%s.pdf" % (profileId, getSpeciesFileName(taxId)) )
    plt.savefig("mfe_v2_40nt_cds_%s_%s.svg" % (profileId, getSpeciesFileName(taxId)) )
    plt.close(fig)


def plotMFEProfileMultiple(taxId, profileId, data, additionalVars, scaleBar=None):
    fig, axes = plt.subplots(2+len(additionalVars), sharex=True)

    data[['native', 'shuffled']].plot(ax=axes[0])
    #data[['shuffled75', 'shuffled25']].plot(ax=axes[0], style='--')

    smfe = data['native']-data['shuffled']
    axes[1].plot([min(data.index), max(data.index)], [0,0], c='black')
    axes[1].plot( data.index, smfe, zorder=10 )
    axes[1].set_ylabel(u"\u0394MFE")

    
    speciesName = getSpeciesName(taxId)

    plt.xlabel('Position (nt, window start, from cds start)')

    axes[0].set_title("Biases for %s" % speciesName)
    axes[0].set_ylabel(r'MFE')
    axes[0].legend(fontsize=8)
    axes[0].grid(True)

    plotRange = []
    for i, var in enumerate(additionalVars):
        currentAxis = i+2
        data[var].plot(ax=axes[currentAxis])
        #axes[currentAxis].set_title(var)
        axes[currentAxis].set_ylabel(var)
        axes[currentAxis].grid(True)

        plotRange = [min(data[var]), max(data[var])]

        yrange = (min(data[var]), max(data[var]))
        warn = ''
        if( len(data[var]) < 500 ):
            warn = ' <!>'

        spearman = spearmanr( smfe, data[var])
        axes[currentAxis].annotate(s="Signed Spearman r: %1.3f (p<%1.2f%s)" % (spearman.correlation, spearman.pvalue, warn),  xy=(max(data.index)*0.8, yrange[0]+(yrange[1]-yrange[0])*0.20),  fontsize=6 )

        spearman2 = spearmanr( abs(smfe), data[var])
        axes[currentAxis].annotate(s="Unsigned Spearman r: %1.3f (p<%1.2f%s)" % (spearman2.correlation, spearman2.pvalue, warn),  xy=(max(data.index)*0.8, yrange[0]+(yrange[1]-yrange[0])*0.05),  fontsize=6 )

    if( not scaleBar is None ):
        # Draw a scale-bar
        scaleBarPosY   = plotRange[0] + 0.75*(plotRange[1]-plotRange[0])  # use the range of the last plot
        scaleBarHeight =               0.075*(plotRange[1]-plotRange[0])  # half-height actually...
        
    
        axes[-1].plot([ 5, 5+scaleBar],          [scaleBarPosY, scaleBarPosY], c='black')
        axes[-1].plot([ 5, 5],                   [scaleBarPosY-scaleBarHeight, scaleBarPosY+scaleBarHeight], c='black')
        axes[-1].plot([ 5+scaleBar, 5+scaleBar], [scaleBarPosY-scaleBarHeight, scaleBarPosY+scaleBarHeight], c='black')

        

    #profileId = "tbd" # str(args.profile.ProfileId).replace(':', '-')
    plt.savefig("mfe_v2_40nt_cds_vars_%s_%s.pdf" % (profileId, getSpeciesFileName(taxId)) )
    plt.savefig("mfe_v2_40nt_cds_vars_%s_%s.svg" % (profileId, getSpeciesFileName(taxId)) )
    plt.close(fig)

    

def plotXY(taxId, profileId, data, xvar, yvar, title):
    fig, ax1 = plt.subplots()

    data[[yvar]].plot(ax=ax1)


    #plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(0.28, 17.0),  fontsize=6 )


    #plt.xlim([.25,.75])
    plt.xlabel(xvar.replace('_',' '))
    plt.ylabel(yvar.replace('_',' '))
    plt.title( title )
    plt.grid(True)
    #plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')

    plt.savefig("mfe_v2_40nt_%s_vs_%s_%s_%s.pdf" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))
    plt.savefig("mfe_v2_40nt_%s_vs_%s_%s_%s.svg" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))
    plt.close(fig)
    
    
def scatterPlot(taxId, profileId, data, xvar, yvar, title):
    data = data.copy()
    data.dropna(subset=(xvar, yvar), inplace=True)

    fig, ax1 = plt.subplots()
    data.plot(x=xvar, y=yvar, ax=ax1, kind='scatter')



    ################################

    xvals = data[xvar]
    yvals = data[yvar]
    #print(data.head())

    # Linear correlation and factors
    pearson = pearsonr(xvals, yvals)
    spearman = spearmanr(xvals, yvals)
    kendall = kendalltau(xvals, yvals)
    l = linregress(xvals, yvals)

    abline_x = np.arange(min(xvals), max(xvals)*1.1, (max(xvals)-min(xvals)) )
    abline_y = abline_x * l.slope + l.intercept
    plt.plot(abline_x, abline_y, '--')

    topr = max(yvals)*1.05
    left = min(xvals)
    scaler = topr/20
    # plot the linear approximation
    plt.annotate(s="Pearson $\\rho$: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                 xy=(left, topr-scaler*1),  fontsize=6 )
    plt.annotate(s="Pearson $r^2$: %1.3f"  % (pearson[0]**2,),                              xy=(left, topr-scaler*2),  fontsize=6 )
    plt.annotate(s="Spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue),  xy=(left, topr-scaler*3),  fontsize=6 )
    plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(left, topr-scaler*4),  fontsize=6 )
    plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(left, topr-scaler*5),  fontsize=6 )





    ################################3


    

    #plt.annotate(s="n= %d"  % len(data),                                                 xy=(0.3, 4),  fontsize=6 )
        

    #plt.xlim([.25,.75])
    plt.xlabel(xvar.replace('_',' '))
    plt.ylabel(yvar.replace('_',' '))
    plt.title( title % getSpeciesName(taxId) )
    plt.grid(True)
    #plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')

    plt.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_%s_%s.pdf" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))
    plt.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_%s_%s.svg" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))
    plt.close(fig)

class CenterPreservingNormlizer(matplotlib.colors.Normalize):
    def __init__(self, negativeRange, positiveRange):
        self._negativeRange = negativeRange
        assert(positiveRange>0.0)
        self._positiveRange = positiveRange
        assert(negativeRange<0.0)


        # Declare to the rest of the API what is our output range
        self.vmin=0
        self.vmax=1

    #def __call__(self, values, clip=None):
    def __call__(self, values):
        outHalfScale = 0.5*(float(self.vmax)-float(self.vmin))
        outCenter = self.vmin + outHalfScale

        out = values.copy()

        
        factor = self._positiveRange/outHalfScale
        #print("+factor: %g" % (factor))
        values[values > 0.0] /= factor
        factor = self._negativeRange/outHalfScale*-1
        #print("-factor: %g" % (factor))
        values[values <= 0.0] /= factor

        values += outCenter

        values = 1-values

        assert(np.all(values >= self.vmin))
        assert(np.all(values <= self.vmax))
        
        # Use the logistic function (https://en.wikipedia.org/wiki/Logistic_function) to increase 'contrast'
        # To test using gnuplot:
        # plot [0:1] 1/(1+exp(-15*(x-0.5)))
        #
        steepness = 15.0
        return 1/(1+np.exp(-steepness*(values-0.5)))


"""
Plot the profile for a single species (contained in 'data'), save into a file, and return its name.
 taxId   - taxId of the species to plot
 data    - profile data for multiple species (including the one specified by taxId...)
 yrange  - y-range scale (this allows all tiles to have the same scale)
"""
def getProfileHeatmapTile(taxId, data, yrange, ticks=False, profileStep=10, phylosignalProfiles=None):
    if not taxId in data:
        print("getProfileHeatmapTile(): taxId {} not found".format(taxId))
        return None

    assert(len(yrange)==2)
    assert(yrange[0] < yrange[1])
    
    #fig, ax = plt.subplots()
    fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    # read profile data
    series = data[taxId]
    
    cmapNormalizer        = CenterPreservingNormlizer(yrange[0], yrange[1])
    phylosignalNormalizer = None
    if not phylosignalProfiles is None:
        phylosignalNormalizer = CenterPreservingNormlizer(np.min(phylosignalProfiles.values), np.max(phylosignalProfiles.values) )

    # read phylosignal profile (if used)
    phylosignalProfile = None
    if( not phylosignalProfiles is None and taxId in phylosignalProfiles.index ):
        phylosignalProfile = phylosignalProfiles.loc[taxId,:]
        
        if( len(series) < len(phylosignalProfile) ):
            phylosignalProfile = phylosignalProfile[:len(series)]

        assert(phylosignalProfile.index[0] == "Profile.1")
        assert(len(phylosignalProfile)==len(series))
        #print( "Got phylosignal for {}".format(taxId) )


    imdata = np.array(series).reshape(1,-1)

    #ax.axis(xmin=series.index[0], xmax=series.index[-1])

    ax1.imshow( imdata, cmap='bwr', aspect='auto', norm=cmapNormalizer, interpolation="bilinear" )

    if( not phylosignalProfile is None ):
        ax2.imshow( phylosignalProfile.values.reshape(1,-1), cmap='bwr', aspect='auto', norm=phylosignalNormalizer, interpolation="bilinear" )


    #pos = list(ax1.get_position().bounds)
    #taxname = getSpeciesFileName(taxId)
    #taxDescriptor = "%s.%s" % (taxname[0], taxname[1:9])
    #fig.text(pos[0]-0.01, pos[1]+pos[3]/2., taxDescriptor, va='center', ha='right', fontsize=8)

    ax1.set_yticks(())
    ax2.set_yticks(())
    if ticks:
        tickValues = range(10, len(series)-10, 10)
        ax1.set_xticks(tickValues)
        ax1.set_xticklabels(["" for x in tickValues])
    else:
        ax1.set_xticks(())
            

    tileFilename = "heatmap_profile_taxid_%d.png" % taxId
    plt.savefig(tileFilename, orientation='portrait', bbox_inches='tight')
    #plt.savefig("heatmap_profile_taxid_%d.svg" % taxId, orientation='portrait')
    plt.close(fig)

    return tileFilename


"""
Plot the profile for a single species (contained in 'data'), save into a file, and return its name.
 taxId   - taxId of the species to plot
 data    - profile data for multiple species (including the one specified by taxId...)
 yrange  - y-range scale (this allows all tiles to have the same scale)
"""
def getLegendHeatmapTile(yrange):
    assert(len(yrange) == 2)
    assert(yrange[0] <= yrange[1])
    
    fig, ax = plt.subplots()

    series = np.linspace( 10**yrange[0], 1 - 10**(-yrange[1]), 100)  # Create a range whose logit image will cover the range yrange...
    #series = np.linspace( yrange[0], yrange[1], 100)   # linear scale
    series = np.log10(series/(1-series))  # Logit function (inverse if logistic function)
    #print(series)
    cmapNormalizer = CenterPreservingNormlizer(yrange[0], yrange[1])

    imdata = series
    imdata = np.vstack((imdata,imdata))  # pretty crude trick borrowed from matplotlib examples...

    ax.imshow( imdata, cmap='bwr', aspect=2.0, norm=cmapNormalizer, interpolation="bilinear" )

    def roundTowardZero(x):
        if x < -0.5:
            return round(x)+1
        elif x > 0.5:
            return round(x)-1
        else:
            return 0

    #ax.set_title(taxId)
    ax.set_yticks(())
    #ax.axis(xmin=yrange[0], xmax=yrange[1])
    tick_values = list(sorted(list(range(int(roundTowardZero(yrange[0])), int(roundTowardZero(yrange[1]))+1 )) + [-0.5, 0.5]))  # Put ticks at integer intervals, plus two special ticks at -.5 and .5 (since this regions spans much of the graph)
    #print(tick_values)
    tick_positions = [bisect_left(series, x) for x in tick_values]  # use bisect to find the approximate position of each tick
    #print(tick_positions)
    
    ax.set_xticks( tick_positions ) # set tick positions
    ax.set_xticklabels( ["%.2g" % x for x in tick_values], size="xx-large" )  # set tick labels

    tileFilename = "heatmap_profile_legend.png"
    plt.savefig(tileFilename, orientation='landscape', bbox_inches='tight', dpi=600)
    plt.close(fig)

    return tileFilename
    


    
def getHeatmaplotProfilesValuesRange(data, dummy1=None, dummy2=None):

    #keysInOrder = data.keys()[:]
    #if not order is None:
    #    keysInOrder.sort( key=lambda x:order(x) )
    keysInOrder = data.keys()
    

    # Find the overall range that will be used to normalize all values
    # Todo - allow some "over-exposure" to expand the mid range at the cost of losing detail at the extremes
    valuesRange = [10.0, -10.0]
    for taxId in keysInOrder:
        series = data[taxId]

        currMin = min(series)
        currMax = max(series)

        if(currMin < valuesRange[0] ):
            valuesRange[0] = currMin
        if(currMax > valuesRange[1] ):
            valuesRange[1] = currMax
    assert(valuesRange[0] < 0.0)
    assert(valuesRange[1] > 0.0)
    maxRange = max(-valuesRange[0], valuesRange[1])
    
    return (-maxRange, maxRange)  # return the normalized range used

    # fig, axes = plt.subplots(nrows=len(data), ncols=2, sharex='col') #, gridspec_kw={'width_ratios': [4,1]})
    # plt.grid(False)
    # cmapNormalizer = CenterPreservingNormlizer(-maxRange, maxRange)
    # for ax, taxId in zip([x[0] for x in axes], keysInOrder):
    #     series = data[taxId]

    #     imdata = np.array(series)
    #     imdata = np.vstack((imdata,imdata))  # pretty crude trick borrowed from matplotlib examples...

    #     #ax.axis(xmin=series.index[0], xmax=series.index[-1])

    #     ax.imshow( imdata, cmap='coolwarm', aspect='auto', norm=cmapNormalizer )
            
    #     pos = list(ax.get_position().bounds)

    #     taxname = getSpeciesFileName(taxId)
    #     taxDescriptor = "%s.%s" % (taxname[0], taxname[1:9])
    #     fig.text(pos[0]-0.01, pos[1]+pos[3]/2., taxDescriptor, va='center', ha='right', fontsize=8)
        
    #     #ax.set_title(taxId)
    #     ax.set_yticks(())
    #     #ax.tick_params

    # #plt.colorbar(im, ax=axes[0], norm=cmapNormalizer, orientation='horizontal')
    # #plt.colorbar(im, cax=axes[-1], orientation='horizontal')
    # cbarRange = np.expand_dims( np.linspace( valuesRange[0]+0.001, valuesRange[1], 200 ), 2)

    # cbx = plt.subplot2grid((len(data),2), (0,1), rowspan=len(data))
    # #cbx.imshow( np.hstack((cbarRange, cbarRange)), aspect='auto', cmap='coolwarm', norm=cmapNormalizer)
    # #cbx.set_xticks(())


    # corrDataForPlotting = corrData.rename( columns={'spearman_smfe_gc_rho':'GC%', 'spearman_smfe_Nc_rho':'ENc', 'spearman_smfe_CAI_rho':'CAI', 'spearman_smfe_Fop_rho':'Fop' } )
    # del corrDataForPlotting['spearman_smfe_gc_pval']
    # del corrDataForPlotting['spearman_smfe_Nc_pval']
    # del corrDataForPlotting['spearman_smfe_CAI_pval']
    # del corrDataForPlotting['spearman_smfe_Fop_pval']
    # orderdf = pd.DataFrame({'order':range(len(keysInOrder))}, index=keysInOrder)
    # df2 = pd.merge(corrDataForPlotting, orderdf, left_index=True, right_index=True, how='inner')
    # df2.sort_values(by=['order'])
    # del df2['order']

    # corrsHeatmap = sns.heatmap(df2, annot=True, fmt=".2g", ax=cbx)
    # #corrsHeatmap.savefig("heatmap_profile_correlations.pdf")
    # #corrsHeatmap.savefig("heatmap_profile_correlations.svg")

    # plt.savefig("heatmap_profile_test.pdf", orientation='portrait')
    # plt.savefig("heatmap_profile_test.svg", orientation='portrait')
    # plt.close(fig)

    # return (-maxRange, maxRange)  # return the normalized range used


def plotCorrelations(data, _labels, group_func=None, order=None):
    #fig, axes = plt.subplots(nrows=len(data), ncols=2, sharex='col') #, gridspec_kw={'width_ratios': [4,1]})

    keysInOrder = data.keys()[:]
    keysInOrder.sort( key=lambda x:order(x) )

    #map  = sns.heatmap()
    #keysInOrder = data.keys()

def scatterPlotWithColor(taxId, profileId, data, xvar, yvar, colorvar, title):
    fig, ax1 = plt.subplots()

    #data.plot(x=xvar, y=yvar, c=colorvar, size=3, ax=ax1, kind='scatter')
    data[(data[colorvar]>=0.25) & (data[colorvar]<=0.75)].plot(x=xvar, y=yvar, c='white', s=3, ax=ax1, kind='scatter')
    data[data[colorvar].isnull()].plot(x=xvar, y=yvar, c='white', s=3, ax=ax1, kind='scatter')
    data[data[colorvar]>0.75].plot(x=xvar, y=yvar, c='red', s=3, ax=ax1, kind='scatter')
    data[data[colorvar]<0.25].plot(x=xvar, y=yvar, c='blue', s=3, ax=ax1, kind='scatter')

    top20 = data.sort_values(by=colorvar, ascending=False).iloc[:20]
    for i in range(20):
        plt.annotate(s=top20.iloc[i]['protid'], xy=( top20.iloc[i][xvar], top20.iloc[i][yvar]), fontsize=2 )

    top20.plot(x=xvar, y=yvar, c='orange', s=10, alpha=0.7, ax=ax1, kind='scatter')

    
    bottom20 = data.sort_values(by=colorvar, ascending=False).iloc[-20:]
    for i in range(20):
        plt.annotate(s=bottom20.iloc[i]['protid'], xy=( bottom20.iloc[i][xvar], bottom20.iloc[i][yvar]), fontsize=2 )

    bottom20.plot(x=xvar, y=yvar, c=(0.2,0.7,1.0), s=10, alpha=0.7, ax=ax1, kind='scatter')
    

    plt.annotate(s="n= %d"  % len(data),                                                 xy=(0.3, 4),  fontsize=6 )
        

    #plt.xlim([.25,.75])
    plt.xlabel(xvar.replace('_',' '))
    plt.ylabel(yvar.replace('_',' '))
    plt.title( title % getSpeciesName(taxId) )
    plt.grid(True)
    #plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')

    plt.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_with_%s_%s_%s.pdf" % (yvar, xvar, colorvar, profileId, getSpeciesFileName(taxId)))
    plt.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_with_%s_%s_%s.svg" % (yvar, xvar, colorvar, profileId, getSpeciesFileName(taxId)))
    plt.close(fig)
    



def plotMFEProfileByPA(taxId, profileId, data):
    fig, ax1 = plt.subplots()

    data[['native', 'shuffled']].plot(ax=ax1)
    data[['native_pa_low', 'native_pa_med', 'native_pa_high']].plot(ax=ax1)
    data[['shuffled75', 'shuffled25']].plot(ax=ax1, style='--')

    speciesName = getSpeciesName(taxId)

    plt.title(speciesName)

    plt.xlabel('Position (nt, window start, from cds start)')

    ax1.set_title("Mean LFE for %s" % speciesName)
    ax1.set_ylabel('Mean LFE')
    ax1.legend(fontsize=8)
    ax1.grid(True)


    #profileId = "tbd" # str(args.profile.ProfileId).replace(':', '-')
    plt.savefig("mfe_v2_40nt_cds_%s_%s_by_pa.pdf" % (profileId, getSpeciesFileName(taxId)) )
    plt.savefig("mfe_v2_40nt_cds_%s_%s_by_pa.svg" % (profileId, getSpeciesFileName(taxId)) )
    plt.close(fig)


def scatterPlotWithKernel(taxId, profileId, data, xvar, yvar, title):
    g = sns.jointplot( xvar, yvar, data=data, kind="kde", dropna=True, ylim=(-5 ,5), xlim=(0.3, 0.7), gridsize=100 )
    g.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_%s_%s.pdf" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))
    g.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_%s_%s.svg" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))


def scatterPlotWithKernel2(taxId, profileId, data, xvar, yvar, title):
    data = data.copy()
    data.dropna(subset=(xvar, yvar), inplace=True)

    fig, ax1 = plt.subplots()
    #data.plot(x=xvar, y=yvar, ax=ax1, kind='scatter')
    ax = sns.kdeplot( data[xvar], data[yvar], cmap="Blues", shade=True, shade_lowest=True, legend=True)
    blue = sns.color_palette("Blues")[-2]



    ################################

    xvals = data[xvar]
    yvals = data[yvar]
    #print(data.head())

    # Linear correlation and factors
    pearson = pearsonr(xvals, yvals)
    spearman = spearmanr(xvals, yvals)
    kendall = kendalltau(xvals, yvals)
    l = linregress(xvals, yvals)

    min_xvals = 0.4
    max_xvals = 0.7
    abline_x = np.arange(min_xvals, max_xvals*1.01, (max_xvals-min_xvals) )

    abline_y = abline_x * l.slope + l.intercept
    plt.plot(abline_x, abline_y, '--')

    #topr = max(yvals)*1.05
    topr = 100
    #left = min(xvals)
    left = min_xvals
    scaler = topr/20
    # plot the linear approximation
    plt.annotate(s="Pearson r: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                 xy=(left, topr-scaler*1),  fontsize=6 )
    plt.annotate(s="Pearson $r^2$: %1.3f"  % (pearson[0]**2,),                              xy=(left, topr-scaler*2),  fontsize=6 )
    plt.annotate(s="Spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue),  xy=(left, topr-scaler*3),  fontsize=6 )
    plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(left, topr-scaler*4),  fontsize=6 )
    plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(left, topr-scaler*5),  fontsize=6 )


    plt.xlim([0.5,0.7])
    plt.ylim([-10,10])
    



    ################################3


    

    #plt.annotate(s="n= %d"  % len(data),                                                 xy=(0.3, 4),  fontsize=6 )
        

    #plt.xlim([.25,.75])
    plt.xlabel(xvar.replace('_',' '))
    plt.ylabel(yvar.replace('_',' '))
    plt.title( title % getSpeciesName(taxId) )
    plt.grid(True)
    #plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')

    plt.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_%s_%s.pdf" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))
    plt.savefig("mfe_v2_40nt_genelevel_%s_vs_%s_%s_%s.svg" % (yvar, xvar, profileId, getSpeciesFileName(taxId)))
    plt.close(fig)


short_names = set()
def shortenTaxName(name):
    currLength=4
    
    if name.startswith("Candidatus "): # drop 'Candidatus' prefix
        name = name[11:]
        
    while(currLength <= len(name)):
        candidate = name[:currLength]
        if not candidate in short_names:
            short_names.add(candidate)
            return candidate
        currLength += 1

    #raise Exception("Failed to shorten name '%s'" % name)
    
    # Try adding numerical indices to resolve ambiguities
    idx=1
    while(True):
        candidate = "%s%d" % (name[:4], idx)
        if not candidate in short_names:
            short_names.add(candidate)
            return candidate
        idx += 1

def calcWilcoxonPvalue_method1(df):
    difs = np.array(df.native - df.shuffled)
    direction = np.sign(np.mean(difs))

    pval = wilcoxon(difs).pvalue
    
    return log10(pval) * direction * -1

def calcWilcoxonPvalue_method2(df2):
    assert(df2.ndim==2)

    df2 = df2[~df2['delta'].isnull()]

    direction = np.sign(np.mean(df2['delta']))
    pval = wilcoxon(df2['delta']).pvalue

    if( pval>0.0 ):
        return log10(pval) * direction * -1
    elif( pval==0.0):    # I think exact comparison to 0.0 is safe with floating point numbers
        return -320.0      * direction * -1
    else:
        assert(False)




def getShortTaxName(taxId):
    return getSpeciesFileName(taxId)

def loadProfileData(files):
    xdata = []
    ydata = []
    ydata_nativeonly = []
    ydata_shuffledonly = []
    labels = []
    groups = []
    filesUsed = 0
    biasProfiles = {}

    dfProfileCorrs = pd.DataFrame( { "spearman_smfe_gc_rho":   pd.Series(dtype='float'),
                                     "spearman_smfe_gc_pval":  pd.Series(dtype='float'),
                                     "spearman_smfe_Nc_rho":   pd.Series(dtype='float'),
                                     "spearman_smfe_Nc_pval":  pd.Series(dtype='float'),
                                     "spearman_smfe_CAI_rho":  pd.Series(dtype='float'),
                                     "spearman_smfe_CAI_pval": pd.Series(dtype='float'),
                                     "spearman_smfe_Fop_rho":  pd.Series(dtype='float'),
                                     "spearman_smfe_Fop_pval": pd.Series(dtype='float') } )

    summaryStatistics = pd.DataFrame({
        'tax_name':pd.Series(dtype='string'),
        'short_tax_name':pd.Series(dtype='string'),
        'tax_id':pd.Series(dtype='int'),
        'genomic_gc':pd.Series(dtype='float'),
        'tax_group':pd.Series(dtype='string'), # TODO: change to categorical data; Categorical([], categories=('Bacteria', 'Archaea', 'Fungi', 'Plants'), ordered=False)
        'CDSs_included':pd.Series(dtype='int'),
        'profileElements':pd.Series(dtype='int'),
        'optimal_temperature':pd.Series(dtype='float'),
        'temperature_range':pd.Categorical([]),
        'mean_delta_lfe':pd.Series(dtype='float'),
        'paired_fraction':pd.Series(dtype='float'),
        'gene_density':pd.Series(dtype='float')
    })

    rl = RateLimit(10)

    for h5 in files:
        with pd.io.pytables.HDFStore(h5) as store:
            for key in store.keys():
                if key[:4] != "/df_":
                    continue

                dfHeader = key.split('_')
                taxId = int(dfHeader[1])
                taxName = getShortTaxName(taxId)
                #taxGroup = data_helpers.getSpeciesTaxonomicGroup(taxId)
                taxGroup = getTaxonomicGroupForSpecies(taxId)
                longTaxName = getSpeciesName(taxId)
                shortTaxName = shortenTaxName(taxName)

                df = store[key]
                df = df.iloc[:-1]  # remove the last value (which is missing)

                deltas_df = store["/deltas_"+key[4:]]
                genes_df = store["/deltas_"+key[4:]]
                summary_df = store["/statistics_"+key[4:]]
                profileCorrelations_df = store["/profiles_spearman_rho_"+key[4:]]


                df['MFEbias'] = pd.Series(df['native']-df['shuffled'], index=df.index)
                dfMFEbias = df['MFEbias']

                biasProfiles[taxId] = dfMFEbias

                meanDeltaLFE = np.mean(dfMFEbias)

                cdsCount = int(summary_df.iloc[0]['cds_count'])
                assert(cdsCount >= 100)
                #meanGC = species_selection_data.findByTaxid(taxId).iloc[0]['GC% (genome)']
                meanGC = getGenomicGCContent(taxId)  # this is actually the genomic GC% (not CDS only)

                # Fetch temperature data for this species (if available)
                optimalTemperatureData = getSpeciesProperty( taxId, 'optimum-temperature')
                optimalTemperature = None
                if not optimalTemperatureData[0] is None:
                    optimalTemperature = float(optimalTemperatureData[0])

                temperatureRangeData = getSpeciesProperty( taxId, 'temperature-range')
                temperatureRange = None
                if not temperatureRangeData[0] is None:
                    temperatureRange = temperatureRangeData[0]
                else:
                    temperatureRange = "Unknown"

                pairedFractionData = getSpeciesProperty( taxId, 'paired-mRNA-fraction')
                pairedFraction = None
                if not pairedFractionData[0] is None:
                    pairedFraction = float(pairedFractionData[0])

                    
                genomeSizeData = getSpeciesProperty( taxId, 'genome-size-mb')
                genomeSize = None
                if not genomeSizeData[0] is None:
                    genomeSize = float(genomeSizeData[0])

                proteinCountData = getSpeciesProperty( taxId, 'protein-count')
                proteinCount = None
                if not proteinCountData[0] is None:
                    proteinCount = int(proteinCountData[0])

                geneDensity = None
                if( (not genomeSize is None) and (not proteinCount is None)  ):
                    geneDensity = float(proteinCount)/genomeSize

                    
                summaryStatistics = summaryStatistics.append(pd.DataFrame({
                    'tax_name':pd.Series([taxName]),
                    'short_tax_name':pd.Series([shortTaxName]),
                    'long_tax_name':pd.Series([longTaxName]),
                    'tax_id':pd.Series([taxId], dtype='int'),
                    'genomic_gc':pd.Series([meanGC]),
                    'tax_group':pd.Series([taxGroup]),
                    'CDSs_included':pd.Series([cdsCount], dtype='int'),
                    'optimal_temperature':pd.Series([optimalTemperature], dtype='float'),
                    'temperature_range':pd.Categorical([temperatureRange]),
                    'mean_delta_lfe':pd.Series([meanDeltaLFE], dtype='float'),
                    'paired_fraction':pd.Series([pairedFraction], dtype='float'),
                    'gene_density':pd.Series([geneDensity], dtype='float')
                }))

                dfProfileCorrs = dfProfileCorrs.append( profileCorrelations_df )

                # Format:

                #         gc  native  position  shuffled
                # 1    0.451  -4.944         1    -5.886
                # 2    0.459  -5.137         2    -6.069
                # 3    0.473  -5.349         3    -6.262
                filesUsed += 1

                #print(df.shape)

                meanGC = np.mean(df.gc)
                xdata.append(meanGC)

                #meanE = np.mean(df.native - df.shuffled)
                #ydata.append(meanE)

                dirpval = calcWilcoxonPvalue_method2(deltas_df)
                #print(dirpval)s
                ydata.append(dirpval)

                #print(df.head())

                #print(df.native)
                #print(df.shuffled)

                meanE_nativeonly = np.mean(df.native)
                ydata_nativeonly.append(meanE_nativeonly)

                meanE_shuffledonly = np.mean(df.shuffled)
                ydata_shuffledonly.append(meanE_shuffledonly)

                labels.append( taxName )

                #groups.append( choice(('Bacteria', 'Archaea', 'Fungi', 'Plants')) )   # Testing only!!!
                groups.append( taxGroup )

                if( rl() ):
                    print("Loaded %d profiles (%.2g%%)" % (filesUsed, float(filesUsed)/len(files)*100))

    return (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics)

def loadPhylosignalProfiles( phylosignalFile ):
    df = None
    with open(phylosignalFile, 'r') as csvfile:
        df = pd.read_csv(csvfile, sep=',', dtype={"":np.int}, index_col=0 )
    return df


def getTaxName(taxId):
    return getSpeciesFileName(taxId)


_shortTaxNames = None
def getSpeciesShortestUniqueNamesMapping_memoized():
    global _shortTaxNames
    if( _shortTaxNames is None ):
        _shortTaxNames = getSpeciesShortestUniqueNamesMapping()
    return _shortTaxNames


def PCAForProfiles(biasProfiles, profileValuesRange, profilesYOffsetWorkaround=0.0, profileScale=1.0, fontSize=7, overlapAction="ignore", showDensity=True, highlightSpecies=None):
    X = np.vstack(biasProfiles.values()) # convert dict of vectors to matrix

    shortNames = getSpeciesShortestUniqueNamesMapping_memoized()

    pca = decomposition.PCA()
    pca.fit(X)

    pca.n_components = 2  # force 2 components
    X_reduced = pca.fit_transform(X)
    print(X_reduced.shape)

    # Assign dimensions to plot axes
    # TODO - test the other values...
    D0 = 0 # D0 - component to show on Y scale 
    D1 = 1 # D1 - component to show on X scale
    assert(D0!=D1)
    

    #fig = matplotlib.figure.Figure(suppressComposite=True)
    #ax = fig.subplots()
    fig, ax = plt.subplots()    

    cmapNormalizer = CenterPreservingNormlizer(profileValuesRange[0], profileValuesRange[1])


    #plt.scatter(X_reduced[:,1], X_reduced[:,0], label=[shortNames[x] for x in biasProfiles.keys()])

    # Calculate general scaling constants
    scaleX = max(X_reduced[:,D1]) - min(X_reduced[:,D1])
    scaleY = max(X_reduced[:,D0]) - min(X_reduced[:,D0])
    ori1 =   max(X_reduced[:,D0])
    ori2 =   max(X_reduced[:,D1])
    #assert(scaleY>scaleX)  # We map the first PCA dimension to the Y axis (since portrait format better fits the figure layout)
    scaleX = scaleY
    imw = scaleX*0.100*profileScale
    imh = scaleY*0.012*profileScale
    imofy = profilesYOffsetWorkaround # 0.0 # 0.45 # TODO - FIX THIS

    spIndex = Index(bbox=(min(X_reduced[:,D1]), min(X_reduced[:,D0]), max(X_reduced[:,D1]), max(X_reduced[:,D0]) ))
    
    #zorders = [x[0] for x in sorted(list(enumerate([sqrt((X_reduced[i,D1]-ori2)**2 + (X_reduced[i,D0]-ori1)**2) for i in range(len(X_reduced))])), key=lambda x:x[D1])]

    dinfo = {}

    highlightPatches = []

    # Plot each profile 
    for i, taxId in enumerate(biasProfiles.keys()):
        # Determine the location for this profile
        x = X_reduced[i,D1]
        y = X_reduced[i,D0]

        showProfile = True
        label = shortNames[taxId]
        
        # is there overlap?
        if overlapAction=="hide":
            bbox = (x-imw, y-imh, x+imw, y+imh)
            print("---"*10)
            print("new: {} {}".format(label, bbox))
            matches = spIndex.intersect( bbox )
            for m in matches:
                print("-X- {} {}".format(m, dinfo[m]))
                
            if matches and taxId not in highlightSpecies:
                print( "Hiding profile at {}".format( (x,y) ) )
                showProfile = False
            else:
                spIndex.insert( label, bbox )
                dinfo[label] = bbox
            
        # plot the profile
        #zorder = zorders[i]
        if showProfile:      # skip this if we decided to hide this profile

            # Show the text label
            if fontSize > 0:
                plt.annotate(label, (x - scaleX*0.03, y + scaleY*0.012), fontsize=fontSize, zorder=5)

            if taxId in highlightSpecies:
                rect = Rectangle((bbox[0]-0.03, bbox[1]+0.03+imh), imw*2+0.06, imh*2+0.06, edgecolor="blue", linewidth=4.0, facecolor="none", fill=False, linestyle="solid")
                highlightPatches.append(rect)
                

            # Show the profile tile
            if True:
                plt.imshow( np.array( biasProfiles[taxId] ).reshape(1,-1), cmap='bwr', norm=cmapNormalizer, extent=(x-imw, x+imw, -y-imh+imofy, -y+imh+imofy ), interpolation='bilinear', zorder=2005 )
        

    # Paint the highlights
    if False:
        ax.add_collection(PatchCollection(highlightPatches))
        
    ax.set_aspect("equal")
    centerY = (max(X_reduced[:,D0]) + min(X_reduced[:,D0])) * 0.5
    centerX = (max(X_reduced[:,D1]) + min(X_reduced[:,D1])) * 0.5
    ax.set_xlim(( centerX - scaleY*0.55, centerX + scaleY*0.55 ))   # Use same scale for both axes
    ax.set_ylim(( centerY - scaleY*0.55, centerY + scaleY*0.55 ))

    if False:
        tlx = max(X_reduced[:,D1])
        tly = min(X_reduced[:,D0])
        bry = max(X_reduced[:,D0])
        cmapNormalizerForPCAvars = CenterPreservingNormlizer(-1, 1)
        plt.imshow( pca.components_[0,:].reshape(1,-1), extent=(tlx-imw*2, tlx,  tly-imh*2, tly-imh*4 ), norm=cmapNormalizerForPCAvars, cmap='coolwarm', interpolation='bilinear', zorder=5 )
        plt.imshow( pca.components_[1,:].reshape(1,-1), extent=(tlx-imw*2, tlx,  tly      , tly-imh*2 ), norm=cmapNormalizerForPCAvars, cmap='coolwarm', interpolation='bilinear', zorder=5 )
        plt.annotate( "V1",                                               (tlx - imw*2.4, bry+imh*2), fontsize=fontSize )
        plt.annotate( "V2",                                               (tlx - imw*2.4, bry),       fontsize=fontSize )
        plt.annotate( "{:.3g}".format( pca.explained_variance_ratio_[0] ), (tlx + imw*0.1, bry+imh*2), fontsize=fontSize )
        plt.annotate( "{:.3g}".format( pca.explained_variance_ratio_[1] ), (tlx + imw*0.1, bry      ), fontsize=fontSize )

    if False:
        loadingScale = 5.4
        plt.annotate( u"\u0394LFE[0nt]",   xy=(0,0), xytext=(pca.components_[D1,0 ]*loadingScale, pca.components_[D0,0 ]*loadingScale), fontsize=fontSize, arrowprops=dict(arrowstyle='<-', alpha=0.4, linewidth=1.5, color='red'), color='red', zorder=200 )
        plt.annotate( u"\u0394LFE[100nt]", xy=(0,0), xytext=(pca.components_[D1,10]*loadingScale, pca.components_[D0,10]*loadingScale), fontsize=fontSize, arrowprops=dict(arrowstyle='<-', alpha=0.4, linewidth=1.5, color='red'), color='red', zorder=200 )


    #plt.scatter(X_reduced[:,1], X_reduced[:,0] )
    
    
    plt.grid(False)
    ax.set_xticks(())
    ax.set_yticks(())
    #plt.xlim( (min(X_reduced[:,D1]) - scaleX*0.1, max(X_reduced[:,D1]) + scaleX*0.1) )
    #plt.ylim( (min(X_reduced[:,D0]) - scaleY*0.1, max(X_reduced[:,D0]) + scaleY*0.1) )

    plt.savefig("pca_profiles.pdf")
    plt.savefig("pca_profiles.png", dpi=300, transparent=True)
    plt.savefig("pca_profiles.svg")
    plt.close(fig)


    # Create separate density plot
    fig, ax = plt.subplots()

    if showDensity:
        sns.kdeplot( X_reduced[:,D1].flatten(), X_reduced[:,D0].flatten(), n_levels=20, cmap="Blues", shade=True, shade_lowest=True, legend=False, ax=ax, zorder=1 )

    #plt.scatter(X_reduced[:,1], X_reduced[:,0] )

    ax.set_aspect("equal")
    centerY = (max(X_reduced[:,D0]) + min(X_reduced[:,D0])) * 0.5
    centerX = (max(X_reduced[:,D1]) + min(X_reduced[:,D1])) * 0.5
    ax.set_xlim(( centerX - scaleY*0.55, centerX + scaleY*0.55 ))   # Use same scale for both axes
    ax.set_ylim(( centerY - scaleY*0.55, centerY + scaleY*0.55 ))
    
    plt.grid(False)
    ax.set_xticks(())
    ax.set_yticks(())
    #plt.xlim( (min(X_reduced[:,D1]) - scaleX*0.1, max(X_reduced[:,D1]) + scaleX*0.1) )
    #plt.ylim( (min(X_reduced[:,D0]) - scaleY*0.1, max(X_reduced[:,D0]) + scaleY*0.1) )

    
    plt.savefig("pca_profiles_density.pdf")
    plt.savefig("pca_profiles_density.png", dpi=300)
    plt.savefig("pca_profiles_density.svg")
    plt.close(fig)
    
    
    print("Explained variance: {}".format(pca.explained_variance_ratio_))
    print("            (Total: {})".format(sum(pca.explained_variance_ratio_)))
    print("Components: {}".format(pca.components_[:2,:]))
    
    a = list(zip(biasProfiles.keys(), X_reduced[:,D0]))
    a.sort(key=lambda x:x[1])
    print("top:")
    print(a[:3])
    print("bottom:")
    print(a[-3:])
    return a


#def drawPCAforTree( profilesAsMap, xdata ):
#    arbitraryProfile =  profilesAsMap.values()[0]
#    assert(arbitraryProfile.ndim == 1)
#    profileLength = arbitraryProfile.shape[0]
#    numProfiles = len(profilesAsMap)#
#
#    profilesArray = np.zeros((numProfiles, profileLength))
#    for i, profile in enumerate(profilesAsMap.values()):
#        profilesArray[i] = profile#
#
#    PCAForProfiles( profilesArray, xdata )

