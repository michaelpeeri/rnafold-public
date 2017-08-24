import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, kendalltau, linregress
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from data_helpers import getSpeciesName, getSpeciesFileName
import seaborn as sns


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


def getProfileHeatmapTile(taxId, data, corrData, yrange):
    if not taxId in data:
        return None

    assert(len(yrange)==2)
    assert(yrange[0] < yrange[1])
    
    fig, ax = plt.subplots()

    series = data[taxId]
    cmapNormalizer = CenterPreservingNormlizer(yrange[0], yrange[1])

    imdata = np.array(series)
    imdata = np.vstack((imdata,imdata))  # pretty crude trick borrowed from matplotlib examples...

    #ax.axis(xmin=series.index[0], xmax=series.index[-1])

    ax.imshow( imdata, cmap='coolwarm', aspect='auto', norm=cmapNormalizer )

    pos = list(ax.get_position().bounds)

    #taxname = getSpeciesFileName(taxId)
    #taxDescriptor = "%s.%s" % (taxname[0], taxname[1:9])
    #fig.text(pos[0]-0.01, pos[1]+pos[3]/2., taxDescriptor, va='center', ha='right', fontsize=8)

    #ax.set_title(taxId)
    ax.set_yticks(())
    #ax.tick_params

    tileFilename = "heatmap_profile_taxid_%d.svg" % taxId
    plt.savefig(tileFilename, orientation='portrait', bbox_inches='tight')
    #plt.savefig("heatmap_profile_taxid_%d.svg" % taxId, orientation='portrait')
    plt.close(fig)

    return tileFilename
    
    

def heatmaplotProfiles(data, dummy1, corrData, dummy2, dummy3=None, order=None):
    fig, axes = plt.subplots(nrows=len(data), ncols=2, sharex='col') #, gridspec_kw={'width_ratios': [4,1]})

    keysInOrder = data.keys()[:]
    if not order is None:
        keysInOrder.sort( key=lambda x:order(x) )
    #keysInOrder = data.keys()
    
    plt.grid(False)

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
    
    cmapNormalizer = CenterPreservingNormlizer(-maxRange, maxRange)

    for ax, taxId in zip([x[0] for x in axes], keysInOrder):
        series = data[taxId]

        imdata = np.array(series)
        imdata = np.vstack((imdata,imdata))  # pretty crude trick borrowed from matplotlib examples...

        #ax.axis(xmin=series.index[0], xmax=series.index[-1])

        ax.imshow( imdata, cmap='coolwarm', aspect='auto', norm=cmapNormalizer )
            
        pos = list(ax.get_position().bounds)

        taxname = getSpeciesFileName(taxId)
        taxDescriptor = "%s.%s" % (taxname[0], taxname[1:9])
        fig.text(pos[0]-0.01, pos[1]+pos[3]/2., taxDescriptor, va='center', ha='right', fontsize=8)
        
        #ax.set_title(taxId)
        ax.set_yticks(())
        #ax.tick_params

    #plt.colorbar(im, ax=axes[0], norm=cmapNormalizer, orientation='horizontal')
    #plt.colorbar(im, cax=axes[-1], orientation='horizontal')
    cbarRange = np.expand_dims( np.linspace( valuesRange[0]+0.001, valuesRange[1], 200 ), 2)

    cbx = plt.subplot2grid((len(data),2), (0,1), rowspan=len(data))
    #cbx.imshow( np.hstack((cbarRange, cbarRange)), aspect='auto', cmap='coolwarm', norm=cmapNormalizer)
    #cbx.set_xticks(())


    corrDataForPlotting = corrData.rename( columns={'spearman_smfe_gc_rho':'GC%', 'spearman_smfe_Nc_rho':'ENc', 'spearman_smfe_CAI_rho':'CAI', 'spearman_smfe_Fop_rho':'Fop' } )
    del corrDataForPlotting['spearman_smfe_gc_pval']
    del corrDataForPlotting['spearman_smfe_Nc_pval']
    del corrDataForPlotting['spearman_smfe_CAI_pval']
    del corrDataForPlotting['spearman_smfe_Fop_pval']
    orderdf = pd.DataFrame({'order':range(len(keysInOrder))}, index=keysInOrder)
    df2 = pd.merge(corrDataForPlotting, orderdf, left_index=True, right_index=True, how='inner')
    df2.sort_values(by=['order'])
    del df2['order']

    corrsHeatmap = sns.heatmap(df2, annot=True, fmt=".2g", ax=cbx)
    #corrsHeatmap.savefig("heatmap_profile_correlations.pdf")
    #corrsHeatmap.savefig("heatmap_profile_correlations.svg")

    plt.savefig("heatmap_profile_test.pdf", orientation='portrait')
    plt.savefig("heatmap_profile_test.svg", orientation='portrait')
    plt.close(fig)

    return (-maxRange, maxRange)  # return the normalized range used


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