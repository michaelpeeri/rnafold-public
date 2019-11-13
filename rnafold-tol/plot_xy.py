from __future__ import print_function
import sys
from random import choice
import numpy as np
import pandas as pd
from math import log10
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress, wilcoxon
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from sklearn import decomposition
import data_helpers
import species_selection_data
from mfe_plots import getHeatmaplotProfilesValuesRange, scatterPlot, loadProfileData, PCAForProfiles, plotMFEvsCUBcorrelation
from fit_profile_params import getEstimatedParams
from ncbi_entrez import getTaxonomicGroupForSpecies

def getTaxName(taxId):
    return data_helpers.getSpeciesFileName(taxId)



def assignColors(data):
    assignment = {}
    availColors = ('green', 'blue', 'cyan', 'red', 'yellow', 'white', 'grey', 'orange', 'pink')
    pos = 0
    for item in set(data):
        assignment[item] = availColors[pos]
        pos += 1

    colors = []
    labels = []
    for d in data:
        labels.append( d )
        colors.append( assignment[d] )

    return colors, labels
        
        
        

def plotXY(xvals, yvals, _labels, groups):
    fig, ax = plt.subplots()
    colors, labels = assignColors(groups)
    #print(colors, labels)

    # Linear correlation and factors
    pearson = pearsonr(xvals, yvals)
    spearman = spearmanr(xvals, yvals)
    kendall = kendalltau(xvals, yvals)
    l = linregress(xvals, yvals)

    abline_x = np.arange(0.25, 0.85, 0.1)
    abline_y = abline_x * l.slope + l.intercept
    plt.plot(abline_x, abline_y, '--')

    plt.plot([0.25,0.50,0.75], [0.0, 0.0, 0.0])

    topr = max(yvals)*1.05
    scaler = topr/20
    # plot the linear approximation
    plt.annotate(s="Pearson r: %1.3f (p<%g)"  % (pearson[0], pearson[1]),                 xy=(0.28, topr-scaler*1),  fontsize=6 )
    plt.annotate(s="Pearson r^2: %1.3f"  % (pearson[0]**2,),                              xy=(0.28, topr-scaler*2),  fontsize=6 )
    plt.annotate(s="Spearman r: %1.3f (p<%g)" % (spearman.correlation, spearman.pvalue),  xy=(0.28, topr-scaler*3),  fontsize=6 )
    plt.annotate(s="Kendall's tau: %1.3f (p<%g)" % (kendall.correlation, kendall.pvalue), xy=(0.28, topr-scaler*4),  fontsize=6 )
    plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(0.28, topr-scaler*5),  fontsize=6 )


    # Species scatter-plot
    for label in sorted(set(labels)):
        xdata = []
        ydata = []
        color = None
        for _x,_y,_c,_l in zip(xvals, yvals, colors, labels):
            if( _l == label ):
                xdata.append(_x)
                ydata.append(_y)
                color = _c
        if( xdata ):
            #print(xdata, ydata)
            plt.scatter(xdata, ydata, c=color, label=label, s=40)

    for x,y,l in zip(xvals, yvals, _labels):
        label = l[:4]
        plt.annotate(label, (x+0.003, y-1.0), fontsize=7)


    # plt.annotate('', xy=(0.60, 0.4), xycoords='data',
    #             xytext=(0, -0.35), textcoords='offset points', color='black',
    #             arrowprops=dict(arrowstyle="->")
    #             )
    # plt.annotate("Stronger", xy=(0.42, 0.4), xycoords='data', color='black', rotation='vertical')
    # plt.annotate("Stronger", xy=(0.42, 0.4), xycoords='data', color='black', rotation='vertical')

    # plt.annotate('', xy=(0.60, -0.4), xycoords='data',
    #             xytext=(0, 0.35), textcoords='offset points', color='black',
    #             arrowprops=dict(arrowstyle="->")
    #             )
    #axis = ax.new_floating_axis(0,1)
    #axis.label.set_text("Stronger")

    #axis = ax.new_floating_axis(0,-1)
    #axis.label.set_text("Weaker")

    plt.xlim([.25,.75])
    plt.xlabel('GC')
    plt.ylabel('Local MFE selection (log10 p-value, Wilcoxon)')
    plt.grid(True)
    plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
    plt.savefig("scatter_xy.pdf")
    plt.savefig("scatter_xy.svg")
    plt.close(fig)

def plotXY_2(xvals, yvals, _labels, groups):
    fig, ax = plt.subplots()
    colors, labels = assignColors(groups)

    for label in sorted(set(labels)):
        xdata = []
        ydata = []
        color = None
        for _x,_y,_c,_l in zip(xvals, yvals, colors, labels):
            if( _l == label ):
                xdata.append(_x)
                ydata.append(_y)
                color = _c
        if( xdata ):
            #print(xdata, ydata)
            plt.scatter(xdata, ydata, c=color, label=label, s=50)

    for x,y,l in zip(xvals, yvals, _labels):
        l = l[:4]
        plt.annotate(l, (x-0.006, y+0.02), fontsize=7)

    plt.xlabel('GC-content')
    plt.ylabel('Native LMFE')
    plt.grid(True)
    plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
    plt.savefig("scatter_xy_nativeonly.pdf")
    plt.savefig("scatter_xy_nativeonly.svg")
    plt.close(fig)


def plotXY_3(xvals, yvals1, yvals2, _labels, groups):
    fig, (ax1,ax2) = plt.subplots(1, 2, sharey=True) #, gridspec_kw={'width_ratios': [1, 1]})

    #fig, ax = plt.subplots()
    colors, labels = assignColors(groups)

    for label in sorted(set(labels)):
        xdata = []
        ydata1 = []
        ydata2 = []
        color = None
        for _x,_y1,_y2,_c,_l in zip(xvals, yvals1, yvals2, colors, labels):
            if( _l == label ):
                xdata.append(_x)
                ydata1.append(_y1)
                ydata2.append(_y2)
                color = _c
        if( xdata ):
            #print(xdata, ydata)
            ax1.scatter(xdata, ydata1, c=color, label=label, s=50)
            ax2.scatter(xdata, ydata2, c=color, label=label, s=50)
            #plt.scatter(xdata, ydata, c=color, label=label, s=80)

    for x,y1,y2,l in zip(xvals, yvals1, yvals2, _labels):
        #plt.annotate(l, (x-0.006, y+0.02))
        ax1.annotate(l, (x-0.006, y1+0.02), fontsize=4)
        ax2.annotate(l, (x-0.006, y2+0.02), fontsize=4)

    plt.xlabel('GC-content')
    ax1.set_ylabel('Native LMFE')
    ax2.set_ylabel('Shuffled LMFE')
    plt.grid(True)
    ax1.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
    plt.savefig("scatter_xy_split.pdf")
    plt.savefig("scatter_xy_split.svg")
    plt.close(fig)


def plotSpeciesVars(data, xvar, yvar, labelvar, groupvar, filename, title=None, showSpeciesLabels=True):
    fig, ax = plt.subplots()

    #data = data.copy()
    data = data.dropna(subset=[xvar,yvar], how='any', axis=0)
    #del data[data[yvar].isnull()]
    #print("///" * 20)
    #print(data)
    #print("///" * 20)
    #print(colors, labels)

    colors, labels = assignColors(data[groupvar])

    _xrange = (min(data[xvar]), max(data[xvar]))
    _yrange = (min(data[yvar]), max(data[yvar]))
    print(_xrange)

    # Linear correlation and factors
    pearson = pearsonr(data[xvar], data[yvar])
    spearman = spearmanr(data[xvar], data[yvar])
    kendall = kendalltau(data[xvar], data[yvar])
    l = linregress(data[xvar], data[yvar])
    abline_x = np.arange(_xrange[0], _xrange[1], 0.1)
    abline_y = abline_x * l.slope + l.intercept
    plt.plot(abline_x, abline_y, '--')

    # plot x-axis
    #plt.plot([_xrange[0],(_xrange[0]+_xrange[1])*0.5,_xrange[1]], [0.0, 0.0, 0.0])

    topr = max(data[yvar])*1.05
    scaler = (_yrange[1]-_yrange[0])*0.036
    plt.annotate(s="Pearson r: %1.3f (p<%g)"     % (pearson[0], pearson[1]),                xy=(_xrange[0], topr-scaler*1),  fontsize=6 )
    plt.annotate(s="Pearson r^2: %1.3f"          % (pearson[0]**2,),                        xy=(_xrange[0], topr-scaler*2),  fontsize=6 )
    plt.annotate(s="Spearman r: %1.3f (p<%g)"    % (spearman.correlation, spearman.pvalue), xy=(_xrange[0], topr-scaler*3),  fontsize=6 )
    plt.annotate(s="Kendall's tau: %1.3f (p<%g)" %  (kendall.correlation, kendall.pvalue),  xy=(_xrange[0], topr-scaler*4),  fontsize=6 )
    plt.annotate(s="n= %d"  % len(data[yvar]),                                              xy=(_xrange[0], topr-scaler*5),  fontsize=6 )


    # Species scatter-plot
    for label in sorted(set(labels)):
        xdata = []
        ydata = []
        color = None
        for _x,_y,_c,_l in zip(data[xvar], data[yvar], colors, labels):
            if( _l == label ):
                xdata.append(_x)
                ydata.append(_y)
                color = _c
        if( xdata ):
            #print(xdata, ydata)
            plt.scatter(xdata, ydata, c=color, label=label, s=20)

    if showSpeciesLabels:
        for x,y,label in zip(data[xvar], data[yvar], data[labelvar]):
            plt.annotate(label, (x+(_xrange[1]-_xrange[0])*0.01, y-(_yrange[1]-_yrange[0])*0.01), fontsize=7)

    w = (_xrange[1]-_xrange[0])*1.05
    plt.xlim((_xrange[1] - w, _xrange[0] + w))
    plt.xlabel(xvar)
    plt.ylabel(yvar)
    if not title is None:
        #plt.title(title)
        ax.annotate(title,
                    xy=(0.95, 0.95), xycoords='figure fraction',
                    xytext=(0.95,0.95), textcoords='figure fraction', horizontalalignment='right')
    plt.grid(True)
    plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
    plt.savefig("scatter_xy_%s.pdf" % filename)
    plt.savefig("scatter_xy_%s.svg" % filename)
    plt.close(fig)


def calcWilcoxonPvalue_method1(df):
    #difs = np.array(df.native - df.shuffled)[30:]
    #indices = [0, 29, 59, 89, 119, 150]
    difs = np.array(df.native - df.shuffled)
    direction = np.sign(np.mean(difs))

    pval = wilcoxon(difs).pvalue
    
    return log10(pval) * direction * -1

def calcWilcoxonPvalue_method2(df2):
    assert(df2.ndim==2)

    #print(df2['delta'].isnull().head())
    df2 = df2[~df2['delta'].isnull()]

    direction = np.sign(np.mean(df2['delta']))
    pval = wilcoxon(df2['delta']).pvalue

    

    #print(df2.shape)
    #print(df2.head())

    #print(pval)

    if( pval>0.0 ):
        return log10(pval) * direction * -1
    elif( pval==0.0):    # I think exact comparison to 0.0 is safe with floating point numbers
        return -320.0      * direction * -1
    else:
        assert(False)




def estimateProfileParams(xs, ys):
    return getEstimatedParams(xs, ys, [0.0, 20.0, 80.0], [1.0, 1.0, 8.0])



def PCAforParams(params):
    X = np.vstack(params.values())
    assert(len(X)==len(xdata))

    pca = decomposition.PCA()
    pca.fit(X)
    print("Explained variance (parametric PCA):")
    print(pca.explained_variance_)

    pca.n_components = 2
    X_reduced = pca.fit_transform(X)
    print(X_reduced.shape)

    fig, ax = plt.subplots()
    plt.scatter(X_reduced[:,0], X_reduced[:,1], label=[getTaxName(x)[:4] for x in biasProfiles.keys()])
    for i, taxId in enumerate(biasProfiles.keys()):
        x = X_reduced[i,0]
        y = X_reduced[i,1]
        label = getTaxName(taxId)[:4]
        plt.annotate(label, (x+0.05, y-0.05), fontsize=7)


    plt.ylim([-3,3])
    plt.xlim([-3,3])
    plt.grid(True)
    plt.savefig("pca_for_params_test.pdf")
    plt.savefig("pca_for_params_test.svg")
    plt.close(fig)

    a = list(zip(biasProfiles.keys(), X_reduced[:,0]))
    print("top:")
    print(a[:3])
    print("bottom:")
    print(a[-3:])
    return a


def standalone():
    import sys
    import os
    from glob import glob
    
    globs = sys.argv[1].split(",")

    files = [[x for x in glob(profilesGlob) if os.path.exists(x)] for profilesGlob in globs]
        
    for group in range(len(files)):
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files[group])
        if group==0:
            print(dfProfileCorrs)
            plotMFEvsCUBcorrelation( biasProfiles, dfProfileCorrs )
    
    return

    #print(dfProfileCorrs)

    summaryStatistics.sort_values(by='tax_group').to_csv("scatter_xy_summary.csv")
    summaryStatistics.sort_values(by='tax_group').to_html(buf=open("scatter_xy_summary.html", "w"))
    summaryStatistics.groupby('tax_group').sum().to_csv("scatter_xy_summary_sums.csv")

    plotXY(xdata, ydata, labels, groups)
    plotXY_2(xdata, ydata_nativeonly, labels, groups)
    plotXY_3(xdata, ydata_nativeonly, ydata_shuffledonly, labels, groups)
    print("%d files included" % filesUsed)

    profileParams = map(lambda x: x[0](x[1]), zip((int,int,str,int), files[0][:-3].split("_")[5:]))
    title = "Profile: %d-%d (step=%d)" % (profileParams[3],profileParams[0], profileParams[1])

    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'optimal_temperature', 'short_tax_name', 'tax_group', 'gc_vs_temperature',   title )
    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'optimal_temperature', 'short_tax_name', 'tax_group', 'gc_vs_temperature_nolabels',   title, showSpeciesLabels=False )
    
    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'mean_delta_lfe',      'short_tax_name', 'tax_group', 'gc_vs_dLFE',          title )
    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'mean_delta_lfe',      'short_tax_name', 'tax_group', 'gc_vs_dLFE_nolabels',          title, showSpeciesLabels=False  )

    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'optimal_temperature', 'short_tax_name', 'tax_group', 'dLFE_vs_temperature', title )
    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'optimal_temperature', 'short_tax_name', 'tax_group', 'dLFE_vs_temperature_nolabels', title, showSpeciesLabels=False  )

    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'paired_fraction',     'short_tax_name', 'tax_group', 'dLFE_vs_paired_fraction', title )
    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'paired_fraction',     'short_tax_name', 'tax_group', 'dLFE_vs_paired_fraction_nolabels', title, showSpeciesLabels=False  )

    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'paired_fraction',     'short_tax_name', 'tax_group', 'gc_vs_paired_fraction',   title )
    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'paired_fraction',     'short_tax_name', 'tax_group', 'gc_vs_paired_fraction_nolabels',   title, showSpeciesLabels=False  )


    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'gene_density',        'short_tax_name', 'tax_group', 'gc_vs_gene_density',   title )
    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'gene_density',        'short_tax_name', 'tax_group', 'gc_vs_gene_density_nolabels',   title, showSpeciesLabels=False  )

    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'gene_density',        'short_tax_name', 'tax_group', 'dLFE_vs_gene_density', title )
    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'gene_density',        'short_tax_name', 'tax_group', 'dLFE_vs_gene_density_nolabels', title, showSpeciesLabels=False  )

    plotSpeciesVars( summaryStatistics, 'gene_density',   'paired_fraction',     'short_tax_name', 'tax_group', 'gene_density_vs_paired_fraction',   title )
    plotSpeciesVars( summaryStatistics, 'gene_density',   'paired_fraction',     'short_tax_name', 'tax_group', 'gene_density_vs_paired_fraction_nolabels',   title, showSpeciesLabels=False  )
    

    order = PCAForProfiles(biasProfiles, xdata)

    orderMap = dict(order)
    #yrange = getHeatmaplotProfilesValuesRange(biasProfiles, dfProfileCorrs, lambda x: orderMap[x] )
    # TODO - restore plotting?

    print("---------")
    params = {}
    for k,v in biasProfiles.items():
        params[k] = estimateProfileParams(v.index, v.values)

    return 0



if __name__=="__main__":
    from sys import exit
    exit( standalone() )
