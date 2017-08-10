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
from mfe_plots import heatmaplotProfiles, getProfileHeatmapTile, scatterPlot
from fit_profile_params import getEstimatedParams
from ncbi_entrez import getTaxonomicGroupForSpecies

def getTaxName(taxId):
    return data_helpers.getSpeciesFileName(taxId)


short_names = set()
# TODO: This t
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
    raise Exception("Failed to shorten name '%s'" % name)


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


def plotSpeciesVars(data, xvar, yvar, labelvar, groupvar, filename, title=None):
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
            plt.scatter(xdata, ydata, c=color, label=label, s=40)

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
    elif( pval==0.0):
        return -320.0      * direction * -1
    else:
        assert(False)



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

    for h5 in files:
        with pd.io.pytables.HDFStore(h5) as store:
            for key in store.keys():
                if key[:4] != "/df_":
                    continue

                dfHeader = key.split('_')
                taxId = int(dfHeader[1])
                taxName = getTaxName(taxId)
                #taxGroup = data_helpers.getSpeciesTaxonomicGroup(taxId)
                taxGroup = getTaxonomicGroupForSpecies(taxId)
                longTaxName = data_helpers.getSpeciesName(taxId)
                shortTaxName = shortenTaxName(taxName)
                print(taxName)

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
                meanGC = data_helpers.getGenomicGCContent(taxId)  # this is actually the genomic GC% (not CDS only)

                # Fetch temperature data for this species (if available)
                optimalTemperatureData = data_helpers.getSpeciesProperty( taxId, 'optimum-temperature')
                optimalTemperature = None
                if not optimalTemperatureData[0] is None:
                    optimalTemperature = float(optimalTemperatureData[0])

                temperatureRangeData = data_helpers.getSpeciesProperty( taxId, 'temperature-range')
                temperatureRange = None
                if not temperatureRangeData[0] is None:
                    temperatureRange = temperatureRangeData[0]
                else:
                    temperatureRange = "Unknown"

                pairedFractionData = data_helpers.getSpeciesProperty( taxId, 'paired-mRNA-fraction')
                pairedFraction = None
                if not pairedFractionData[0] is None:
                    pairedFraction = float(pairedFractionData[0])

                    
                genomeSizeData = data_helpers.getSpeciesProperty( taxId, 'genome-size-mb')
                genomeSize = None
                if not genomeSizeData[0] is None:
                    genomeSize = float(genomeSizeData[0])

                proteinCountData = data_helpers.getSpeciesProperty( taxId, 'protein-count')
                proteinCount = None
                if not proteinCountData[0] is None:
                    proteinCount = int(proteinCountData[0])

                geneDensity = None
                if( (not genomeSize is None) and (not proteinCount is None)  ):
                    geneDensity = float(proteinCount)/genomeSize
                print(geneDensity)

                    
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
                #print(dirpval)
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

    return (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics)


def PCA(biasProfiles, xdata):
    X = np.vstack(biasProfiles.values())
    assert(len(X)==len(xdata))

    pca = decomposition.PCA()
    pca.fit(X)
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

    ax.set_aspect("equal")

    plt.grid(True)
    plt.savefig("pca_test.pdf")
    plt.savefig("pca_test.svg")
    plt.close(fig)
    
    a = list(zip(biasProfiles.keys(), X_reduced[:,0]))
    a.sort(key=lambda x:x[1])
    print("top:")
    print(a[:3])
    print("bottom:")
    print(a[-3:])
    return a


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
    files = sys.argv[1:]
    assert(len(files)>0)

    (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)

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
    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'mean_delta_lfe',      'short_tax_name', 'tax_group', 'gc_vs_dLFE',          title )
    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'optimal_temperature', 'short_tax_name', 'tax_group', 'dLFE_vs_temperature', title )

    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'paired_fraction',     'short_tax_name', 'tax_group', 'dLFE_vs_paired_fraction', title )
    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'paired_fraction',     'short_tax_name', 'tax_group', 'gc_vs_paired_fraction',   title )


    plotSpeciesVars( summaryStatistics, 'genomic_gc',     'gene_density',        'short_tax_name', 'tax_group', 'gc_vs_gene_density',   title )
    plotSpeciesVars( summaryStatistics, 'mean_delta_lfe', 'gene_density',        'short_tax_name', 'tax_group', 'dLFE_vs_gene_density', title )
    plotSpeciesVars( summaryStatistics, 'gene_density',   'paired_fraction',     'short_tax_name', 'tax_group', 'gene_density_vs_paired_fraction',   title )
    

    order = PCA(biasProfiles, xdata)

    orderMap = dict(order)
    yrange = heatmaplotProfiles(biasProfiles, 'MFEbias', dfProfileCorrs, [], None, lambda x: orderMap[x] )

    #tileFilename = getProfileHeatmapTile(511145, biasProfiles, dfProfileCorrs, yrange)
    #print(tileFilename)

    #order2 = PCAforParams(params)

    print("---------")
    params = {}
    for k,v in biasProfiles.items():
        params[k] = estimateProfileParams(v.index, v.values)

    return 0



if __name__=="__main__":
    from sys import exit
    exit( standalone() )
