from __future__ import print_function
import sys
from random import choice
import numpy as np
import pandas as pd
import matplotlib
from math import log10
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress, wilcoxon
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
import data_helpers


files = sys.argv[1:]

def getTaxName(taxId):
    return data_helpers.getSpeciesFileName(taxId)

def assignColors(data):
    assignment = {}
    availColors = ('red', 'blue', 'white', 'green', 'yellow', 'cyan', 'grey')
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
        plt.annotate(label, (x+0.003, y-0.003), fontsize=7)


    rmin = min(min(xvals), min(yvals))-0.1
    rmax = max(max(xvals), max(xvals))+0.1
    plt.plot([rmin,0.0,rmax], [0.0, 0.0, 0.0], color='black')
    plt.plot([0.0,0.0,0.0], [rmin, 0.0, rmax], color='black')
    plt.plot([rmin,0.0,rmax], [rmin, 0.0, rmax], '--', color='black')

    plt.annotate(s="n= %d"  % len(yvals),                                                 xy=(rmin+0.05, rmax-0.05 ),  fontsize=6 )
        

    plt.xlim([rmin,rmax])
    plt.ylim([rmin,rmax])
    plt.xlabel('Local MFE selection (native-shuffled), CDS pos < 150nt')
    plt.ylabel('Local MFE selection (native-shuffled), CDS pos >= 400nt')
    plt.grid(True)
    plt.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
    plt.savefig("mfe_cds_begin_vs_rest.pdf")
    plt.savefig("mfe_cds_begin_vs_rest.svg")
    plt.close(fig)

    
xdata = []
ydata = []
ydata_nativeonly = []
ydata_shuffledonly = []
labels = []
groups = []
filesUsed = 0

for h5 in files:
    with pd.io.pytables.HDFStore(h5) as store:
        for key in store.keys():
            if key[:4] != "/df_":
                continue

            dfHeader = key.split('_')
            taxId = int(dfHeader[1])
            taxName = getTaxName(taxId)
            print(taxName)

            df = store[key]
            df = df.iloc[:-1]
            # Format:

            #         gc  native  position  shuffled
            # 1    0.451  -4.944         1    -5.886
            # 2    0.459  -5.137         2    -6.069
            # 3    0.473  -5.349         3    -6.262
            filesUsed += 1

            #print(df.shape)

            dfBegin = df.iloc[:15]
            dfRest = df.iloc[40:]

            meanE_begin = np.mean(dfBegin.native - dfBegin.shuffled)
            xdata.append(meanE_begin)

            meanE_end = np.mean(dfRest.native - dfRest.shuffled)
            ydata.append(meanE_end)



            #meanE_nativeonly = np.mean(df.native)
            #ydata_nativeonly.append(meanE_nativeonly)

            #meanE_shuffledonly = np.mean(df.shuffled)
            #ydata_shuffledonly.append(meanE_shuffledonly)
            
            labels.append( taxName )

            #groups.append( choice(('Bacteria', 'Archaea', 'Fungi', 'Plants')) )
            groups.append( data_helpers.getSpeciesTaxonomicGroup(taxId) )

plotXY(xdata, ydata, labels, groups)
print("%d files included" % filesUsed)

