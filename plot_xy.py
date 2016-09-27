from __future__ import print_function
import sys
from random import choice
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
import data_helpers


files = sys.argv[1:]

def getTaxName(taxId):
    return data_helpers.getSpeciesFileName(taxId)

def assignColors(data):
    assignment = {}
    availColors = ('red', 'blue', 'green', 'yellow', 'cyan')
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
            plt.scatter(xdata, ydata, c=color, label=label, s=80)

    for x,y,l in zip(xvals, yvals, _labels):
        plt.annotate(l, (x-0.006, y+0.02))


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


    plt.xlabel('GC% near CDS start')
    plt.ylabel('Delta LMFE (Native - Shuffled)')
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
            plt.scatter(xdata, ydata, c=color, label=label, s=80)

    for x,y,l in zip(xvals, yvals, _labels):
        plt.annotate(l, (x-0.006, y+0.02))

    plt.xlabel('GC% near CDS start')
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
            ax1.scatter(xdata, ydata1, c=color, label=label, s=80)
            ax2.scatter(xdata, ydata2, c=color, label=label, s=80)
            #plt.scatter(xdata, ydata, c=color, label=label, s=80)

    for x,y1,y2,l in zip(xvals, yvals1, yvals2, _labels):
        #plt.annotate(l, (x-0.006, y+0.02))
        ax1.annotate(l, (x-0.006, y1+0.02))
        ax2.annotate(l, (x-0.006, y2+0.02))

    plt.xlabel('GC% near CDS start')
    ax1.set_ylabel('Native LMFE')
    ax2.set_ylabel('Shuffled LMFE')
    plt.grid(True)
    ax1.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
    plt.savefig("scatter_xy_split.pdf")
    plt.savefig("scatter_xy_split.svg")
    plt.close(fig)


for h5 in files:
    xdata = []
    ydata = []
    ydata_nativeonly = []
    ydata_shuffledonly = []
    labels = []
    groups = []
    with pd.io.pytables.HDFStore(h5) as store:
        for key in store.keys():
            if key[:4] != "/df_":
                continue

            taxId = int(key[4:])
            taxName = getTaxName(taxId)

            df = store[key]

            print(taxName)

            meanGC = np.mean(df.gc)
            xdata.append(meanGC)

            meanE = np.mean(df.native - df.shuffled)
            ydata.append(meanE)

            meanE_nativeonly = np.mean(df.native)
            ydata_nativeonly.append(meanE_nativeonly)

            meanE_shuffledonly = np.mean(df.shuffled)
            ydata_shuffledonly.append(meanE_shuffledonly)
            
            labels.append( taxName )

            #groups.append( choice(('Bacteria', 'Archaea', 'Fungi', 'Plants')) )
            groups.append( data_helpers.getSpeciesTaxonomicGroup(taxId) )

    plotXY(xdata, ydata, labels, groups)
    plotXY_2(xdata, ydata_shuffledonly, labels, groups)
    plotXY_3(xdata, ydata_nativeonly, ydata_shuffledonly, labels, groups)

            

                

