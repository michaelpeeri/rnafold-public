from __future__ import print_function
import sys
from random import choice
from math import log10
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress, wilcoxon
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
import data_helpers
import species_selection_data


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
        
        
        



def plotBars(data, varnames):
    fig, ax1 = plt.subplots()
    print(data['tax_name'])

    #for varname in varnames:
    #    data.plot(y=varname, x='short_tax_name', ax=ax1, kind='bar')

    data = data.sort_values(by=varnames[0], ascending=False)

    #colors, _ = assignColors(data['tax_group'])
    #data.plot(y=varnames, x='short_tax_name', color=colors, ax=ax1, kind='bar')
    data.plot(y=varnames, x='short_tax_name', ax=ax1, kind='bar')

    #plt.xlabel('GC-content')
    #ax1.set_ylabel('Native LMFE')
    #ax2.set_ylabel('Shuffled LMFE')
    #plt.grid(True)
    ax1.legend(loc=(0,1), scatterpoints=1, ncol=3, fontsize='small')
    plt.savefig("bars_%s.pdf" % '_'.join(varnames))
    plt.savefig("bars_%s.svg" % '_'.join(varnames))
    plt.close(fig)
    
    
    

xdata = []
ydata = []
ydata_nativeonly = []
ydata_shuffledonly = []
labels = []
groups = []
filesUsed = 0

short_names = set()
# TODO: This t
def shortenTaxName(name):
    currLength=4
    while(currLength <= len(name)):
        candidate = name[:currLength]
        if not candidate in short_names:
            short_names.add(candidate)
            return candidate
        currLength += 1
    raise Exception("Failed to shorten name '%s'" % name)

#SMFEcorrelations = pd.DataFrame({'spearman:smfe:gc':pd.Series(dtype='float'), 'spearman:smfe:cds_length_nt':pd.Series(dtype='float'), 'tax_name':pd.Series(dtype='string'), 'tax_id':pd.Series(dtype='int'), 'genomic_gc':pd.Series(dtype='float'), 'tax_group':pd.Series(dtype='category')})
SMFEcorrelations = pd.DataFrame({
    'spearman:smfe:gc':pd.Series(dtype='float'),
    'spearman:smfe:cds_length_nt':pd.Series(dtype='float'),
    'tax_name':pd.Series(dtype='string'),
    'short_tax_name':pd.Series(dtype='string'),
    'tax_id':pd.Series(dtype='int'),
    'genomic_gc':pd.Series(dtype='float'),
    'tax_group':pd.Series(dtype='string') # TODO: change to categorical data; Categorical([], categories=('Bacteria', 'Archaea', 'Fungi', 'Plants'), ordered=False)
})


for h5 in files:
    with pd.io.pytables.HDFStore(h5) as store:
        for key in store.keys():
            if not key.startswith( "/spearman_rho_" ):
                continue

            dfHeader = key.split('_')
            taxId = int(dfHeader[2])
            taxName = getTaxName(taxId)
            shortTaxName = shortenTaxName(taxName)

            


            print(taxName)

            dfSpearman = store[key]

            #df = store["/df_"+key[14:]]
            #df = df.iloc[:-1]  # remove the last value (which is missing)
            #meanGC = np.mean(df.gc)
            meanGC = species_selection_data.findByTaxid(taxId).iloc[0]['GC% (genome)']
            assert(meanGC > 10.0 and meanGC < 90.0)

            #dfSpearman = df.iloc[:-1]

            #print(dfSpearman)
            #print(dfSpearman['logpval']['gc'])

            SMFE_corr_gc = None
            SMFE_corr_cds_length = None
            try:
                SMFE_corr_gc = dfSpearman['logpval']['gc']
            except Exception:
                pass
            
            try:
                SMFE_corr_cds_length = dfSpearman['logpval']['cds_length_nt']
            except Exception:
                pass

            taxGroup = data_helpers.getSpeciesTaxonomicGroup(taxId)
            print(taxGroup)

            
            SMFEcorrelations = SMFEcorrelations.append(pd.DataFrame({
                'spearman:smfe:gc':pd.Series([SMFE_corr_gc]),
                'spearman:smfe:cds_length_nt': pd.Series([SMFE_corr_cds_length]),
                'tax_name':pd.Series([taxName]),
                'short_tax_name':pd.Series([taxName[:4]]),
                'tax_id':pd.Series([taxId]),
                'genomic_gc':pd.Series([meanGC]),
                'tax_group':pd.Series([taxGroup])
            }))

            # Format:

            #         gc  native  position  shuffled
            # 1    0.451  -4.944         1    -5.886
            # 2    0.459  -5.137         2    -6.069
            # 3    0.473  -5.349         3    -6.262
            filesUsed += 1

            #print(df.shape)

            #meanGC = np.mean(df.gc)
            #xdata.append(meanGC)

            #meanE = np.mean(df.native - df.shuffled)
            #ydata.append(meanE)

            #dirpval = calcWilcoxonPvalue_method2(deltas_df)
            #print(dirpval)
            #ydata.append(dirpval)


            #print(df.native)
            #print(df.shuffled)

            #meanE_nativeonly = np.mean(df.native)
            #ydata_nativeonly.append(meanE_nativeonly)

            #meanE_shuffledonly = np.mean(df.shuffled)
            #ydata_shuffledonly.append(meanE_shuffledonly)
            

#plotXY(xdata, ydata, labels, groups)
#plotXY_2(xdata, ydata_nativeonly, labels, groups)
#plotXY_3(xdata, ydata_nativeonly, ydata_shuffledonly, labels, groups)

print(SMFEcorrelations)

plotBars(SMFEcorrelations, ['spearman:smfe:gc', 'spearman:smfe:cds_length_nt'])

print("%d files included" % filesUsed)

