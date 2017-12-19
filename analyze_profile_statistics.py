from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, relfreq
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from mysql_rnafold import Sources
from process_series_data import convertResultsToMFEProfiles, readSeriesResultsForSpecies, sampleProfilesFixedIntervals
from data_helpers import allSpeciesSource



# ------------------------------------------------------------------------------------
# Configuration
confWindowWidth = 40


# ------------------------------------------------------------------------------------
# Command-line args

def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert
    
def parseProfileSpec():
    def convert(value):
        o = value.split(':')
        assert(len(o) >= 3 and len(o) <= 4)
        
        o[0] = int(o[0])
        assert(o[0]>0)
        
        o[1] = int(o[1])
        assert(o[1]>0)
        
        assert(o[2]=="begin" or o[2]=="end")

        if( len(o) == 4 ):
            o[3] = int(o[3])
        else:
            o.append(0)
        
        return (o[0], o[1], o[2], o[3])
    return convert


def plot2dProfile(profileData, taxid):
    fig, ax = plt.subplots(nrows=2, ncols=2, #sharey=True, #sharex=True,
                           gridspec_kw={'height_ratios': [1, 3], 'width_ratios': [3, 1]} )


    x = np.apply_along_axis( lambda x: np.histogram(x[~np.isnan(x)], bins=100, range=(-5,5), density=True), 0, profileData)
    #print(x.shape)
    freqs = np.vstack(x[0,])

    
    #ax.imshow( profileData, cmap='coolwarm', aspect='auto' )
    ax[1][0].imshow( freqs.T, cmap='gist_heat', aspect='auto' )


    rr = profileData.ravel()
    print(np.histogram(rr[~np.isnan(rr)], bins=10))
    hist, bins = np.histogram(rr[~np.isnan(rr)], bins=100, range=(-5,5), density=True)

    print(np.mean(rr))
    print(wilcoxon(rr))

    ax[1][1].plot( hist, bins[:-1] )
    ax[1][1].set_ylim((-5,5))
    
    plt.savefig("2dprofile_taxid_%d.png" % taxid, orientation='portrait')
    plt.savefig("2dprofile_taxid_%d.svg" % taxid, orientation='portrait')
    plt.close(fig)




def calculate2dProfile(args):

    maxLength = 300
    profileStep = 10

    taxids = []
    if args.all_species:
        taxids = [x for x in allSpeciesSource()]
    else:
        taxids = args.taxid

    for taxid in taxids:
        count = 0
        nativeArrayData  = []
        controlArrayData = []

        testt = dict(map(lambda x: (x,0), range(21)))
        
        for result in sampleProfilesFixedIntervals(
                convertResultsToMFEProfiles(
                    readSeriesResultsForSpecies( args.computation_tag, taxid, args.num_shuffles, args.num_shuffles )
                    , args.num_shuffles)
                , startPosition=0, endPosition=maxLength, interval=profileStep):
            profileData = result["profile-data"]

            # Check the sequence-id
            seqId = result["content"][0]["id"]
            if( seqId.find(":") != -1 ):
                seqId = seqId.replace(":", "/")
            shuffleId = int(seqId.split("/")[3])  # the first result should belong to shuffle-id -1 (i.e., the native sequence)
            assert(shuffleId==-1)

            expectedProfileLength = min( len( result["content"][0]["MFE-profile"]) / profileStep, maxLength / profileStep )
            profileLength = profileData.shape[1]
            assert( abs( profileLength - expectedProfileLength) <= 1 )
            
            if( profileData.shape[0] != args.num_shuffles + 1 ):  # we require one vector per suffled sequence, plus one for the native sequence
                print("Warning: ignoring record '%s' containing %d records" % ( seqId, profileData.shape[0]))
                continue

            nativeDiffs  = profileData[0, ]  - profileData[1:, ]   # Calculate NativeLFE - ShuffledLFE (for each of the shuffles, and for each window)
            #controlDiffs = profileData[-1,] - profileData[:-1,]   # Calculate NativeLFE - ShuffledLFE (for each of the shuffles, and for each window)
            controlDiffs = profileData[8, ] - profileData[1:, ]   # Calculate NativeLFE - ShuffledLFE (for each of the shuffles, and for each window)
            assert( nativeDiffs.shape[0]  == args.num_shuffles)
            assert( controlDiffs.shape[0] == args.num_shuffles)

            for i in range(21):
                deltas = profileData[i,] - profileData
                T, pval = wilcoxon(deltas.ravel())
                if( pval < 0.05 ):
                    testt[i] += 1
                    #print("%d - pval: %g" % (i, pval))

            direction        = np.sign( np.apply_along_axis( np.mean, 0, nativeDiffs) )  # TODO - prove this is equivalent to checking whether the sign of the sum of ranks
            controlDirection = np.sign( np.apply_along_axis( np.mean, 0, controlDiffs) )  # TODO - prove this is equivalent to checking whether the sign of the sum of ranks

            wilc = np.apply_along_axis( wilcoxon, 0, nativeDiffs )   # The wilcoxon test is performed separately for each window (with N=args.num_shuffles, typically = 20)

                                                               # Note that Nr <= N (because ties are ignored when using the default settings), and Nr should be at least 10 or 20 for the distribution to approach normal.
                                                               # See:
                                                               # Explanation of Python impl. (using T statistic):  https://stackoverflow.com/a/18966286
                                                               # Wilcoxon signed-rank test tutorial:               http://vassarstats.net/textbook/ch12a.html
            assert( wilc.shape == (2, profileLength) )
            #controlWilc = np.apply_along_axis( wilcoxon, 0, controlDiffs )   # The wilcoxon test is performed separately for each window (with N=args.num_shuffles, typically = 20)
            controlWilc = np.apply_along_axis( np.mean, 0, controlDiffs )   # The wilcoxon test is performed separately for each window (with N=args.num_shuffles, typically = 20)

            #print("--"*10)
            #print(direction)
            #print(np.mean(wilc[0,]))
            #Nr = sum( np.abs(nativeDiffs) > 1e-6 )
            #print(Nr)
            #sigma = np.sqrt(Nr * (Nr+1) * (2*Nr+1) / 6)   # =SD of W
            assert( np.all( wilc[0,]        >= 0.0 ))  # test statistic T (not W) - The sum of the ranks of the differences above or below zero, whichever is smaller
            #assert( np.all( controlWilc[0,] >= 0.0 ))  # test statistic T (not W) - The sum of the ranks of the differences above or below zero, whichever is smaller
            #S = Nr * (Nr+1) / 2.0  # Sum of all ranks
            #W = S - 2*wilc[0,]
            #Z = W / sigma
            #print(wilc[1,])
            assert( np.all( ((wilc[1,] >= 0.0) & (wilc[1,] <= 1.0)) | np.isnan(wilc[1,]) ) )   # P-values
            #print(Z * direction)

            #wilc.resize((2, maxLength / profileStep))  # pad with zeros
            #arrayData.append( Zwilc[1,] )
            #out = np.resize( Z*direction, (2, maxLength / profileStep))
            out        = np.resize( np.log10(wilc[1,])        * direction        * -1, (2, maxLength / profileStep))
            nativeArrayData.append( out )

            #controlOut = np.resize( np.log10(controlWilc[1,]) * controlDirection * -1, (2, maxLength / profileStep))
            controlOut = np.resize( controlWilc[1,], (2, maxLength / profileStep))
            controlArrayData.append( controlOut )

            count += 1

        if( not nativeArrayData ):
            print("Warning: no data found for taxid=%d" % taxid)
            continue
        
        nativeAr  = np.vstack( nativeArrayData )
        controlAr = np.vstack( controlArrayData )
        #print(ar.shape)
        #print(ar[0,])

        #x = np.apply_along_axis( lambda x: relfreq(x[~np.isnan(x)], numbins=100, defaultreallimits=(-5,5)), 0, ar)
        #x = np.apply_along_axis( lambda x: np.histogram(x[~np.isnan(x)], bins=100, range=(-5,5), density=True), 0, nativeAr)
        #print(x.shape)
        #nativeFreqs = np.vstack(x[0,])

        #y = np.apply_along_axis( lambda x: np.histogram(x[~np.isnan(x)], bins=100, range=(-5,5), density=True), 0, controlAr)
        #print(x.shape)
        #controlFreqs = np.vstack(y[0,])
        
        #print( np.apply_along_axis( np.sum, 1, controlFreqs ) )
        #assert( np.allclose( np.apply_along_axis( np.sum, 0, freqs ), 1.0 ) )
        #print(freqs.shape)
        #print(freqs[0])
        
        #plot2dProfile(nativeFreqs, taxid)
        print(testt)
        plot2dProfile(controlAr, taxid)

        print(count)
    return 0


def standalone():
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--taxid", type=parseList(int))
    argsParser.add_argument("--all-species", action="store_true", default=False)
    #argsParser.add_argument("--profile", type=parseProfileSpec())
    argsParser.add_argument("--computation-tag", type=int, default=Sources.RNAfoldEnergy_SlidingWindow40_v2)
    argsParser.add_argument("--num-shuffles", type=int, default=20)
    argsParser.add_argument("--pax-db", type=str, required=False)
    #argsParser.add_argument("--codonw", type=bool, default=False)
    args = argsParser.parse_args()
    return calculate2dProfile(args)
    

if __name__=="__main__":
    import sys
    sys.exit(standalone())
