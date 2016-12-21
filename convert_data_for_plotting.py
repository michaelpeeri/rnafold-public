# Enqueue failed entries (i.e. entries missing the calculation results) for repeated processing.
# Scan all CDS entries from a given taxid; Find entries missing a computation result; Re-insert them into the queue for repeated processing
# Input - taxid1,taxid2,taxid3
#       - computationTag
#       - randomFraction
# Example:
# python2 requeue_sequences_missing_energies_for_sliding_window.py 3055,556484 rna-fold-window-40-0 10
# 
# TODO: This script will requeue sequences that have already in the queue but haven't been completed yet.
# TODO: Add support for step-size >1
import sys
import codecs
from random import randint
from math import floor
import json
from datetime import datetime
from collections import Counter
from bz2 import BZ2File
import numpy as np
import pandas as pd
from scipy import stats
import config
import mysql_rnafold as db
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
from data_helpers import CDSHelper, countSpeciesCDS, getSpeciesName, getSpeciesFileName, SpeciesCDSSource, numItemsInQueue, getAllComputedSeqsForSpecies
from runningstats import RunningStats, OfflineStats
from rate_limit import RateLimit



# command-line arguments
species = map(int, sys.argv[1].split(","))
#computationTag = sys.argv[2]
#if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0
#randomFraction = int(sys.argv[3])

# Configuration
#queueKey = "queue:tag:awaiting-%s:members" % computationTag
#nativeCdsSeqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
#seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
windowWidth = 40
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
numWindows = 150
numShuffledGroups = 20
requiredNumShuffledGroups = 20
#computationTag = "rna-fold-window-40-0"
# TODO: Add support for step-size >1

# Establish DB connections
#r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
#session = db.Session()


rl = RateLimit(60)
rl2 = RateLimit(300)


#match = {}

def buildProfile(length, type="online"):
    out = []
    for i in range(length):
        if( type=="online"):
            out.append(RunningStats())
        elif( type=="offline"):
            out.append(OfflineStats())
        else:
            assert(False)
    return out

def addMargins(x0, x1, margin=0.1):
    if( x0 > x1 ):
        (x0,x1) = (x1,x0)
    assert(x1 >= x0)

    _margin = abs(x1-x0) * 0.1
    y0, y1 = 0, 0
    y0 = x0 - _margin

    y1 = x1 + _margin

    return (y0,y1)
    
        


# xaxis = np.arange(0,150,1)
# def plotProfile_(taxId, native, shuffle):
#     fig = plt.figure()
#     print(xaxis.shape, native.shape, shuffle.shape)

#     plt.plot(xaxis, native, label="Native")
#     plt.plot(xaxis, shuffle, label="Shuffled (mean)")
#     plt.xlabel('Position (nt, window start, from cds)')
#     plt.ylabel('Mean LFE')

#     yrange = addMargins(min(np.amin(native), np.amin(shuffle)), max(np.amax(native), np.amax(shuffle)))
#     plt.axis([0, 150, yrange[0], yrange[1]])

#     plt.title(getSpeciesName(taxId))
#     plt.legend()
#     plt.grid(True)
#     plt.savefig("mfe_40nt_cds_%s.pdf" % getSpeciesFileName(taxId) )
#     plt.close(fig)


plt.style.use('ggplot') # Use the ggplot style

def plotProfile(taxId, data):
    fig, (ax1,ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    data[['native', 'shuffled']].plot(ax=ax1)

    #plt.title(getSpeciesName(taxId))

    plt.xlabel('Position (nt, window start, from cds)')

    ax1.set_title("Mean LFE for %s" % getSpeciesName(taxId))
    ax1.set_ylabel('Mean LFE')
    ax1.legend()
    ax1.grid(True)


    data['gc'].plot(ax=ax2)
    ax2.set_title("GC%")
    ax2.set_ylabel('GC% (in window)')
    ax2.grid(True)


    plt.savefig("mfe_40nt_cds_%s.pdf" % getSpeciesFileName(taxId) )
    plt.savefig("mfe_40nt_cds_%s.svg" % getSpeciesFileName(taxId) )
    plt.close(fig)


def plotXY(xvals, yvals, _labels):
    fig = plt.figure()
    plt.scatter(xvals, yvals )

    for x,y,l in zip(xvals, yvals, _labels):
        plt.annotate(l, (x-0.1, y+0.2))

    plt.xlabel('window start (nt, from cds)')
    plt.ylabel('MFE')
    plt.grid(True)
    plt.savefig("scatter_xy.pdf")
    plt.close(fig)


def plotGCContentHist( medianGCContent, taxId ):
    histBins = np.arange(0.24, 0.76, 0.02)
    hist, _ = np.histogram( medianGCContent, bins = histBins )
    print(len(medianGCContent))
    print(hist)
    print(histBins)

    fig = plt.figure()
    #plt.hist( hist, histBins, normed=0 )

    plt.bar( np.arange(0.24, 0.74, 0.02), hist, width=0.02)
    plt.xlabel('Median GC% (5\' cds section)')
    plt.ylabel('CDS count')
    plt.xticks(np.arange(0.24, 0.76, 0.04))

    plt.savefig("mediangc_%s.pdf" % getSpeciesFileName(taxId) )
    plt.savefig("mediangc_%s.svg" % getSpeciesFileName(taxId) )
    plt.close(fig)
    



def printOutput(taxId, nativeProfile, shuffleProfiles, GCProfile, medianGCContent):
    store = pd.io.pytables.HDFStore("gcdata_taxid_%d.h5" % taxId)
    nativeMean =    ["Native_mean"]
    nativeStdev =   ["Native_stdev"]
    nativeMin =     ["Native_min"]
    nativeMax =     ["Native_max"]

    shuffledMean =  ["Shuffled_mean"]
    shuffledStdev = ["Shuffled_stdev"]
    shuffledMin =   ["Shuffled_min"]
    shuffledMax =   ["Shuffled_max"]

    gcMean =        ["GCcontent_mean"]

    #build a "profile-of-profiles" for the shuffled groups
    aggregateProfile = buildProfile(numWindows, "offline")
    for n in range(numShuffledGroups):
        for i in range(numWindows):
            if( shuffleProfiles[n][i].count() ):
                aggregateProfile[i].push( shuffleProfiles[n][i].mean() )

    if( not nativeProfile[i].count() ):
        return

    for i in range(numWindows):
        nativeMean.append(    "%.4g" % nativeProfile[i].mean() )
        nativeStdev.append(   "%.4g" % nativeProfile[i].stdev() )
        nativeMin.append(     "%.3g" % nativeProfile[i].min() )
        nativeMax.append(     "%.3g" % nativeProfile[i].max() )

        shuffledMean.append(  "%.4g" % aggregateProfile[i].mean() )
        shuffledStdev.append( "%.4g" % aggregateProfile[i].stdev() )
        shuffledMin.append(   "%.3g" % aggregateProfile[i].min() )
        shuffledMax.append(   "%.3g" % aggregateProfile[i].max() )

        gcMean.append( "%.3g" % GCProfile[i].mean() )
                
    for d in (nativeMean,nativeStdev,nativeMin,nativeMax,shuffledMean,shuffledStdev,shuffledMin,shuffledMax, gcMean):
        print(",".join(d))

    #plotProfile(taxId, np.asarray(nativeMean[1:], dtype="float"), np.asarray(shuffledMean[1:], dtype="float"))

    # TODO - optimize this to save multiple reverse conversions!

    native = np.asarray(nativeMean[1:], dtype="float")
    shuffled = np.asarray(shuffledMean[1:], dtype="float")
    gc = np.asarray(gcMean[1:], dtype="float")
    xrange = np.arange(1,151,1)
    df = pd.DataFrame( { "native": native, "shuffled": shuffled, "gc": gc, "position": xrange}, index=xrange )
    plotProfile(taxId, df)
    store["df_%d" % taxId] = df


    plotGCContentHist( medianGCContent, taxId )
    medianGCContentHist, binEdges = np.histogram( medianGCContent, bins = np.arange(24, 76, 2) )
    
    store.flush()

#refGamma = stats.gamma(a=13.5) # Test reference distribution

def calcGCcontent(cdsSequence):
    #windowWidth = 40
    #seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
    #numWindows = 150
    if( len(cdsSequence) < numWindows+windowWidth-1 ):
        return

    out = []

    # DEBUG ONLY ##### DEBUG ONLY ##### DEBUG ONLY ##### DEBUG ONLY ##### DEBUG ONLY #
    #return out
    # DEBUG ONLY ##### DEBUG ONLY ##### DEBUG ONLY ##### DEBUG ONLY ##### DEBUG ONLY #

    for start in range(numWindows):
        end = start+windowWidth
        fragment = cdsSequence[start:end]
        assert(len(fragment)==windowWidth)

        freqs = Counter(fragment.lower())
        # Calc GC content; Don't count ambiguous symbols
        # Note: It is also possible to count N's as 1/4 of each...
        totalCount = freqs['a'] + freqs['c'] + freqs['g'] + freqs['t']
        gcCount = freqs['c'] + freqs['g']

        if( totalCount > 0 ):
            out.append(float(gcCount) / totalCount)
        else:
            print("Warning: no data for window %d" % start)
            return []

    assert(len(out)==numWindows)
    return out


for taxIdForProcessing in species:
    print("Procesing %d sequences for tax-id %d (%s)..."
          % (countSpeciesCDS(taxIdForProcessing),
             taxIdForProcessing,
             getSpeciesName(taxIdForProcessing)))

    computed = getAllComputedSeqsForSpecies(seriesSourceNumber, taxIdForProcessing)
    print("Collecting data from %d computation results..." % len(computed))



    #fskew = open('skew_taxid_%d.csv' % taxIdForProcessing, 'w')
    #fskewpval = open('skewpval_taxid_%d.csv' % taxIdForProcessing, 'w')
    #fnorm = open('norm_taxid_%d.csv' % taxIdForProcessing, 'w')
    #fnormpval = open('normpval_taxid_%d.csv' % taxIdForProcessing, 'w')
    #fshapiro = open('shapiro_taxid_%d.csv' % taxIdForProcessing, 'w')
    #fshapiropval = open('shapiropval_taxid_%d.csv' % taxIdForProcessing, 'w')
    ##fanderson = open('anderson_taxid_%d.csv', % taxIdForProcessing, 'w')
    ##fandersonpval = open('andersonpval_taxid_%d.csv', % taxIdForProcessing, 'w')
    #fkurtosis = open('kurtosis_taxid_%d.csv' % taxIdForProcessing, 'w')
    #frefgamma = open('refgamma_taxid_%d.csv' % taxIdForProcessing, 'w')
    #frawdata = BZ2File('rawdata_taxid_%d.csv.bz2' % taxIdForProcessing, 'w', 2048)
    fdebug_zeros_scores = BZ2File('zeros_taxid_%d.csv.bz2' % taxIdForProcessing, 'w', 2048)
    fdebug_zeros_fasta = BZ2File('zeros_taxid_%d.fna.bz2' % taxIdForProcessing, 'w', 2048)
    
    fdebug_sample_scores = BZ2File('sample_taxid_%d.csv.bz2' % taxIdForProcessing, 'w', 2048)
    fdebug_sample_fasta = BZ2File('sample_taxid_%d.fna.bz2' % taxIdForProcessing, 'w', 2048)

    fwilcoxon_dist = BZ2File('wilcoxon_taxid_%d.csv.bz2' % taxIdForProcessing, 'w', 2048)

    
    skipped = 0
    selected = 0
    alreadyCompleted = 0

    # build an average profile for the native proteins
    nativeProfile = buildProfile(numWindows)

    # build a average profiles for each of the shuffled groups
    shuffleProfiles = []
    for i in range(numShuffledGroups):
        shuffleProfiles.append( buildProfile(numWindows) )

    # build an average profile for the GC content
    GCProfile = buildProfile(numWindows)

    medianGCContent = []

    # Iterate over all CDS entries for this species
    for protId in SpeciesCDSSource(taxIdForProcessing):
        cds = CDSHelper(taxIdForProcessing, protId)

        seqLength = cds.length()
        if( not seqLength is None ):
            # Skip sequences that are too short
            if(seqLength < numWindows + windowWidth + 1 ):
                skipped += 1
                continue
        else:
            print("Warning: Could not find CDS length entry for taxid=%d, protid=%s" % (taxIdForProcessing, protId) )
            skipped += 1
            continue

        requiredNumWindows = seqLength - windowWidth + 1

        cdsSeqId = cds.seqId()

        #match[cdsSeqId] = (protId,-1)

        shuffledIds = cds.shuffledSeqIds()

        computedShufflesCount = 0

        for n in range(min(numShuffledGroups, len(shuffledIds))):
            #match[shuffledIds[n]] = (protId,n)
            if shuffledIds[n] in computed:
                computedShufflesCount += 1

        if( computedShufflesCount<requiredNumShuffledGroups or (not cdsSeqId in computed) ):
            #print("%s - found only %d groups, skipping" % (protId, computedShufflesCount))
            skipped += 1
            continue

        # TODO - do something here...
        #cdsResults = cds.getCalculationResult( seriesSourceNumber, -1 )
        #print(cdsResults[:20])

        cdsSequence = cds.sequence()
        gcContent = calcGCcontent(cdsSequence)
        #print(gcContent)
        if( len(gcContent) ):
            for i in range(numWindows):  # Limit to 150; TODO - treat this generically?
                GCProfile[i].push( gcContent[i] )
                medianGCContent.append( np.median( gcContent ) )

        #print(GCProfile[i].mean())

        results = cds.getCalculationResult2( seriesSourceNumber, range(-1,numShuffledGroups), True )
        del cds

        if( results is None or len(filter(lambda x: not x is None, results)) < requiredNumShuffledGroups ):
            print("Not enough results found for %s" % protId)
            skipped += 1
            continue

        shuffles = pd.Index(range(len(results)-1))
        df = pd.DataFrame(index=shuffles, columns=range(151))  # TODO - set the profile lengths according to the data...

        for shuffleId, data in zip(range(-1,numShuffledGroups), results):
            if( data is None ):
                #print("Warning: Missing data for protein %s, shuffle-id %d" % (protId, shuffleId))
                continue

            # Make sure we are seeing the correct record
            recordIdentifier = None
            if( data["id"].find('/') == -1 ):
                recordIdentifier = data["id"].split(":")
            else:
                recordIdentifier = data["id"].split("/")
            
            assert(int(recordIdentifier[0]) == taxIdForProcessing)
            assert(recordIdentifier[1] == protId )
            assert(int(recordIdentifier[3]) == shuffleId )

            profile = data["MFE-profile"]  # MFE profile for all calculated positions in protId, shuffleId
            del data

            if( shuffleId>=0 ):
                df.iloc[shuffleId,:] = profile

            # Check all values are non-positive
            # Note: The energy values reaches 0.0 at some points (although I would've thought it should always be negative), so I ignore such cases.
            numNonNegativeResults = len(filter(lambda x: x >= 0.0, profile))
            if( numNonNegativeResults > 0 and shuffleId == -1):  # Ignore shuffled seqs for now (concentrate on native seqs)
                #print("Warning: %d stored values are not negative for taxId=%d, protId=%s, shuffleId=%d" % (numNonNegativeResults, taxIdForProcessing, protId, shuffleId))
                if( numNonNegativeResults > 15 ):
                    #print(profile)
                    #if(shuffleId >= 0):
                    #    print(cds.getShuffledSeq(shuffleId))
                    #else:
                    #    print(cds.sequence())
                    
                    # TODO - consider changing to use proper csv and fasta output...
                    fdebug_zeros_scores.write(",".join( ["%s_%s" % (taxIdForProcessing, protId)] + list(map( lambda x: "%.3g"%x, profile ))) )
                    fdebug_zeros_scores.write("\n")
                    #
                    fdebug_zeros_fasta.write(">%s_%s\n" % (taxIdForProcessing, protId))
                    fdebug_zeros_fasta.write(cdsSequence)
                    fdebug_zeros_fasta.write("\n")
                    
                assert(len(filter(lambda x: x > 0.0, profile)) == 0) # Positive results aren't valid
                
            if(shuffleId<0):
                for i in range(numWindows):
                    assert(profile[i] <= 0.0)
                    nativeProfile[i].push( profile[i] )
            else:
                for i in range(numWindows):
                    assert(profile[i] <= 0.0)
                    shuffleProfiles[shuffleId][i].push( profile[i] )
            #del profile

        #skew = []
        #skewpval = []
        #norm = []
        #normpval = []
        #shapiro = []
        #shapiropval = []
        #kurtosis = []
        #refgamma = []

        
        #anderson = []
        #andersonpval = []
        
        for i in range(150):
            s = df[i]

            # window is ready
            
            if(s.var() == 0):
                print("Skipping %d" % i)
                #for a in (skew, skewpval, norm, normpval, shapiro, shapiropval, kurtosis):
                #    a.append(None)
                continue

            if( len(s) < requiredNumShuffledGroups ):
                print("Error: only %d elements found, skipping" % len(s))
                continue

            #fittedDist = stats.gamma.fit( -s, floc=0.0 )

            
            #skewtest = None
            #normaltest = None
            #shapirotest = None
            #kurtosisval = None
            #refgammaval = None
            #try:
            #    skewtest = stats.skewtest(s)
            #    normaltest = stats.normaltest(s)
            #    shapirotest = stats.shapiro(s)
            #    kurtosisval = stats.kurtosis(s)
            #    refgammaval = stats.kstest(-1.0*np.array(s, dtype=np.dtype('float32')), 
            #                               refGamma.cdf).pvalue
            #    #andersontest = stats.anderson(s)
            #except Exception as e:
            #    print(e)
            #    print("While processing: ")
            #    print(protId)
            #    print(",".join(map(lambda x: "%.4g" % x, s)))
            #    raise
            
            #skew.append(skewtest.statistic)
            #skewpval.append(skewtest.pvalue)
            #norm.append(normaltest.statistic)
            #normpval.append(normaltest.pvalue)
            #shapiro.append(shapirotest[0])
            #shapiropval.append(shapirotest[1])
            #kurtosis.append(kurtosisval)
            #refgamma.append(refgammaval)

            #frawdata.write("%s,%d,%s\n" % (protId, i, ','.join(map(lambda x: "%.4g"%x, s))))
            ##anderson.append(shapirotest.statistic)
            ##andersonpval.append(shapirotest.pvalue)


        #fskew.write("%s,%s\n" % (protId, ','.join(map(str, skew))))
        #fskewpval.write("%s,%s\n" % (protId, ','.join(map(str, skewpval))))
        #fnorm.write("%s,%s\n" % (protId, ','.join(map(str, norm))))
        #fnormpval.write("%s,%s\n" % (protId, ','.join(map(str, normpval))))
        #fshapiro.write("%s,%s\n" % (protId, ','.join(map(str, shapiro))))
        #fshapiropval.write("%s,%s\n" % (protId, ','.join(map(str, shapiropval))))
        #fkurtosis.write("%s,%s\n" % (protId, ','.join(map(str, kurtosis))))
        #frefgamma.write("%s,%s\n" % (protId, ','.join(map(str, refgamma))))

        # Perform Wilcoxon signed-rank test
        #print("DF: ", df.shape)  # (20,151)
        #print("s: ", s.shape)    # (20,)
        #print("u: ", df.iloc[1,:].shape) # (151,)
        #print("mean(u): ", df.mean(axis=0).shape)
        #print(np.array(profile).shape) #(151,)


        # Method 1 -- Wilcoxon, native vs. mean(shuffled), window step=1nt, per gene
        meanOfShufflesProfileForW = df.mean(axis=0)
        nativeProfileForW = np.array(profile)
        assert(meanOfShufflesProfileForW.shape == nativeProfileForW.shape)

        wilx1 = stats.wilcoxon(nativeProfileForW, meanOfShufflesProfileForW)
        #print(wilx)
        #fwilcoxon_dist.write("%s,%.4g,%.3g\n" % (protId, wilx1.statistic, wilx1.pvalue))


        # Method 2 -- Wilcoxon, native vs. shuffled, window step=30nt, per gene
        startPositions = [int(floor(i*29.9)) for i in range(6)] # 0-based; Use 29.9 to force the last window to 150 (vs. 151); TODO - Change this...
        assert(len(startPositions)==6)
        assert(np.all(np.array([(startPositions[i+1]-startPositions[i]) >= 29 for i in range(5)])))
        shuffleScores = df.iloc[:,startPositions]
        assert(shuffleScores.shape == (20,len(startPositions)))
        #
        assert(nativeProfileForW.shape == (151,))
        nativeProfiles = nativeProfileForW.take(startPositions)
        assert(nativeProfiles.shape == (len(startPositions),))

        x = nativeProfiles - shuffleScores
        assert(x.shape == (20,len(startPositions)))
        wilx2 = stats.wilcoxon( x.unstack().values )
        fwilcoxon_dist.write("%s,%.4g,%.3g,%.4g,%.3g,%.3g,%.3g,%.3g\n" % (protId, wilx1.statistic, wilx1.pvalue, wilx2.statistic, wilx2.pvalue, np.median(nativeProfiles), np.mean( gcContent ), np.median(x.unstack().values) ) )
        
        if(rl()):
            print("# %s - %d records included, %d records skipped" % (datetime.now().isoformat(), selected, skipped))
            if( nativeProfile[0].count() > 1005 and rl2()):
                printOutput(taxIdForProcessing, nativeProfile, shuffleProfiles, GCProfile, medianGCContent )
                
            #for f in (fskew, fskewpval, fnorm, fnormpval, fshapiro, fshapiropval, fkurtosis):
            #    f.flush()
        
        selected += 1
        alreadyCompleted += 1

    #print(len(match))


    print("#%d selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
    #print("queue contains %d items" % numItemsInQueue(computationTag))

    print("# Counts: native %d, shuffled %d" % (nativeProfile[0].count(), shuffleProfiles[0][0].count()))
    printOutput(taxIdForProcessing, nativeProfile, shuffleProfiles, GCProfile, medianGCContent)

    
    for f in (fdebug_zeros_scores, fdebug_zeros_fasta, fdebug_sample_scores, fdebug_sample_fasta, fwilcoxon_dist): # (fskew, fskewpval, fnorm, fnormpval, fshapiro, fshapiropval, fkurtosis, frawdata):
        f.close()
