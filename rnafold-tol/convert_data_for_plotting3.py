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
import json
from datetime import datetime
from collections import Counter
#import redis
import numpy as np
import pandas as pd
import config
import mysql_rnafold as db
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
from data_helpers import CDSHelper, countSpeciesCDS, getSpeciesName, getSpeciesFileName, SpeciesCDSSource, numItemsInQueue, RateLimit, getAllComputedSeqs
from runningstats import RunningStats, OfflineStats



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
numShuffledGroups = 25
#computationTag = "rna-fold-window-40-0"
# TODO: Add support for step-size >1

# Establish DB connections
#r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
#session = db.Session()


rl = RateLimit(60)
rl2 = RateLimit(300)

computed = getAllComputedSeqs(seriesSourceNumber)
print("Collecting data from %d computation results..." % len(computed))


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


store = pd.io.pytables.HDFStore("gcdata.h5")

def printOutput(taxId, nativeProfile, shuffleProfiles, GCProfile):
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
    store.flush()


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

        out.append(float(gcCount) / totalCount)

    assert(len(out)==numWindows)
    return out




for taxIdForProcessing in species:
    print("Procesing %d sequences for tax-id %d (%s)..."
          % (countSpeciesCDS(taxIdForProcessing),
             taxIdForProcessing,
             getSpeciesName(taxIdForProcessing)))

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

    # Iterate over all CDS entries for this species
    for protId in SpeciesCDSSource(taxIdForProcessing):
        #protId = codecs.decode(protId)
        # Filtering


        # Skip sequences with partial CDS annotations
        #if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        #    skipped += 1
        #    continue

        #if( not r.exists(nativeCdsSeqIdKey % (taxIdForProcessing, protId)) ):
        #    skipped +=1
        #    continue

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

        for n in range(numShuffledGroups):
            #match[shuffledIds[n]] = (protId,n)
            if shuffledIds[n] in computed:
                computedShufflesCount += 1

        if( computedShufflesCount<numShuffledGroups-5 or (not cdsSeqId in computed) ):
            skipped += 1
            continue

        # TODO - do something here...
        #cdsResults = cds.getCalculationResult( seriesSourceNumber, -1 )
        #print(cdsResults[:20])

        cdsSequence = cds.sequence()
        gcContent = calcGCcontent(cdsSequence)
        #print(gcContent)
        for i in range(numWindows):  # Limit to 150; TODO - treat this generically?
            GCProfile[i].push( gcContent[i] )
        #print(GCProfile[i].mean())

        results = cds.getCalculationResult2( seriesSourceNumber, range(-1,numShuffledGroups) )
        if( results is None or len(results) < numShuffledGroups - 5):
            print("Not enough results found for %s" % protId)
            skipped += 1
            continue

        for shuffleId, content in zip(range(-1,numShuffledGroups), results):
            if( content is None ):
                print("Warning: Missing data for protein %s, shuffle-id %d" % (protId, shuffleId))
                continue
            data = json.loads(content.replace('id=', '"id":').replace('seq-crc=', '"seq-crc":').replace('MFE-profile=','"MFE-profile":').replace('MeanMFE=','"Mean-MFE":'))
            # Make sure we are seeing the correct record
            recordIdentifier = data["id"].split(":")
            assert(int(recordIdentifier[0]) == taxIdForProcessing)
            assert(recordIdentifier[1] == protId )
            assert(int(recordIdentifier[3]) == shuffleId )

            profile = data["MFE-profile"]

            if(shuffleId<0):
                for i in range(numWindows):
                    assert(profile[i] <= 0.0)
                    nativeProfile[i].push( profile[i] )
            else:
                for i in range(numWindows):
                    assert(profile[i] <= 0.0)
                    shuffleProfiles[shuffleId][i].push( profile[i] )
                

        if(rl()):
            print("# %s - %d records included, %d records skipped" % (datetime.now().isoformat(), selected, skipped))
            if( nativeProfile[0].count() > 1005 and rl2()):
                printOutput(taxIdForProcessing, nativeProfile, shuffleProfiles, GCProfile )
        
        selected += 1
        alreadyCompleted += 1

    #print(len(match))


    print("#%d selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
    #print("queue contains %d items" % numItemsInQueue(computationTag))

    print("# Counts: native %d, shuffled %d" % (nativeProfile[0].count(), shuffleProfiles[0][0].count()))
    printOutput(taxIdForProcessing, nativeProfile, shuffleProfiles, GCProfile)
