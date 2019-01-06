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
#import redis
import config
import mysql_rnafold as db
#from sqlalchemy import sql
#from sqlalchemy.sql.expression import func
from data_helpers import CDSHelper, countSpeciesCDS, getSpeciesName, SpeciesCDSSource, numItemsInQueue, RateLimit, getAllComputedSeqs
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
numShuffledGroups = 50
#computationTag = "rna-fold-window-40-0"
# TODO: Add support for step-size >1

# Establish DB connections
#r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
#session = db.Session()

skipped = 0
selected = 0
alreadyCompleted = 0


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

# build an average profile for the native proteins
nativeProfile = buildProfile(numWindows)

# build a average profiles for each of the shuffled groups
shuffleProfiles = []
for i in range(numShuffledGroups):
    shuffleProfiles.append( buildProfile(numWindows) )


def printOutput():
    nativeMean =    ["Native_mean"]
    nativeStdev =   ["Native_stdev"]
    nativeMin =     ["Native_min"]
    nativeMax =     ["Native_max"]

    shuffledMean =  ["Shuffled_mean"]
    shuffledStdev = ["Shuffled_stdev"]
    shuffledMin =   ["Shuffled_min"]
    shuffledMax =   ["Shuffled_max"]

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
                
    for d in (nativeMean,nativeStdev,nativeMin,nativeMax,shuffledMean,shuffledStdev,shuffledMin,shuffledMax):
        print(",".join(d))


for taxIdForProcessing in species:
    print("Procesing %d sequences for tax-id %d (%s)..."
          % (countSpeciesCDS(taxIdForProcessing),
             taxIdForProcessing,
             getSpeciesName(taxIdForProcessing)))

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

        if(rl()):
            print("# %s - %d records included, %d records skipped" % (datetime.now().isoformat(), selected, skipped))
            if( nativeProfile[0].count() > 1005 and rl2()):
                printOutput()

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

        if( computedShufflesCount<numShuffledGroups or (not cdsSeqId in computed) ):
            skipped += 1
            continue

        # TODO - do something here...
        #cdsResults = cds.getCalculationResult( seriesSourceNumber, -1 )
        #print(cdsResults[:20])

        results = cds.getCalculationResult2( seriesSourceNumber, range(-1,numShuffledGroups) )

        for shuffleId, content in zip(range(-1,numShuffledGroups), results):
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
                
        
        selected += 1
        alreadyCompleted += 1

#print(len(match))


print("#%d selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
#print("queue contains %d items" % numItemsInQueue(computationTag))

print("# Counts: native %d, shuffled %d" % (nativeProfile[0].count(), shuffleProfiles[0][0].count()))
printOutput()
