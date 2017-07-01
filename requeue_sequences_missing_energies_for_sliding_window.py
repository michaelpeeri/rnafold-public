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
import argparse
from random import randint
#import redis
import config
import mysql_rnafold as db
#from sqlalchemy import sql
#from sqlalchemy.sql.expression import func
from data_helpers import CDSHelper, countSpeciesCDS, getSpeciesName, SpeciesCDSSource, numItemsInQueue
import _distributed
import dask
from store_new_shuffles import storeNewShuffles
from rnafold_sliding_window import calculateMissingWindowsForSequence

scheduler = _distributed.open()


def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert


argsParser = argparse.ArgumentParser()
argsParser.add_argument("--species", type=parseList(int), required=True)
argsParser.add_argument("--computation-tag", default="rna-fold-window-40-0")
argsParser.add_argument("--random-fraction", type=int, default=1)
argsParser.add_argument("--window-step", type=int, default=10)
argsParser.add_argument("--from-shuffle", type=int, default=-1)
argsParser.add_argument("--to-shuffle", type=int, default=20)
args = argsParser.parse_args()

# command-line arguments
species = args.species
computationTag = args.computation_tag
if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0
randomFraction = args.random_fraction
windowStep = args.window_step


# Configuration
#queueKey = "queue:tag:awaiting-%s:members" % computationTag
#nativeCdsSeqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
#seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
windowWidth = 40
lastWindowStart = 2000
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40_v2
expectedNumberOfShuffles = 20
fromShuffle = args.from_shuffle
toShuffle = args.to_shuffle
# TODO: Add support for step-size >1

# Establish DB connections
#r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
#session = db.Session()


skipped = 0
selected = 0
alreadyCompleted = 0
totalMissingResults = 0

queuedDelayedCalls = []

for taxIdForProcessing in species:
    print("Procesing %d sequences for tax-id %d (%s)..."
          % (countSpeciesCDS(taxIdForProcessing),
             taxIdForProcessing,
             getSpeciesName(taxIdForProcessing)))

    # Iterate over all CDS entries for this species
    # TODO - preloading all sequences and results should optimize this
    for protId in SpeciesCDSSource(taxIdForProcessing):
        #protId = codecs.decode(protId)
        # Filtering

        # Only process 1/N of the sequences, selected randomly (N=randomFraction)
        # (if randomFraction==1, all sequences will be processed)
        if( randint(1, randomFraction) != 1 ):
            skipped += 1
            continue

        # Skip sequences with partial CDS annotations
        #if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        #    skipped += 1
        #    continue

        #if( not r.exists(nativeCdsSeqIdKey % (taxIdForProcessing, protId)) ):
        #    skipped +=1
        #    continue

        cds = CDSHelper(taxIdForProcessing, protId)

        seqLength = cds.length()
        if seqLength is None:
            print("Warning: Could not find CDS length entry for taxid=%d, protid=%s" % (taxIdForProcessing, protId) )
            skipped += 1
            continue

        assert( seqLength > 3 )
        # Skip sequences with length <40nt (window width)
        if(seqLength < windowWidth + 1 ):
            print("short seq")
            skipped += 1
            continue

        requiredWindows = list(range(0, min(seqLength - windowWidth - 1, lastWindowStart), windowStep))
        requiredShuffles = [-1] # Check the native profile, regardless of the requested range
        requiredShuffles.extend(range(fromShuffle, toShuffle+1))

        existingResults = None
        try:
            missingResults = cds.checkCalculationResultWithWindows( seriesSourceNumber, requiredShuffles, requiredWindows )
        except IndexError as e:
            print("Missing sequences for %s, skipping..." % protId)
            skipped += 1
            continue

        #print(missingResults)

        shufflesWithMissingWindows = [requiredShuffles[n] for n,v in enumerate(missingResults) if v ] # get the indices (not shuffle-ids!) of existing shuffles with missing positions
        print("Existing shuffles with missing windows: %s" % shufflesWithMissingWindows)
        completelyMissingShuffles = [requiredShuffles[n] for n,v in enumerate(missingResults) if v is None]
        print("Missing shuffles: %s" % completelyMissingShuffles)

        # Submit all missing shuffles for processing in a single task

        ret = None
        if( completelyMissingShuffles ):
            try:
                #ret = storeNewShuffles(cds.getTaxId(), cds.getProtId(), completelyMissingShuffles)

                ret = scheduler.submit(storeNewShuffles, cds.getTaxId(), cds.getProtId(), completelyMissingShuffles)
                newIds = ret.result()
                print("Created new seqs:")
                print(zip(completelyMissingShuffles, newIds))

                # reload cds helper data
                del cds
                cds = CDSHelper(taxIdForProcessing, protId)
                print("(done with new seqs)")

            except Exception as e:
                print("Error creating new seqs")
                print(e)
                skipped += 1
                continue

            

        lastwin = requiredWindows[-1]
        if(shufflesWithMissingWindows):
            #cds.enqueueForProcessing(computationTag, shufflesWithMissingWindows, lastwin, windowStep)#

            shuffleIdsToProcess = sorted(shufflesWithMissingWindows + completelyMissingShuffles)
            allSeqIds = cds.shuffledSeqIds()

            def shuffleIdToSeqId(shuffleId):
                if shuffleId==-1:
                    return cds.seqId()
                else:
                    return allSeqIds[shuffleId]
                
            requiredSeqIds = list(map(shuffleIdToSeqId, shuffleIdsToProcess))

            queueItem = "%d/%s/%s/%s/%d/%d" % (cds.getTaxId(), cds.getProtId(), ",".join(map(str, requiredSeqIds)), ",".join(map(str, shufflesWithMissingWindows + completelyMissingShuffles)), lastwin, windowStep)
            print(queueItem)

            # To maximize node utilization, we will delay the main part of the calcualtion, the energy calculation, until after
            # we finished creating all necessary sequences (Otherwise, both types of calculations are interleaved and we may
            # be unable to generate enough work when it is needed).
            #
            # An even better alternative might be to interleave submission of both types of tasks, but give the "loading"
            # tasks higher priority.
            
            delayedCall = dask.delayed( calculateMissingWindowsForSequence )(taskDescription=queueItem) # create a delayed call for the calculations needed
            queuedDelayedCalls.append( delayedCall ) # store the call for later submission

            #print("%s: enqueued %d additional results" % (protId, len(shufflesWithMissingWindows)))
            totalMissingResults += len(shufflesWithMissingWindows)
        else:
            print("No pending shuffles, skipping...")
            skipped += 1
            continue
                

print("Added %d additional calculations" % totalMissingResults)
print("%d sequences selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
#print("queue contains %d items" % numItemsInQueue(computationTag))
#print("%d proteins queued" % len(queuedResults))

futures = scheduler.compute(queuedDelayedCalls) # submit all delayed calculations; obtain futures immediately

_distributed.progress(futures) # wait for all calculations to complete
print("\n")

results = scheduler.gather(futures) # get the results of the completed calcualtions
print(results)
# todo -- recover partial results when errors occur

del scheduler # free connections?
