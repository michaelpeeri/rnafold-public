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
#import redis
import config
import mysql_rnafold as db
#from sqlalchemy import sql
#from sqlalchemy.sql.expression import func
from data_helpers import CDSHelper, countSpeciesCDS, getSpeciesName, SpeciesCDSSource, numItemsInQueue


# command-line arguments
species = map(int, sys.argv[1].split(","))
computationTag = sys.argv[2]
if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0
randomFraction = int(sys.argv[3])

# Configuration
#queueKey = "queue:tag:awaiting-%s:members" % computationTag
#nativeCdsSeqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
#seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
windowWidth = 40
numWindows = 150
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
expectedNumberOfShuffles = 50
# TODO: Add support for step-size >1

# Establish DB connections
#r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
#session = db.Session()

skipped = 0
selected = 0
alreadyCompleted = 0
totalMissingResults = 0

for taxIdForProcessing in species:
    print("Procesing %d sequences for tax-id %d (%s)..."
          % (countSpeciesCDS(taxIdForProcessing),
             taxIdForProcessing,
             getSpeciesName(taxIdForProcessing)))

    # Iterate over all CDS entries for this species
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
        if( not seqLength is None ):
            # Skip sequences with length <40nt (window width)
            if(seqLength < windowWidth + numWindows - 1 ):
                skipped += 1
                continue
        else:
            print("Warning: Could not find CDS length entry for taxid=%d, protid=%s" % (taxIdForProcessing, protId) )
            skipped += 1
            continue

        existingResults = cds.checkCalculationResult( seriesSourceNumber, range(-1, expectedNumberOfShuffles) )
        assert(len(existingResults) == expectedNumberOfShuffles + 1 )
        
        missingResults = 0

        for shuffleId, resultOk in enumerate(existingResults):
            if not resultOk:
                cds.enqueueForProcessing(computationTag, shuffleId)
                missingResults += 1

        if missingResults:
            print("%s: enqueued %d additional results" % (protId, missingResults))
            totalMissingResults += missingResults    
                
        #    selected += 1
        #else:
        #    alreadyCompleted += 1

print("Added %d additional calculations" % totalMissingResults)
print("%d sequences selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
print("queue contains %d items" % numItemsInQueue(computationTag))
