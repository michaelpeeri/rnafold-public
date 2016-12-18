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
#argsParser.add_argument("--gff")
#argsParser.add_argument("--variant", type=parseOption(set(("yeastgenome","NCBI","Ensembl")), "variant"))
#argsParser.add_argument("--output-gene-ids")
#argsParser.add_argument("--areas", type=parseList(int))
#argsParser.add_argument("--transl-table", type=int)
#argsParser.add_argument("--alt-protein-ids", type=parseOption(set(("locus_tag",)), "alt-protein-id"))
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
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
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
        if seqLength is None:
            print("Warning: Could not find CDS length entry for taxid=%d, protid=%s" % (taxIdForProcessing, protId) )
            skipped += 1
            continue

        assert( seqLength > 3 )
        # Skip sequences with length <40nt (window width)
        if(seqLength < windowWidth + 1 ):
            skipped += 1
            continue

        #requiredShuffles = [-1] + list(range(0, min(seqLength, lastWindowStart), windowStep))

        requiredShuffles = range(fromShuffle, toShuffle)


        #existingResults = cds.checkCalculationResult( seriesSourceNumber, range(-1, expectedNumberOfShuffles) )
        existingResults = None
        try:
            existingResults = cds.checkCalculationResult( seriesSourceNumber, requiredShuffles )
        except IndexError as e:
            print("Missing sequences for %s, skipping..." % protId)
            skipped += 1
            continue 
        #assert(len(existingResults) == expectedNumberOfShuffles + 1 )
        
        pendingShuffles = []
        for shuffleId, resultOk in enumerate(existingResults):
            if not resultOk:
                pendingShuffles.append( shuffleId-1 )

        if(pendingShuffles):
            lastwin = min( lastWindowStart, seqLength-windowWidth-1)
            lastwin -= lastwin%windowStep
            cds.enqueueForProcessing(computationTag, pendingShuffles, lastwin, windowStep)

            print("%s: enqueued %d additional results" % (protId, len(pendingShuffles)))
            totalMissingResults += len(pendingShuffles)
                
        #    selected += 1
        #else:
        #    alreadyCompleted += 1

print("Added %d additional calculations" % totalMissingResults)
print("%d sequences selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
print("queue contains %d items" % numItemsInQueue(computationTag))
