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
import redis
import config
import mysql_rnafold as db
from sqlalchemy import sql
from sqlalchemy.sql.expression import func

# command-line arguments
species = map(int, sys.argv[1].split(","))
computationTag = sys.argv[2]
if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0
randomFraction = int(sys.argv[3])

# Configuration
queueKey = "queue:tag:awaiting-%s:members" % computationTag
shuffledSeqIdsKey = "CDS:taxid:%d:protid:%s:shuffled-seq-ids-v2"
seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
windowWidth = 40
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
# TODO: Add support for step-size >1

# Establish DB connections
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
session = db.Session()

skipped = 0
selected = 0
alreadyCompleted = 0

for taxIdForProcessing in species:
    print("Procesing %d sequences for tax-id %d (%s)..."
          % (r.scard("species:taxid:%d:CDS" % taxIdForProcessing),
             taxIdForProcessing,
             r.get("species:taxid:%d:name" % taxIdForProcessing)))

    # Iterate over all CDS entries for this species
    for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
        protId = codecs.decode(protId)
        # Filtering

        # Only process 1/N of the sequences, selected randomly (N=randomFraction)
        # (if randomFraction==1, all sequences will be processed)
        if( randint(1, randomFraction) != 1 ):
            skipped += 1
            continue

        # Skip sequences with partial CDS annotations
        if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
            skipped += 1
            continue

        if( not r.exists(shuffledSeqIdsKey % (taxIdForProcessing, protId)) ):
            skipped +=1
            continue

        seqLength = r.get(seqLengthKey % (taxIdForProcessing, protId))
        if( not seqLength is None ):
            # Skip sequences with length <40nt (window width)
            seqLength = int(seqLength)
            if(seqLength < windowWidth ):
                skipped += 1
                continue
        else:
            print("Warning: Could not find CDS length entry for taxid=%d, protid=%s" % (taxIdForProcessing, protId) )
            skipped += 1
            continue

        requiredNumWindows = seqLength - windowWidth + 1

        shuffledSeqIds = map(int, r.lrange(shuffledSeqIdsKey % (taxIdForProcessing, protId), 0, -1))

        anySeriesSelected = False
        print(shuffledSeqIds)
        for idx, currentSeqId in enumerate(shuffledSeqIds):
            # Count the number of window results stored for this series (in the current sequence)
            windowCount = db.connection.execute( sql.select(( sql.func.count('*'),)).select_from(db.sequence_series).where(
                    sql.and_(
                        db.sequence_series.c.sequence_id==currentSeqId,
                        db.sequence_series.c.source==seriesSourceNumber,
                        db.sequence_series.c.ext_index==0
                        )) ).scalar()
            # 
            if( windowCount < requiredNumWindows ):
                anySeriesSelected = True
                # enqueue this sequence for processing
                # Note: all windows will be computed (even if partial results exist)
                #print("%d:%s:%d" % (taxIdForProcessing, protId, currentSeqId))
                r.rpush(queueKey, "%d:%s:%d" % (taxIdForProcessing, protId, currentSeqId))

        if( anySeriesSelected ):  
            selected += 1
        else:
            alreadyCompleted += 1


print("%d selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
print("queue contains %d items" % r.llen(queueKey))
