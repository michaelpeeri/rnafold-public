# Worker process for calculating RNA folding energy profiles, using a sliding window
# Read sequences queued for processing, do the processing for each, and store the result in the sequence's entry
import sys
import redis
import RNA
import config
import mysql_rnafold as db
from sqlalchemy import sql
from sqlalchemy.sql.expression import func

# Command-line arguments
#computationTag = sys.argv[1]
#if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0

# Configuration
#queueTag = "queue:tag:awaiting-%s:members" % computationTag
seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
windowWidth = 40
windowStep = 1

# Connections
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
session = db.Session()

#assert(r.exists(queueTag))


class XException(object):
    pass

a = ["3055:XP_001696195.1:8490"]

while(True):
    try:
        itemToProcess = None
        # Remove an item from the queue.
        # Note: Failed sequences (i.e. those removed from the queue but not process successfully) will be requeued using
        # the requeue_sequences scripts.
        itemToProcess = a.pop() # r.lpop(queueTag)
        if( itemToProcess is None):
            print("No more items to process...")
            break

        taxId, protId, seqId = itemToProcess.split(":")
        taxId = int(taxId)
        seqId = int(seqId)

        # Get the sequence for this entry
        seq = db.connection.execute( sql.select( (db.sequences.c.sequence,)).select_from(db.sequences).where(
                        db.sequences.c.id==seqId )
                        ).scalar()

        expectedSeqLength = r.get(seqLengthKey % (taxId, protId))
        if( not expectedSeqLength is None ):
            # Skip sequences with length <40nt (window width)
            expectedSeqLength = int(expectedSeqLength)
            if( expectedSeqLength != len(seq) ):
                print("taxid=%d, protid=%s, seqid=%d - unexpected length %d (expected: %d)" % (taxId, protId, seqId, len(seq), expectedSeqLength) )
                continue

        if( len(seq) < windowWidth ):
            # Sequence is shorter than required window; skip
            continue

        # Each entry in the queue is in the format "taxid:protid"
        print("Processing item %s (length=%d, %d windows)..." % (itemToProcess, len(seq), len(seq) - windowWidth + 1))

        # TODO - Remove any old value stored in this key?
        #prevWindowCount = db.connection.execute( sql.select(( sql.func.count('*'),)).select_from(db.sequence_series).where(
        #        sql.and_(
        #            db.sequence_series.c.sequence_id==seqId,
        #            db.sequence_series.c.source==seriesSourceNumber,
        #            db.sequence_series.c.ext_index==0
        #            )) ).scalar()
        #print("prevWindows = %d" % (prevWindowCount,))
        #if( prevWindowCount >= len(seq)-windowWidth+1 ):
        #    # Sufficient data seems to exist. Skip...
        #    continue

        for start in range(len(seq)-windowWidth+1):
            fragment = seq[start:(start+windowWidth)]
            assert(len(fragment)==windowWidth)

            # Calculate the RNA folding energy. This is the computation-heavy part.
            strct, energy = RNA.fold(fragment)
            assert(energy <= 0.0)
            print("(%d) %.3g" % (start, energy))

            # Store the calculation result
            #print("%d:%s --> %f" % (taxId, protId, energy))

            #s1 = db.SequenceSeries(sequence_id=seqId, source=seriesSourceNumber, ext_index=0, index=start, value=energy)
            #session.add(s1)

        #session.commit()

        
        # Verify the number of results we stored matches the number of windows
        #print(r.llen(computationResultTag % (taxId, protId)))
        #print(len(seq)-windowWidth+1)
        #assert(r.llen(computationResultTag % (taxId, protId)) == len(seq) - windowWidth + 1)


    except XException:
        print(str(Exception))
        break

print("Done!")
