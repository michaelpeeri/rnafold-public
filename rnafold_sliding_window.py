# Worker process for calculating RNA folding energy profiles, using a sliding window
# Read sequences queued for processing, do the processing for each, and store the result in the sequence's entry
import sys
import os
import subprocess
import re
import gzip
from time import sleep
from datetime import datetime, timedelta
import mysql_rnafold as db
from data_helpers import CDSHelper, RateLimit, countSpeciesCDS, calcCrc, QueueSource, createWorkerKey, updateJobStatus_ItemCompleted
from runningstats import RunningStats
 

# Command-line arguments
computationTag = sys.argv[1]
if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0

class RunUntil(object):
    def __init__(self, duration=timedelta(0, 0, 0, 0, 5) ):
        self._endTime = datetime.now() + duration
        print("Program will exit at %s" % self._endTime)

    def isExpired(self):
        return datetime.now() > self._endTime

runUntil = None
if( len(sys.argv) > 2 ):
    hours, mins = map(int, sys.argv[2].split(':'))
    assert(mins<60)
    runUntil = RunUntil(timedelta(0, 0, 0, 0, mins, hours))
    

# Configuration
#queueTag = "queue:tag:awaiting-%s:members" % computationTag
#seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
windowWidth = 40
windowStep = 1

# Connections
#r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
#session = db.Session()

#assert(r.exists(queueTag))



reMFEScore = re.compile(".*\n.*[(]\s*([\d.-]+)[)]\n")
def RNAfold_direct(seq):
    out = subprocess.check_output("echo %s | ~/anaconda2_node/bin/RNAfold --noPS" % seq, shell=True)
    score = float(reMFEScore.match(out).group(1))
    assert(score<=0.0)
    return score


myWorkerKey = createWorkerKey(computationTag)
print("My worker-key is: %s" % myWorkerKey )

updateJobStatus_ItemCompleted( myWorkerKey )


with gzip.open('worker_%s_%d' % (subprocess.check_output('hostname')[:-1], os.getpid()), 'wb') as f:
    for itemToProcess in QueueSource(computationTag):
        taxId, protId, seqId, shuffleId = itemToProcess.split(":")
        taxId = int(taxId)
        seqId = int(seqId)
        shuffleId = int(shuffleId)

        cds = CDSHelper( taxId, protId )

        seq = None
        annotatedSeqId = None
        # Get the sequence for this entry
        if( shuffleId < 0 ):
            seq = cds.sequence()
            annotatedSeqId = cds.seqId()
        else:
            seq = cds.getShuffledSeq(shuffleId)
            annotatedSeqId = cds.getShuffledSeqId(shuffleId)

        if( annotatedSeqId != seqId ):
            print("Error: SeqId specified in queue item %s does not match annotated seq-id %d" % (itemToProcess, annotatedSeqId))
            continue

        expectedSeqLength = cds.length()
        if( not expectedSeqLength is None ):
            if( expectedSeqLength != len(seq) ):
                print("taxid=%d, protid=%s, seqid=%d - unexpected length %d (expected: %d)" % (taxId, protId, seqId, len(seq), expectedSeqLength) )
                continue

        if( len(seq) < windowWidth ):
            # Sequence is shorter than required window; skip
            continue

        # Each entry in the queue is in the format "taxid:protid:shuffleid"
        print("Processing item %s (length=%d, %d windows)..." % (itemToProcess, len(seq), len(seq) - windowWidth + 1))

        # TODO - Remove any old value stored in this key?

        if( cds.isCalculationDone( seriesSourceNumber, shuffleId )):
            # Sufficient data seems to exist. Skip...
            print("Item %s appears to be already completed, skipping..." % itemToProcess)
            continue

        print(seq[:50])

        MFEProfile = []
        stats = RunningStats()

        l = 0

        for start in range(len(seq)-windowWidth+1):
            fragment = seq[start:(start+windowWidth)]
            assert(len(fragment)==windowWidth)

            # Calculate the RNA folding energy. This is the computation-heavy part.
            #strct, energy = RNA.fold(fragment)
            energy = RNAfold_direct(fragment)
            assert(energy <= 0.0)

            # Store the calculation result
            #print("%d:%s --> %f" % (taxId, protId, energy))

            stats.push(energy)
            MFEProfile.append(energy)
            l += 1
            if( l>150 ):
                break

        # Format
        crc = calcCrc(seq)
        result = """{id="%s",seq-crc=%d,MFE-profile=[%s],MeanMFE=%.6g}""" % (itemToProcess, crc, ",".join(map(lambda x: "%.3g" % x, MFEProfile)), stats.mean())
        print(result[:100])
        print(len(result))

        cds.saveCalculationResult2( seriesSourceNumber, result, seqId )

        updateJobStatus_ItemCompleted( myWorkerKey, taxId )

        if( not runUntil is None):
            if( runUntil.isExpired() ):
                print("Time limit reached, exiting...")
                break
        
        
        # Verify the number of results we stored matches the number of windows
        #print(r.llen(computationResultTag % (taxId, protId)))
        #print(len(seq)-windowWidth+1)
        #assert(r.llen(computationResultTag % (taxId, protId)) == len(seq) - windowWidth + 1)


print("Done!")
