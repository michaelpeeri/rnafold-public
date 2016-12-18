# Worker process for calculating RNA folding energy profiles, using a sliding window
# Read sequences queued for processing, do the processing for each, and store the result in the sequence's entry
import sys
import os
import subprocess
import re
import gzip
import json
from time import sleep
from datetime import datetime, timedelta
import mysql_rnafold as db
from data_helpers import CDSHelper, countSpeciesCDS, calcCrc, QueueSource, createWorkerKey, updateJobStatus_ItemCompleted
from runningstats import RunningStats
from rate_limit import RateLimit 

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

def round4(x):
    if x is None:
        return None
    else:
        return round(x,4)

with gzip.open('worker_%s_%d_log.gz' % (subprocess.check_output('hostname')[:-1], os.getpid()), 'wb') as f:
    for itemToProcess in QueueSource(computationTag):
        parts = []
        if( itemToProcess.find('/') == -1 ):
            parts = itemToProcess.split(":")
        else:
            parts = itemToProcess.split("/")

        # Mandatory fields
        taxId, protId, seqId_, shuffleIds_ = parts[:4]
        taxId = int(taxId)
        seqIds = set(map(int, seqId_.split(",")))
        requestedShuffleIds = map(int, shuffleIds_.split(","))
        assert(len(requestedShuffleIds) < 50 ) # Maximum number of shuffle-ids allowed to be combined into a single work record

        # Optional fields
        lastWindowStart = 150
        windowStep = 1
        if( len(parts)>4 ):
            lastWindowStart = int(parts[4])
            #assert(lastWindowStart >= 150)
            windowStep = int(parts[5])
            assert(windowStep > 0 and windowStep <= 30)

        # We will process all listed shuffle-ids for the following protein record
        cds = CDSHelper( taxId, protId )

        if( cds.length() < windowWidth ):
            f.write("Refusing to process item %s because the sequence length (%d nt) is less than the window size (%d nt)" % (itemToProcess, cds.length(), windowWidth))
            continue

        # Create a list of the windows we need to calculate for this CDS
        requestedWindowStarts = list(range(0, min(lastWindowStart, cds.length()-windowWidth-1), windowStep ))
        assert(len(requestedWindowStarts) > 0)

        # First, read available results (for all shuffle-ids) in JSON format
        # Array is indexed by shuffle-id, so results not requested will be represented by None
        #existingResults = cds.getCalculationResult2( seriesSourceNumber, requestedShuffleIds, True )
        #assert(len(existingResults) >= len(requestedShuffleIds))
        #existingResults = [None] * (max(requestedShuffleIds)+1)
        #print(existingResults)

        # Check for which of the required shuffle-ids there are values missing
        shuffleIdsToProcess = {}

        # (TODO - complete and test this)
        #
        #for n, r in enumerate(existingResults):
        #    alreadyProcessedWindowStarts = set( [i for i,x in enumerate(r["MFE-profile"] ) if x is not None] ) # Get the indices (=window starts) of all non-None values
        #    missingWindows = requestedWindowStarts - alreadyProcessedWindows # Are there any requested windows that are not already computed?
        #    if( missingWindows ): 
        #        shuffleIdsToProcess[n-1] = missingWindows
        #if( not shuffleIdsToProcess):
        #    f.write("All requested shuffle-ids in %s seems to have already been processed. Skipping..." % (itemToProcess,))
        #    continue
        #f.write("DEBUG: shuffleIdsToProcess: %s" % shuffleIdsToProcess)

        # (For now, assume all records are missing...)
        #
        # Initialize records for missing results
        for shuffleId in requestedShuffleIds:
            shuffleIdsToProcess[shuffleId] = { "id": "%s/%s/%d" % (taxId, protId, shuffleId), "seq-crc": None, "MFE-profile": [], "MeanMFE": None, "v": 2 }

        print(shuffleIdsToProcess)
        
        # Load the sequences of all shuffle-ids we need to work on
        # TODO - combine loading of multiple sequences into one DB operation

        for shuffleId, record in shuffleIdsToProcess.iteritems():
            seq = None
            annotatedSeqId = None
            # Get the sequence for this entry
            if( shuffleId < 0 ):
                seq = cds.sequence()
                annotatedSeqId = cds.seqId()
            else:
                seq = cds.getShuffledSeq(shuffleId)
                annotatedSeqId = cds.getShuffledSeqId(shuffleId)

            if( annotatedSeqId not in seqIds ):
                f.write("Error: SeqId specified in queue item %s does not match annotated seq-id %d" % (itemToProcess, annotatedSeqId))
                continue

            expectedSeqLength = cds.length()
            if( not expectedSeqLength is None ):
                if( expectedSeqLength != len(seq) ):
                    f.write("taxid=%d, protid=%s, seqid=%d - unexpected length %d (expected: %d)" % (taxId, protId, seqId, len(seq), expectedSeqLength) )
                    continue

            if( len(seq) < windowWidth ):
                # Sequence is shorter than required window; skip
                f.write("Warning: skipping sequence because it is shorter than the requested window...")
                continue
            
            f.write("Processing item %s, shuffle=%d (length=%d, %d windows)..." % (itemToProcess, shuffleId, len(seq), len(requestedWindowStarts)))

            # TODO - Remove any old value stored in this key?

            # Skip this for now
            # This will be made redundant by completing the "updating" implementation
            #
            #if( cds.isCalculationDone( seriesSourceNumber, shuffleId )):
            #    # Sufficient data seems to exist. Skip...
            #    f.write("Item %s appears to be already completed, skipping..." % itemToProcess)
            #    continue

            f.write(seq[:50])

            MFEprofile = record["MFE-profile"]

            # Make sure the profile array contains enough entries for all new windows (and possibly, if windows are non-contiguous, entries between them that we are not going to compute right now)
            if( len(MFEprofile) < requestedWindowStarts[-1] ):
                entriesToAdd = requestedWindowStarts[-1] - len(MFEprofile) + 1
                MFEprofile.extend( [None] * entriesToAdd )
            assert(len(MFEprofile) >= requestedWindowStarts[-1] )
                
            stats = RunningStats()
            stats.extend([x for x in MFEprofile if x is not None])

            for start in requestedWindowStarts:
                fragment = seq[start:(start+windowWidth)]
                assert(len(fragment)==windowWidth)

                # Calculate the RNA folding energy. This is the computation-heavy part.
                #strct, energy = RNA.fold(fragment)
                energy = RNAfold_direct(fragment)
                assert(energy <= 0.0)

                # Store the calculation result
                #print("%d:%s --> %f" % (taxId, protId, energy))

                stats.push(energy)
                MFEprofile[start] = energy

            # Format
            crc = calcCrc(seq)
            #result = """{"id":"%s","seq-crc":%d,"MFE-profile":[%s],"MeanMFE":%.6g,v:2}""" % (itemToProcess, crc, ",".join(map(lambda x: "%.3g" % x, MFEprofile)), stats.mean())
            record["seq-crc"] = crc
            record["MFE-profile"] = [round4(x) for x in MFEprofile] # Round items down to save space (these are not exact numbers anyway)
            record["MeanMFE"] = stats.mean()
            result = json.dumps(record)
            
            f.write(result)
            f.write("\n")

            cds.saveCalculationResult2( seriesSourceNumber, result, annotatedSeqId )

        updateJobStatus_ItemCompleted( myWorkerKey, taxId )

        if( not runUntil is None):
            if( runUntil.isExpired() ):
                f.write("Time limit reached, exiting...")
                break

    f.write("Done!")
