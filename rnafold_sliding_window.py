# Create mRNA folding energy profiles for sequences, by calculation energy for windows missing from the profile
# (Read sequences and existing results, do the processing for missing values, and store the result in the sequence's entry)
# This module can be used as a stand-alone worker process for calculating RNA folding energy profiles, or expose this functionality for use elsewhere
import os
import subprocess
import gzip
import json
from time import sleep
from datetime import datetime, timedelta
import argparse
import mysql_rnafold as db
from data_helpers import CDSHelper, countSpeciesCDS, calcCrc, QueueSource, createWorkerKey, updateJobStatus_ItemCompleted
from runningstats import RunningStats
from rate_limit import RateLimit
from rnafold_vienna import RNAfold_direct
import logging


# hold configuration arguments (for either stand-alone or module use)
args = None


def parseOption(possibleValues, name):
    def checkOption(value):
        if value in possibleValues:
            return value
        else:
            raise argparse.ArgumentTypeError("Unknown %s '%s', allowed values: %s" % (name, value, ",".join(possibleValues)))
    return checkOption


"""
Check if a timer has expired (used to implement running-time limit)
"""
class RunUntil(object):
    def __init__(self, duration=timedelta(0, 0, 0, 0, 5) ):
        self._endTime = datetime.now() + duration
        print("Program will exit at %s" % self._endTime)

    def isExpired(self):
        return datetime.now() > self._endTime


def round4(x):
    if x is None:
        return None
    else:
        return round(x,4)



class RNAFold(object):
    def __init__(self, windowWidth=40, logfile=None, debugDoneWriteResults=False, computationTag="rna-fold-window-40-0", seriesSourceNumber=db.Sources.RNAfoldEnergy_SlidingWindow40_v2):
        if logfile is None:
            self._logfile = open('/dev/null', 'w')
        else:
            self._logfile = logfile

        self._debugDoneWriteResults = debugDoneWriteResults

        self._compuatationTag = computationTag
        self._seriesSourceNumber = seriesSourceNumber
        self._windowWidth = windowWidth


            
    def calculateMissingWindowsForSequence(self, taxId, protId, seqIds, requestedShuffleIds, firstWindow, lastWindowStart, windowStep, reference="begin"):

        logging.info("Parameters: %d %s %s %s %d %d %s" % (taxId, protId, seqIds, requestedShuffleIds, lastWindowStart, windowStep, reference))
        f = self._logfile

        assert(len(seqIds)>0)
        assert(len(seqIds)==len(requestedShuffleIds))

        if( reference != "begin"):
            raise Exception("Specificed profile reference '%s' is not supported!" % reference)

        # We will process all listed shuffle-ids for the following protein record
        cds = CDSHelper( taxId, protId )

        if( cds.length() < self._windowWidth ):
            e = "Refusing to process item %s because the sequence length (%d nt) is less than the window size (%d nt)\n" % (itemToProcess, cds.length(), self._windowWidth)
            f.write(e)
            logging.error(e)
            raise Exception(e)

        # Create a list of the windows we need to calculate for this CDS
        requestedWindowStarts = frozenset(range(0, min(lastWindowStart, cds.length()-self._windowWidth-1), windowStep ))
        if( len(requestedWindowStarts) == 0):
            e = "No windows exist for calculation taxid=%d, protId=%s, CDS-length=%d, lastWindowStart=%d, windowStep=%d, windowWidth=%d - Skipping...\n" % (taxId, protId, cds.length(), lastWindowStart, windowStep, self._windowWidth)
            f.write(e)
            logging.error(e)
            raise Exception(e)

        # First, read available results (for all shuffle-ids) in JSON format
        # Array is indexed by shuffle-id, so results not requested will be represented by None (as will requested items that have no results yet).
        logging.info("DEBUG: requestedShuffleIds (%d items): %s\n" % (len(requestedShuffleIds), requestedShuffleIds))
        existingResults = cds.getCalculationResult2( self._seriesSourceNumber, requestedShuffleIds, True )
        #assert(len(existingResults) >= len(requestedShuffleIds))  # The returned array must be at least as large as the requested ids list
        assert(len(existingResults) == max(requestedShuffleIds)+1)
        #existingResults = [None] * (max(requestedShuffleIds)+1)
        logging.info("DEBUG: existingResults (%d items): %s\n" % (len(existingResults), existingResults))

        # Check for which of the requested shuffle-ids there are values missing
        shuffleIdsToProcess = {}
        for n, r in enumerate(existingResults):
            shuffleIdForThisN = n-1
            if r is None:
                # There are no existing results for shuffled-id n. If it was requested, it should be calculated now (including all windows)
                if n in requestedShuffleIds:
                    shuffleIdsToProcess[shuffleIdForThisN] = list(requestedWindowStarts)
                    
                # ------------------------------------------------------------------------------------
                continue   # TODO - verify this line; should we abort this sequence by throwing????
                # ------------------------------------------------------------------------------------
            
            # Check the existing results for this shuffle
            alreadyProcessedWindowStarts = frozenset( [i for i,x in enumerate(r["MFE-profile"] ) if x is not None] ) # Get the indices (=window starts) of all non-None values
            missingWindows = requestedWindowStarts - alreadyProcessedWindowStarts # Are there any requested windows that are not already computed?
            if( missingWindows ): 
                shuffleIdsToProcess[shuffleIdForThisN] = missingWindows

        if( not shuffleIdsToProcess):
            e = "All requested shuffle-ids in (taxId: %d, protId: %s, seqs: %s) seems to have already been processed. Skipping...\n" % (taxId, protId, str(list(zip(seqIds, requestedShuffleIds))) )
            logging.error(e)
            raise Exception(e)
        logging.info("DEBUG: shuffleIdsToProcess (%d items): %s\n" % (len(shuffleIdsToProcess), shuffleIdsToProcess))

        logging.info("DEBUG: Before (%d items): %s\n" % (len(existingResults), existingResults))
        # Initialize new results records
        for shuffleId in shuffleIdsToProcess.keys():
            if existingResults[shuffleId+1] is None:
                logging.info(seqIds)
                logging.info(requestedShuffleIds)
                logging.info(shuffleId)
                thisSeqId = seqIds[ requestedShuffleIds.index(shuffleId) ]
                existingResults[shuffleId] = { "id": "%s/%s/%d/%d" % (taxId, protId, thisSeqId, shuffleId), "seq-crc": None, "MFE-profile": [], "MeanMFE": None, "v": 2 }
        logging.info("DEBUG: existingResults (%d items): %s\n" % (len(existingResults),existingResults) )

        # Load the sequences of all shuffle-ids we need to work on
        # TODO - combine loading of multiple sequences into one DB operation
        for shuffleId, record in [(u-1,v) for u,v in enumerate(existingResults) if not v is None]:
            seq = None
            annotatedSeqId = None
            # Get the sequence for this entry
            if( shuffleId < 0 ):
                seq = cds.sequence()
                annotatedSeqId = cds.seqId()
            else:
                seq = cds.getShuffledSeq(shuffleId)
                annotatedSeqId = cds.getShuffledSeqId(shuffleId)

            if( seq is None or (not seq is None and len(seq)==0 )):
                seq2 = cds.getShuffledSeq2( annotatedSeqId )
                seq3 = cds._fetchSequence( annotatedSeqId )
                seq4 = cds._cache.get("%d:seq"%annotatedSeqId)
                if not seq4 is None:
                    del cds._cache["%d:seq"%annotatedSeqId]
                seq5 = cds.getShuffledSeq2( annotatedSeqId )
                e = "Got empty sequence for shuffleId=%d, seqId=%d, taxId=%d, protId=%s, numShuffled=%d, ids[%d:%d]=%s, len(seq2)=%d, len(seq3)=%d, len(seq4)=%d, len(seq5)=%d" % (shuffleId, annotatedSeqId, taxId, protId, len(cds.shuffledSeqIds()), shuffleId-2, shuffleId+2, cds.shuffledSeqIds()[shuffleId-2:shuffleId+2], len(seq2) if not seq2 is None else -1, len(seq3) if not seq3 is None else -1, len(seq4) if not seq4 is None else -1, len(seq5) if not seq5 is None else -1 )
                logging.error(e)
                raise Exception(e)

            #
            # Disabled - calculation needn't include the native sequence...
            #
            #if( annotatedSeqId not in seqIds ):
            #    e = "Error: SeqId specified in queue item %s does not match annotated seq-id %d\n" % (itemToProcess, annotatedSeqId)
            #    f.write(e)
            #    f.write("Current shuffle-id: %d\n" % shuffleId)
            #    f.write("Ids in existing results:\n")
            #    for shuffleId, record in enumerate(existingResults):
            #        f.write(" %d) %s\n" % (shuffleId, record['id']))
            #    f.write("Debug info:\n")
            #    f.write("\n".join(cds.getDebugInfo()))
            #    f.write("\n")
            #    f.write("Skipping...\n")
            #    print("Skipping...")
            #    raise Exception(e)

            expectedSeqLength = cds.length()
            if( not expectedSeqLength is None ):
                if( expectedSeqLength != len(seq) ):
                    e = "Warning: taxid=%d, protid=%s, seqid=%d - unexpected length %d (expected: %d)\n" % (taxId, protId, annotatedSeqId, len(seq), expectedSeqLength)
                    f.write(e)
                    logging.error(e)
                    raise Exception(e)

            if( len(seq) < self._windowWidth ):
                # Sequence is shorter than required window; skip
                e = "Warning: skipping sequence because it is shorter than the requested window...\n"
                f.write(e)
                logging.error(e)
                raise Exception(e)

            logging.info("DEBUG: Processing item taxId=%d, protId=%s, shuffle=%d (length=%d, %d windows)...\n" % (taxId, protId, shuffleId, len(seq), len(requestedWindowStarts)))

            # TODO - Remove any old value stored in this key?

            # Skip this for now
            # This will be made redundant by completing the "updating" implementation
            #
            #if( cds.isCalculationDone( seriesSourceNumber, shuffleId )):
            #    # Sufficient data seems to exist. Skip...
            #    f.write("Item %s appears to be already completed, skipping..." % itemToProcess)
            #    continue

            logging.info(seq[:50])
            #f.write("\n")

            MFEprofile = record["MFE-profile"]
            #f.write("Profile: %s\n" % MFEprofile)

            # Make sure the profile array contains enough entries for all new windows (and possibly, if windows are non-contiguous, entries between them that we are not going to compute right now)
            if( len(MFEprofile) < max(requestedWindowStarts) ):
                entriesToAdd = max(requestedWindowStarts) - len(MFEprofile) + 1
                MFEprofile.extend( [None] * entriesToAdd )
            assert(len(MFEprofile) >= max(requestedWindowStarts))

            stats = RunningStats()
            stats.extend([x for x in MFEprofile if x is not None])

            for start in requestedWindowStarts:
                fragment = seq[start:(start+self._windowWidth)]
                assert(len(fragment)==self._windowWidth)

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

            if( not self._debugDoneWriteResults):
                cds.saveCalculationResult2( self._seriesSourceNumber, result, annotatedSeqId, False )

        
        if( not self._debugDoneWriteResults):
            cds.commitChanges()


def parseTaskDescription(taskDescription):
    parts = []
    if( taskDescription.find('/') == -1 ):
        parts = taskDescription.split(":")
    else:
        parts = taskDescription.split("/")

    # Mandatory fields
    taxId, protId, seqId_, shuffleIds_ = parts[:4]
    taxId = int(taxId)
    seqIds = map(int, seqId_.split(","))
    requestedShuffleIds = list(map(int, shuffleIds_.split(",")))
    assert(len(requestedShuffleIds) < 50 ) # Maximum number of shuffle-ids allowed to be combined into a single work record
    assert(seqIds)
    assert(len(seqIds)==len(requestedShuffleIds))

    # Optional fields
    lastWindowStart = 150
    windowStep = 1
    if( len(parts)>4 ):
        lastWindowStart = int(parts[4])
        #assert(lastWindowStart >= 150)
        windowStep = int(parts[5])
        assert(windowStep > 0 and windowStep <= 30)

    return (taxId, protId, seqIds, requestedShuffleIds, lastWindowStart, windowStep)



def calculateMissingWindowsForSequence(taskDescription):
    logging.info("Processing task %s" % taskDescription)
    #def __init__(self, windowWidth=40, logfile=None, debugDoneWriteResults=False, computationTag="rna-fold-window-40-0", seriesSourceNumber=db.Sour
    (taxId, protId, seqIds, requestedShuffleIds, lastWindowStart, windowStep) = parseTaskDescription(taskDescription)

    windowWidth = 40 # Todo: get this
    windowRef = "begin"
    firstWindowStart = 0
    
    f = RNAFold(windowWidth)
    #def calculateMissingWindowsForSequence(self, taxId, protId, seqIds, requestedShuffleIds, firstWindow, lastWindowStart, windowStep, reference="begin"):
    ret = f.calculateMissingWindowsForSequence(taxId, protId, seqIds, requestedShuffleIds, firstWindowStart, lastWindowStart, windowStep, windowRef)

    return taskDescription


"""
Main function for stand-alone worker process -- read and process queued items from redis queue
"""            
def standaloneMainWithRedisQueue():
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--computation-tag")
    argsParser.add_argument("--limit-time")
    argsParser.add_argument("--dry-run", type=parseOption(frozenset(("no", "yes")), "dry-run") )
    args = argsParser.parse_args()

    # Command-line arguments
    computationTag = args.computation_tag
    if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
    # e.g. rna-fold-0

    runUntil = None
    if( args.limit_time ):
        hours, mins = map(int, args.limit_time.split(':'))
        assert(mins<60)
        runUntil = RunUntil(timedelta(0, 0, 0, 0, mins, hours))


    # Configuration
    #queueTag = "queue:tag:awaiting-%s:members" % computationTag
    #seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
    seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40_v2
    windowWidth = 40
    windowStep = 1


    myWorkerKey = createWorkerKey(computationTag)

    updateJobStatus_ItemCompleted( myWorkerKey )

    
    with gzip.open('worker_%s_%d_log.gz' % (subprocess.check_output('hostname')[:-1], os.getpid()), 'wb') as f:
        print("My worker-key is: %s" % myWorkerKey )
        f.write("My worker-key is: %s\n" % myWorkerKey )

        rnafold = RNAFold(windowWidth, f, False, computationTag, seriesSourceNumber)

        for itemToProcess in QueueSource(computationTag):
            parts = []
            if( itemToProcess.find('/') == -1 ):
                parts = itemToProcess.split(":")
            else:
                parts = itemToProcess.split("/")

            # Mandatory fields
            taxId, protId, seqId_, shuffleIds_ = parts[:4]
            taxId = int(taxId)
            seqIds = frozenset(map(int, seqId_.split(",")))
            requestedShuffleIds = list(map(int, shuffleIds_.split(",")))
            assert(len(requestedShuffleIds) < 50 ) # Maximum number of shuffle-ids allowed to be combined into a single work record

            # Optional fields
            lastWindowStart = 150
            windowStep = 1
            if( len(parts)>4 ):
                lastWindowStart = int(parts[4])
                #assert(lastWindowStart >= 150)
                windowStep = int(parts[5])
                assert(windowStep > 0 and windowStep <= 30)

            try:
                rnafold.calculateMissingWindowsForSequence(taxId, protId, seqId_, shuffleIds_, 0, lastWindowStart, windowStep, "begin", f)
                updateJobStatus_ItemCompleted( myWorkerKey, taxId )
                
            except Exception as e:
                print(e)


            if( not runUntil is None):
                if( runUntil.isExpired() ):
                    f.write("Time limit reached, exiting...\n")
                    break

        f.write("Done!\n\n")
        
    return 0
    

    

if __name__=="__main__":
    import sys
    #sys.exit(standaloneMainWithRedisQueue())
    #testTask = "436017/ABO99737/7086688,20745703,7868392,20366322,7856723,7875007,21062723,21062724,21062725,21062726,21062727,21062728,21062729,21062730,21062731,21062732,21062733/-1,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160/780/10"
    #  shuffleId=150, seqId=21063218, taxId=436017, protId=ABO96972
    #testTask = "436017/ABO96972/21063218/150/780/10"
    #testTask = "436017/ABO94655/7087318,8008675,7880166,20829673,20600299,20676590,21064374,21064375,21064376,21064377,21064378,21064379,21064380,21064381,21064382,21064383,21064384/-1,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160/350/10"
    #testTask = "436017/ABP00033/7076680,20044248,20045786,7638494,20042209,20129253,21064286,21064287,21064288,21064289,21064290,21064291,21064292,21064293,21064294,21064295,21064296/-1,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160/640/10"
    testTask = "436017/ABP00428/7087197,20371424,20676675,20750822,20600300,7810032,21064946,21064947,21064948,21064949,21064950,21064951,21064952,21064953,21064954,21064955,21064956/-1,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160/370/10"
    calculateMissingWindowsForSequence(testTask)
    sys.exit(0)
else:
    # TODO - define args when used as a module
    args = {}
