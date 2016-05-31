# Worker process for calculating RNA folding energy profiles, using a sliding window
# Read sequences queued for processing, do the processing for each, and store the result in the sequence's entry
import sys
import redis
import RNA
import config

computationTag = sys.argv[1]
if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0

# Configuration
queueTag = "queue:tag:awaiting-%s:members" % computationTag
sequenceTag = "CDS:taxid:%d:protid:%s:seq" # TODO - think how to configure this...
computationResultTag = "CDS:taxid:%%d:protid:%%s:computed:%s:value" % computationTag
print("\"%s\"" % computationResultTag)
windowWidth = 40
windowStep = 1

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

class XException(object):
    pass

while(True):
    try:
        itemToProcess = None
        # Remove an item from the queue.
        # Note: Failed sequences (i.e. those removed from the queue but not process successfully) will be requeued using
        # the requeue_sequences scripts.
        itemToProcess = r.lpop(queueTag)
        if( itemToProcess is None):
            print("No more items to process...")
            break

        taxId, protId = itemToProcess.split(":")
        taxId = int(taxId)

        # Get the sequence for this entry
        seq = r.get(sequenceTag % (taxId, protId))

        if( len(seq) < windowWidth ):
            # Sequence is shorter than required window; skip
            continue

        # Each entry in the queue is in the format "taxid:protid"
        print("Processing item %s (length=%d, %d windows)..." % (itemToProcess, len(seq), len(seq) - windowWidth + 1))

        # Remove any old value stored in this key
        r.delete(computationResultTag % (taxId, protId))

        for start in range(len(seq)-windowWidth+1):
            fragment = seq[start:(start+windowWidth)]
            assert(len(fragment)==windowWidth)

            # Calculate the RNA folding energy. This is the computation-heavy part.
            strct, energy = RNA.fold(fragment)
            assert(energy <= 0.0)

            # Store the calculation result
            #print("%d:%s --> %f" % (taxId, protId, energy))
            r.lpush(computationResultTag % (taxId, protId), energy)
            #print("(%d,%d) -> %f" % (start, start+windowWidth, energy))
        
        # Verify the number of results we stored matches the number of windows
        #print(r.llen(computationResultTag % (taxId, protId)))
        #print(len(seq)-windowWidth+1)
        assert(r.llen(computationResultTag % (taxId, protId)) == len(seq) - windowWidth + 1)


    except XException:
        print(str(Exception))
        break

print("Done!")
