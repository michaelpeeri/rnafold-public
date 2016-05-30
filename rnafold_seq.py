# Worker process for the CDS sequences queue
# Read sequences queued for processing, do the processing for each, and store the result in the sequence's entry
import redis
import RNA
import config

# Configuration
queueTag = "queue:tag:awaiting-rna-fold-0:members"
sequenceTag = "CDS:taxid:%d:protid:%s:seq"
computationResultTag = "CDS:taxid:%d:protid:%s:computed:rna-fold-0:energy"

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

        # Each entry in the queue is in the format "taxid:protid"
        print("Processing item %s..." % itemToProcess)
        taxId, protId = itemToProcess.split(":")
        taxId = int(taxId)

        # Get the sequence for this entry
        seq = r.get(sequenceTag % (taxId, protId))
        # Calculate the RNA folding energy. This is the computation-heavy part.
        strct, energy = RNA.fold(seq)

        # Store the calculation result
        print("%d:%s --> %f" % (taxId, protId, energy))
        r.set(computationResultTag % (taxId, protId), energy)


    except XException:
        print(str(Exception))
        break

print("Done!")
