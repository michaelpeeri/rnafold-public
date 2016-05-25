import redis
import RNA
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

queueTag = "queue:tag:awaiting-rna-fold-for-shuffled-0:members"

class XException(object):
    pass

while(True):
    try:
        itemToProcess = None
        # Remove an item from the queue, and mark it as "at-risk"
        itemToProcess = r.lpop(queueTag)
        if( itemToProcess is None):
            print("No more items to process...")
            break

        print("Processing item %s..." % itemToProcess)
        taxId, protId = itemToProcess.split(":")
        taxId = int(taxId)

        seq = r.get("CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq" % (taxId, protId))
        strct, energy = RNA.fold(seq)
        print("%d:%s --> %f" % (taxId, protId, energy))
        r.set("CDS:taxid:%d:protid:%s:computed:rna-fold-for-shuffled-0:energy" % (taxId, protId), energy)

        #r.srem(atRiskTag, itemToProcess)

    except XException:
        print(str(Exception))
        break


print("Done!")
