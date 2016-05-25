import redis
import RNA
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

queueTag = "queue:tag:awaiting-rna-fold-0:members"
atRiskTag = "queue:tag:awaiting-rna-fold-0:risk"

class XException(object):
    pass

while(True):
    try:
        itemToProcess = None
        # Remove an item from the queue, and mark it as "at-risk"
        with r.pipeline() as pipe:
            pipe.multi()
            pipe.lpop(queueTag)
            pipe.sadd(atRiskTag, itemToProcess)
            itemToProcess, _ = pipe.execute()
            if( itemToProcess is None):
                print("No more items to process...")
                break

        print("Processing item %s..." % itemToProcess)
        taxId, protId = itemToProcess.split(":")
        taxId = int(taxId)

        seq = r.get("CDS:taxid:%d:protid:%s:seq" % (taxId, protId))
        strct, energy = RNA.fold(seq)
        print("%d:%s --> %f" % (taxId, protId, energy))
        r.set("CDS:taxid:%d:protid:%s:computed:rna-fold-0:energy" % (taxId, protId), energy)

        r.srem(atRiskTag, itemToProcess)

    except XException:
        print(str(Exception))
        break


atRiskCount = r.scard(atRiskTag)
if( atRiskCount > 0 ):
    print("WARNING: at-risk queue contains %d items!" % atRiskCount)

print("Done!")
