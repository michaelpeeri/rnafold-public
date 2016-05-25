# Input - taxid
# Enqueue failed entries for repeated processing.
# Scan all CDS entries from a given taxid; Find entries missing a computation result; Re-insert them into the queue for repeated processing
import sys
import redis
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

taxIdForProcessing = int(sys.argv[1])
print("Procesing %d sequences for tax-id %d (%s)..."
    % (r.scard("species:taxid:%d:CDS" % taxIdForProcessing),
    taxIdForProcessing,
    r.get("species:taxid:%d:name" % taxIdForProcessing)))

skipped = 0
selected = 0
queueKey = "queue:tag:awaiting-rna-fold-for-shuffled-0:members"
for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
    if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    # Is the computed folding energy missing (or not a non-positive float?)
    existingEnergy = r.get("CDS:taxid:%d:protid:%s:computed:rna-fold-for-shuffled-0:energy" % (taxIdForProcessing, protId))
    if(existingEnergy is None or not( float(existingEnergy) <= 0.0) ):
        selected += 1
        r.rpush(queueKey, "%s:%s" % (taxIdForProcessing, protId))

    skipped += 1

print("%d selected, %d skipped (%d total)" % (selected, skipped, selected+skipped))
print("queue contains %d items" % r.llen(queueKey))
