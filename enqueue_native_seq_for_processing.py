# Input: taxid
# Enqueue all complete CDSs for a particular organism into a processing queue
import sys
import redis
import config

# Configuration
# Key name for this queue
queueKey = "queue:tag:awaiting-rna-fold-0:members"

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

taxIdForProcessing = int(sys.argv[1])
print("Procesing %d sequences for tax-id %d (%s)..."
    % (r.scard("species:taxid:%d:CDS" % taxIdForProcessing),
    taxIdForProcessing,
    r.get("species:taxid:%d:name" % taxIdForProcessing)))

skipped = 0
selected = 0
# Iterate over all CDS entries for this species
for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
    # Filters
    # Skip partial CDSs
    if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    selected += 1
    # Insert entry into the processing queue
    r.rpush(queueKey, "%s:%s" % (taxIdForProcessing, protId))

print("%d selected, %d skipped (%d total)" % (selected, skipped, selected+skipped))
print("queue contains %d items" % r.llen(queueKey))
