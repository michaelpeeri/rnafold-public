# Enqueue failed entries (i.e. entries missing the calculation results) for repeated processing.
# Scan all CDS entries from a given taxid; Find entries missing a computation result; Re-insert them into the queue for repeated processing
# Input - taxid
#       - computationTag
# Example:
# python2 requeue_sequences_missing_energies_for_sliding_window.py 3055 rna-fold-window-40-0
# TODO: Add support for step-size >1
import sys
import codecs
import redis
import config

computationTag = sys.argv[2]
if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0

# Configuration
queueKey = "queue:tag:awaiting-%s:members" % computationTag
calculationResultKey = "CDS:taxid:%%d:protid:%%s:computed:%s:value" % computationTag
print(calculationResultKey)
windowWidth = 40
# TODO: Add support for step-size >1

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
    protId = codecs.decode(protId)
    # Filtering
    # Skip sequences with partial CDS annotations
    if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    if(not r.exists("CDS:taxid:%d:protid:%s:seq" % (taxIdForProcessing, protId))):
        print("Warning: couldn't find sequence for protid '%s'" % protId )
        skipped += 1
        continue

    seqLength = r.strlen("CDS:taxid:%d:protid:%s:seq" % (taxIdForProcessing, protId))
    # Skip sequences with length <40nt (window width)
    if(seqLength < windowWidth ):
        skipped += 1
        continue

    numWindows = seqLength - windowWidth + 1

    # Do the computed folding energies exist (and contain a value for each window?)
    if( r.exists( calculationResultKey % (taxIdForProcessing, protId)) and
        r.llen( calculationResultKey % (taxIdForProcessing, protId)) == numWindows ):
        skipped += 1
        continue

    # Enqueue this protein for processing
    r.rpush(queueKey, "%s:%s" % (taxIdForProcessing, protId))
    selected += 1


print("%d selected, %d skipped (%d total)" % (selected, skipped, selected+skipped))
print("queue contains %d items" % r.llen(queueKey))
