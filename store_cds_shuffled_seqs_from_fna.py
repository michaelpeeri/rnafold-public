import sys
import redis
from Bio import SeqIO
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)


taxId = int(sys.argv[1])
f = open(sys.argv[2], 'r')

visitedProteinIds = set()


assert(r.exists("species:taxid:%d:name" % taxId))

cdsCount = 0
notFoundCount = 0
skippedCount = 0
for record in SeqIO.parse(f, "fasta"):
    proteinId = record.id

    # verify there are no duplicates entries
    assert(proteinId not in visitedProteinIds)
    visitedProteinIds.add(proteinId)

    if( not r.exists("CDS:taxid:%d:protid:%s:seq" % (taxId, proteinId))):
        notFoundCount += 1
        print("Warning: Couldn't find existing entry for protein %s" % proteinId)
        continue

    cdsLength = r.strlen("CDS:taxid:%d:protid:%s:seq" % (taxId, proteinId))
    if( len(record.seq) != cdsLength ):
        skippedCount += 1
        print("Warning: Found entry for protein %d, but the original CDS length (%d) is different than the shuffled CDS length (%d)" % (proteinId, len(record.seq), cdsLength))
        continue


    # Store the shuffled CDS sequence
    r.set('CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq' % (taxId, proteinId), record.seq)
    cdsCount += 1

if( notFoundCount + skippedCount > 0):
    print("Warning: %d entries skipped and %d entries not found" % (skippedCount, notFoundCount))

print("Processed %d CDS entries" % (cdsCount,))
print("(out of %d CDS entries for this specied)" % (r.scard("species:taxid:%d:CDS" % (taxId,))))
