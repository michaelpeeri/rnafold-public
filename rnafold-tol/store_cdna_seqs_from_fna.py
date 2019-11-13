import sys
import redis
from Bio import SeqIO
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)


#f = open('/home/michael/rnafold/Creinherdtii.json', 'r')
#f = open('/home/michael/rnafold/Ptricornutum.json', 'r')
taxId = int(sys.argv[1])
f = open(sys.argv[2], 'r')

visitedProteinIds = set()


# Store a species entry for this species
#r.set('species:taxid:%d:name' % (taxId,), species)
#r.set('species:name:%s:taxid' % (species,), taxId)

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
    if( len(record.seq) < cdsLength ): # The cDNA must be at least as long as the CDS
        skippedCount += 1
        print("Warning: Found entry for protein %d, but the cDNA length (%d) is shorter than the CDS length (%d)" % (proteinId, len(record.seq), cdsLength))
        continue


    # Store the cDNA sequence
    r.set('CDS:taxid:%d:protid:%s:cdna-seq' % (taxId, proteinId), record.seq)
    cdsCount += 1

if( notFoundCount + skippedCount > 0):
    print("Warning: %d entries skipped and %d entries not found" % (skippedCount, notFoundCount))

print("Processed %d CDS entries" % (cdsCount,))
print("(out of %d CDS entries for this specied)" % (r.scard("species:taxid:%d:CDS" % (taxId,))))
