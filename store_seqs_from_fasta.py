# Command-line args: <taxId> <fastaFile> <shuffledSeriesIndex>
# Read a fasta file containing the shuffled versions of the CDS for a given species;
# Store the shuffled sequences in MySql, and add the sequence-ids to the metadata in redis.
#
from __future__ import print_function
import sys
import redis
from Bio import SeqIO
import config
import mysql_rnafold as db
from sqlalchemy import sql
from sqlalchemy.sql.expression import func

# command-line arguments
taxId = int(sys.argv[1])
f = open(sys.argv[2], 'r')
shuffledSeriesIndex = int(sys.argv[3]) # 0-based

# configuration


# establish connections
# metadata server (redis)
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
# sequences server (mysql)
session = db.Session()

# redis keys
seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
cdsSeqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
shuffledSeqIdsKey = "CDS:taxid:%d:protid:%s:shuffled-seq-ids-v2"
partialCDSKey = "CDS:taxid:%d:protid:%s:partial"

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

    # Skip sequences that don't have an existing entry
    if( not r.exists(cdsSeqIdKey % (taxId, proteinId))):
        notFoundCount += 1
        print("Warning: Couldn't find existing entry for protein %s" % proteinId)
        continue

    # Skip partial sequences
    if(r.exists(partialCDSKey % (taxId, proteinId))):
        skippedCount += 1
        continue

    # Skip entries that have a length different than the shuffled length
    cdsLength = r.get(seqLengthKey % (taxId, proteinId))
    if( not cdsLength is None ):
        cdsLength = int(cdsLength)
        if( len(record.seq) != cdsLength ):
            skippedCount += 1
            print("Warning: Found entry for protein %s, but the original CDS length (%d) is different than the shuffled CDS length (%d)" % (proteinId, len(record.seq), cdsLength))
            continue
    else:
        print("Warning: Could not find CDS length entry for protein %s" % proteinId )
        # attempt go on with processing...
        

    # If the existing list for this protein
    nextSeriesIndex = r.llen(shuffledSeqIdsKey % (taxId, proteinId)) # Note - this will work if the list doesn't exist (e.g. when inserting the first items...)
    if( nextSeriesIndex < shuffledSeriesIndex ):
        print("ERROR: list for (taxid=%d, protId=%s) is missing values for previous series (shuffledSeriesIndex=%d, llen=%d" % (taxId, proteinId, shuffledSeriesIndex, nextSeriesIndex))
        print("Aborting...")
        break

    # Store the shuffled CDS sequence
    s1 = db.Sequence(sequence=record.seq, taxid=taxId, alphabet=db.Alphabets.RNA, source=db.Sources.ShuffleCDSv2)
    session.add(s1)
    session.commit()
    
    if( nextSeriesIndex == shuffledSeriesIndex ):
        r.rpush( shuffledSeqIdsKey % (taxId, proteinId), s1.id )
    else:
        r.lset( shuffledSeqIdsKey % (taxId, proteinId), shuffledSeriesIndex, s1.id )

    cdsCount += 1

if( notFoundCount + skippedCount > 0):
    print("Warning: %d entries skipped and %d entries not found" % (skippedCount, notFoundCount))

print("Processed %d CDS entries" % (cdsCount,))
print("(out of %d CDS entries for this species)" % (r.scard("species:taxid:%d:CDS" % (taxId,))))
