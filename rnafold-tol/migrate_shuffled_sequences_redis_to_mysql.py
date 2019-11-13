# Transfer main sequences from redis to MySQL
# Input: 
#
from __future__ import print_function
import sys
import redis
import config
import mysql_rnafold as db

# Configuration
species = (3055, 556484)

# Connect to redis (source DB)
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
# Connect to MySQL (destination DB)
session = db.Session()


skipped = 0
selected = 0

for taxIdForProcessing in species:
    # Iterate over all CDS entries for this taxId
    for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
        # Filters
        # (no filters defined)

        # Extract the CDS sequence
        cds = r.get("CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq" % (taxIdForProcessing, protId))
        if( cds is None):
            continue
        if( len(cds) == 0 ):
            continue
            

        selected += 1
        # Store the CDS sequence in destination db; get id of new record
        s1 = db.Sequence(sequence=cds, taxid=taxIdForProcessing, alphabet=db.Alphabets.RNA, source=db.Sources.Computed)
        session.add(s1)
        session.commit()

        # Store record-id in redis
        r.set("CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq-id" % (taxIdForProcessing, protId), s1.id)
        # Remove the sequence itself from redis
        r.delete("CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq" % (taxIdForProcessing, protId))

print("Moved %d sequences." % selected)
