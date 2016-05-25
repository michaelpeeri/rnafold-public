from __future__ import print_function
import sys
import redis
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import NucleotideAlphabet
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)


taxIdForProcessing = int(sys.argv[1])


skipped = 0
selected = 0
for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
    # Filters
    # Skip partial CDSs
    if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    selected += 1
    # Insert into the processing queue
    cds = r.get("CDS:taxid:%d:protid:%s:seq" % (taxIdForProcessing, protId))
    record = SeqRecord(Seq(cds, NucleotideAlphabet), id=protId)
    print(record.format("fasta"), end='')
